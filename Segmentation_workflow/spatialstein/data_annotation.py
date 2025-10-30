import numpy as np
from numpy import random as rd
from .data_preprocessing import get_msi_shape
from masserstein import Spectrum, estimate_proportions

def annotate_msi(msi, reference_spectra, MTD, MTD_th, control_ids=None, max_pixels=-1):
    """
    We will keep track of the total estimated proportion for each lipid,
    and the number of pixels in which this lipid had a non-zero proportion.
    This saves memory compared to keeping the deconvolved ion image for each lipid
    at this stage.
    However, to validate the estimation,
    we will save deconvolved ion images for selected ions.
    Note that since we estimate the proportions through
    linear spectral deconvolution,
    the correctness of estimation depends on the correctness of deconvolution.
    """
    
    proportions = np.zeros(len(reference_spectra))
    number_of_nonzero_pixels = np.zeros(len(reference_spectra), dtype='int')
    msi_shape = get_msi_shape(msi)
    if control_ids is not None:
        control_images = np.zeros(msi_shape + (len(control_ids), )) 

    # Check the m/z range of the reference spectra
    reference_mass_range = [
        min(s.confs[0][0] for s in reference_spectra) - 1,
        max(s.confs[-1][0] for s in reference_spectra) + 1
        ]
    
    # Total signal images for normalization purposes. 
    # global = whole spectra; local = restricted to the analyzed mass range.  
    # global_tic_array = np.zeros(msi_shape)
    # local_tic_array = np.zeros(msi_shape) 

    if max_pixels > -1 and max_pixels < len(msi_coordinates):
        selected_pixels = set(rd.choice(len(msi.coordinates),
                                    max_pixels,
                                    replace=False))
    else:
        selected_pixels = None
    analyzed_pixels = 0
    for idx, (xcoord,ycoord,zcoord) in enumerate(msi.coordinates):
        if selected_pixels is not None:
            if idx not in selected_pixels:
                continue
        else:
            if not idx % 10000:
                print('Processing pixel number', idx)
        analyzed_pixels += 1
        mz, intsy = msi.getspectrum(idx)
        
        ##  Get the total signal before truncation: 
        global_tic = np.sum(intsy)
        # global_tic_array[ycoord-1, xcoord-1] = global_tic
        
        # Truncate to the selected mass range to speed up computations:
        selected_range = (reference_mass_range[0] <= mz)*(mz <= reference_mass_range[1])
        mz = mz[selected_range]
        intsy = intsy[selected_range]
        
        ##  Get the total signal in the truncated region:
        local_tic = np.sum(intsy)
        # local_tic_array[ycoord-1, xcoord-1] = local_tic
        
        # Create the pixel spectrum object and normalize
        pixel_spectrum = Spectrum(confs=list(zip(mz, intsy)))
        pixel_spectrum.normalize()
        
        # Deconvolve the pixel spectrum
        regression = estimate_proportions(pixel_spectrum, reference_spectra, 
                                          MTD=MTD, MTD_th=MTD_th,
                                          MDC=1e-12,
                                          MMD=100*MTD, progress=False)
        pr_array = np.array(regression['proportions'])
        
        # Normalize the proportions so that they are equal the lipid's proportion
        # in the whole spectrum, not just in the analyzed mass range.
        # If the true lipid signal is S, then its proportion in the whole spectrum 
        # is S/G, which is equal to (S/L)*(L/G), with S/L estimated by masserstein
        pr_array = pr_array * local_tic / global_tic
        proportions += pr_array
        number_of_nonzero_pixels += pr_array > 0
        control_images[ycoord-1, xcoord-1, ...] = pr_array[control_ids]
    proportions /= analyzed_pixels
    return(proportions, number_of_nonzero_pixels, control_images, analyzed_pixels)
    
### Parallel implementation
def _annotate_pixel(pixel_id,
                    pixel_tuple,
                    reference_mass_range, reference_spectra,
                    MTD, MTD_th):
    ## Get the pixel
    mz, intsy = pixel_tuple

    ##  Get the total signal before truncation: 
    global_tic = np.sum(intsy)
    # global_tic_array[ycoord-1, xcoord-1] = global_tic
    
    # Truncate to the selected mass range to speed up computations:
    selected_range = (reference_mass_range[0] <= mz)*(mz <= reference_mass_range[1])
    mz = mz[selected_range]
    intsy = intsy[selected_range]
    
    ##  Get the total signal in the truncated region:
    local_tic = np.sum(intsy)
    # local_tic_array[ycoord-1, xcoord-1] = local_tic
    
    # Create the pixel spectrum object and normalize
    pixel_spectrum = Spectrum(confs=list(zip(mz, intsy)))
    pixel_spectrum.normalize()
    
    # Deconvolve the pixel spectrum
    regression = estimate_proportions(pixel_spectrum, reference_spectra, 
                                      MTD=MTD, MTD_th=MTD_th,
                                      MDC=1e-12,
                                      MMD=100*MTD, progress=False)
    pr_array = np.array(regression['proportions'])
    
    # Normalize the proportions so that they are equal the lipid's proportion
    # in the whole spectrum, not just in the analyzed mass range.
    # If the true lipid signal is S, then its proportion in the whole spectrum 
    # is S/G, which is equal to (S/L)*(L/G), with S/L estimated by masserstein
    pr_array = pr_array * local_tic / global_tic
    return (pixel_id, pr_array)

def annotate_msi_parallel(msi, reference_spectra, MTD, MTD_th,
                          control_ids=None, max_pixels=-1,
                          max_processes = -1):
    """
    We will keep track of the total estimated proportion for each lipid,
    and the number of pixels in which this lipid had a non-zero proportion.
    This saves memory compared to keeping the deconvolved ion image for each lipid
    at this stage.
    However, to validate the estimation,
    we will save deconvolved ion images for selected ions.
    Note that since we estimate the proportions through
    linear spectral deconvolution,
    the correctness of estimation depends on the correctness of deconvolution.
    """
    from joblib import Parallel, delayed, parallel_config, cpu_count
    if max_processes == -1:
        max_processes = cpu_count()
        
    proportions = np.zeros(len(reference_spectra))
    number_of_nonzero_pixels = np.zeros(len(reference_spectra), dtype='int')
    msi_shape = get_msi_shape(msi)
    if control_ids is not None:
        control_images = np.zeros(msi_shape + (len(control_ids), )) 

    # Check the m/z range of the reference spectra
    reference_mass_range = [
        min(s.confs[0][0] for s in reference_spectra),
        max(s.confs[-1][0] for s in reference_spectra)
        ]
    
    # Total signal images for normalization purposes. 
    # global = whole spectra; local = restricted to the analyzed mass range.  
    # global_tic_array = np.zeros(msi_shape)
    # local_tic_array = np.zeros(msi_shape) 

    if max_pixels > -1 and max_pixels < len(msi.coordinates):
        selected_pixels = list(rd.choice(len(msi.coordinates),
                                    max_pixels,
                                    replace=False))
    else:
        selected_pixels = list(range(len(msi.coordinates)))

    workers = (delayed(_annotate_pixel)(pixel_id, msi.getspectrum(pixel_id),
                               reference_mass_range, reference_spectra,
                               MTD, MTD_th) for pixel_id in selected_pixels)
    
    parallel = Parallel(backend='loky',
                        n_jobs=max_processes,
                        return_as='generator_unordered')
    annotator = parallel(workers)
    
    for pixel_id, pr_array in annotator:
        proportions += pr_array
        number_of_nonzero_pixels += pr_array > 0
        xcoord, ycoord, zcoord = msi.coordinates[pixel_id]
        control_images[ycoord-1, xcoord-1, ...] = pr_array[control_ids]
    proportions /= len(selected_pixels)
    return(proportions, number_of_nonzero_pixels, control_images)
    
