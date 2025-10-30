import numpy as np
from numpy import random as rd
from matplotlib import pyplot as plt
from collections import Counter
from masserstein import Spectrum, estimate_proportions
from .data_preprocessing import get_msi_shape

def deconvolve_msi(msi,
                   reference_spectra,
                   MTD, MTD_th,
                   mask=None,
                   max_processes=1,
                   verbose=False, 
                   masserstein_kwargs = {}):
    """
    Performs a Wasserstein regression of the reference spectra against all pixel spectra of the image,
    ignoring coordinates for which mask = 0 (if supplied). 
    Returns a 3D array of estimated proportions. 
    The image is assumed to be in centroided mode. 
    All spectra are TIC-normalized prior to deconvolution.  
    """
    analyzed_mass_range = [
        min(s.confs[0][0] for s in reference_spectra) - 1,
        max(s.confs[-1][0] for s in reference_spectra) + 1
        ]
    if max_processes == 1:
        return deconvolve_msi_singlecore(
            msi = msi,
            reference_spectra = reference_spectra,
            MTD = MTD,
            MTD_th = MTD_th,
            mask = mask,
            analyzed_mass_range = analyzed_mass_range,
            verbose = verbose,
            masserstein_kwargs = masserstein_kwargs
            )
    else:
        from joblib import cpu_count
        if max_processes == -1:
            max_processes = cpu_count()
        return deconvolve_msi_parallel(
            msi = msi,
            reference_spectra = reference_spectra,
            MTD = MTD,
            MTD_th = MTD_th,
            mask = mask,
            analyzed_mass_range = analyzed_mass_range,
            max_processes = max_processes,
            verbose = verbose,
            masserstein_kwargs = masserstein_kwargs
            )
        

def deconvolve_msi_singlecore(msi, reference_spectra,
                              MTD, MTD_th, mask=None,
                              analyzed_mass_range = None,
                              verbose=False, 
                              masserstein_kwargs = {}):
    """
    Performs a Wasserstein regression of the reference spectra against all pixel spectra of the image,
    ignoring coordinates for which mask = 0 (if supplied). 
    Returns a 3D array of estimated proportions. 
    The image is assumed to be in centroided mode. 
    All spectra are TIC-normalized prior to deconvolution.  
    """
    # Get the image shape:
    i_coord_max = max(ycoord for xcoord,ycoord,zcoord in msi.coordinates) 
    i_coord_min = min(ycoord for xcoord,ycoord,zcoord in msi.coordinates) 
    i_range = i_coord_max - i_coord_min + 1
    j_coord_max = max(xcoord for xcoord,ycoord,zcoord in msi.coordinates) 
    j_coord_min = min(xcoord for xcoord,ycoord,zcoord in msi.coordinates) 
    j_range = j_coord_max - j_coord_min + 1
    
    # Get the mass range of the reference spectra 
    min_mz = min(s.confs[0][0] for s in reference_spectra) - 1
    max_mz = max(s.confs[-1][0] for s in reference_spectra) + 2
    
    # Deconvolve
    proportion_images = np.zeros((i_range, j_range, len(reference_spectra)))
    for idx, (xcoord,ycoord,zcoord) in enumerate(msi.coordinates):
        if verbose and not idx % 1000:
            print('Processing pixel number', idx)
        if mask is not None and mask[ycoord-1, xcoord-1] == 0: 
            continue
        # Get the pixel spectrum data 
        mz, intsy = msi.getspectrum(idx)
        
        # Truncate the pixel spectrum to the range of the reference spectra
        if analyzed_mass_range is not None:
            full_tic = sum(intsy)
            selected_range = (min_mz <= mz)*(mz <= max_mz)
            mz = mz[selected_range]
            intsy = intsy[selected_range]
            partial_tic = sum(intsy) 
            
        # Create a spectrum object, normalize for regression purposes
        pixel_spectrum = Spectrum(confs=list(zip(mz, intsy)))
        pixel_spectrum.normalize()
        
        # Regression
        regression = estimate_proportions(pixel_spectrum, reference_spectra, 
                                          MTD=MTD, MTD_th=MTD_th, 
                                          progress=False, **masserstein_kwargs)
        pr_array = np.array(regression['proportions'])
        
        if analyzed_mass_range is not None:
            # Rescale the proportions so that they reflect the proportion of the signal
            # in the full spectrum, not just the truncated part
            pr_array *= partial_tic/full_tic
        proportion_images[ycoord-1, xcoord-1, ...] = pr_array
    return proportion_images

def deconvolve_msi_parallel(msi, reference_spectra,
                            MTD, MTD_th, mask,
                            analyzed_mass_range,
                            max_processes,
                            verbose, 
                            masserstein_kwargs = {}):
    """
    Performs a Wasserstein regression of the reference spectra against all pixel spectra of the image,
    ignoring coordinates for which mask = 0 (if supplied). 
    Returns a 3D array of estimated proportions. 
    The image is assumed to be in centroided mode. 
    All spectra are TIC-normalized prior to deconvolution.  
    """
    from joblib import Parallel, delayed
    msi_shape = get_msi_shape(msi)
    proportion_images = np.zeros((msi_shape[0], msi_shape[1], len(reference_spectra)))

    # Create workers
    def worker_generator(analyzed_mass_range, reference_spectra, MTD, MTD_th):
        for idx, (xcoord,ycoord,zcoord) in enumerate(msi.coordinates):
            if mask is not None and mask[ycoord-1, xcoord-1] == 0: 
                continue
            mz, intsy = msi.getspectrum(idx)
            worker = delayed(_deconvolve_pixel)
            yield worker(idx,
                         mz,
                         intsy,
                         analyzed_mass_range,
                         reference_spectra,
                         MTD,
                         MTD_th)
            
    # Create the parallel object
    parallel = Parallel(backend='loky',
                        n_jobs=max_processes,
                        return_as='generator_unordered')
    # Run the workers
    workers = worker_generator(analyzed_mass_range, reference_spectra, MTD, MTD_th)
    annotator = parallel(workers)
    
    # Get the results
    for pixel_id, pr_array in annotator:
        xcoord, ycoord, zcoord = msi.coordinates[pixel_id]
        proportion_images[ycoord-1, xcoord-1, ...] = pr_array
    return proportion_images
        

def kappa_correlation_analysis(msi, mtd_values, mtd_th_values,
                               control_spectra, control_images,
                               number_of_pixels=1000,
                               max_processes=1):
    """
    Generate a heatmap of correlations between single-ion and deconvolved ion images
    for the control spectra and their single-ion images.
    Parameter combinations for which MTD_th < MTD are automatically ignored.  
    """
    assert len(control_spectra) == control_images.shape[2], 'Different numbers of control spectra and images'
    nspectra = len(control_spectra)
    msi_shape = get_msi_shape(msi)
    fit_mask = np.zeros(msi_shape, dtype='bool')
    sampled_pixel_ids = rd.choice(len(msi.coordinates), number_of_pixels, replace=False)
    for pxid in sampled_pixel_ids:
        rowcoord = pxid // msi_shape[1]
        colcoord = pxid % msi_shape[1]
        fit_mask[rowcoord,colcoord] = 1
        
    correlation_heatmaps = np.zeros((len(mtd_values), len(mtd_th_values), nspectra))
    for i, MTD in enumerate(mtd_values):
        for j, MTD_th in enumerate(mtd_th_values):
            if MTD_th < MTD: continue  
            deconvolved_images = deconvolve_msi(msi, control_spectra, mask=fit_mask,
                                               MTD=MTD, MTD_th = MTD_th,
                                               max_processes=max_processes,
                                               verbose=False,
                                               masserstein_kwargs = {'MDC':1e-12, 'MMD':.5,})

  
            # Below a manual calculation of the correlation that takes into account
            # the masking - masked pixels are ignored
            for lid in range(nspectra):
                masked_control = control_images[..., lid]*fit_mask
                deconv = deconvolved_images[..., lid]
                exy = np.sum(deconv*masked_control, axis=(0,1))/number_of_pixels
                ex = np.sum(deconv, axis=(0,1))/number_of_pixels
                ey = np.sum(masked_control, axis=(0,1))/number_of_pixels
                ex2 = np.sum(deconv**2, axis=(0,1))/number_of_pixels
                ey2 = np.sum(masked_control**2, axis=(0,1))/number_of_pixels
                covxy = exy - ex*ey
                sdx = np.sqrt(ex2 - ex**2)
                sdy = np.sqrt(ey2 - ey**2)
                correlation_heatmaps[i,j, lid] = covxy/(sdx*sdy)
    return correlation_heatmaps
            
def _deconvolve_pixel(pixel_id,
                      mz, intsy,
                      analyzed_mass_range, reference_spectra,
                      MTD, MTD_th,
                      masserstein_kwargs = {}):
    # Truncate the pixel spectrum to the range of the reference spectra
    if analyzed_mass_range is not None:
        min_mz, max_mz = analyzed_mass_range
        full_tic = sum(intsy)
        selected_range = (min_mz <= mz)*(mz <= max_mz)
        mz = mz[selected_range]
        intsy = intsy[selected_range]
        partial_tic = sum(intsy) 
            
    # Create a spectrum object, normalize for regression purposes
    pixel_spectrum = Spectrum(confs=list(zip(mz, intsy)))
    pixel_spectrum.normalize()
        
    # Regression
    regression = estimate_proportions(pixel_spectrum, reference_spectra, 
                                      MTD=MTD, MTD_th=MTD_th, 
                                      progress=False, **masserstein_kwargs)
    pr_array = np.array(regression['proportions'], dtype='float32')
        
    if analyzed_mass_range is not None:
        # Rescale the proportions so that they reflect the proportion of the signal
        # in the full spectrum, not just the truncated part
        pr_array *= partial_tic/full_tic
    return (pixel_id, pr_array)     

def segment_average_spectra(image, segment_mask, mass_axis):
    """
    Return a list of average spectra (vectors of intensities) corresponding to segments.
    segment_mask needs to be an array of the same shape as the image.  
    Values of this array denote segment IDs.
    image needs to be xy-indexed starting from 1 (i.e. (1, 1) is the bottom-left corner),
    while segment_mask needs to be array-indexed (i.e. segment_mask[0, 0] is the bottom-left corner)
    This is the default behavior if image is from ImzMLParser and segment_mask is plotted using plt.imshow.
    Image needs to be in profile mode. Spectra will be resampled in the points of the mass_axis.
    Spectra are normalized by TIC before summing within segments, and normalized by segment counts subsequently. 
    So basically everything is normalized correctly so that the output is bona fide average in segments.
    """
    segment_counter = Counter(segment_mask.flatten())
    segments = sorted(set(segment_counter))
    average_spectra = [np.zeros(len(mass_axis)) for _ in range(len(segments))]
    # generate a mapping between coordinates of image and segment_mask
    xy_coord_to_segment = {}
    for ar_i in range(segment_mask.shape[0]):
        for ar_j in range(segment_mask.shape[1]):
            xy_coord = ar_j + 1, ar_i + 1
            xy_coord_to_segment[xy_coord] = segment_mask[ar_i, ar_j]
    segment_name_to_id = {sg_n: sg_i for sg_i, sg_n in enumerate(segments)}
    for idx, (xcoord,ycoord,zcoord) in enumerate(image.coordinates):
        mz, intsy = image.getspectrum(idx)
        intsy = intsy / np.trapz(intsy, mz)
        intsy = np.interp(mass_axis, mz, intsy)
        sg_name = xy_coord_to_segment[(xcoord, ycoord)]
        average_spectra[segment_name_to_id[sg_name]] += intsy
    for sg in segments:
        average_spectra[segment_name_to_id[sg_name]] /= segment_counter[sg]
    return average_spectra


if __name__ == '__main__':
    from pyimzml.ImzMLParser import ImzMLParser
    from masserstein import Spectrum
    from masserstein import estimate_proportions
    cerebellum_centroided_image = ImzMLParser('MSimages/cerebellum_centroided.imzML')
    ctr_ion = centroided_ion_image(cerebellum_centroided_image, 798.541)
    cerebellum_profile_image = ImzMLParser('MSimages/test_POS.imzML')
    ion = profile_ion_image(cerebellum_profile_image, 798.541)
    segm = (ion > 0)*1 + (ion>3e06)*1
    mass_axis = np.arange(600, 800, step=0.001)
    av_sptra = segment_average_spectra(cerebellum_profile_image, segm, mass_axis)
    av_S = [Spectrum(confs=list(zip(mass_axis, intsy))) for intsy in av_sptra]
    for s in av_S: s.normalize()
    Spectrum.plot_all(av_S, profile=True)

##    TS1 = total_spectrum(cerebellum_profile_image, mass_axis, segm==1)
##    TS2 = total_spectrum(cerebellum_profile_image, mass_axis, segm==2)
##    av_sptra2 = [Spectrum(confs=list(zip(mass_axis, tsp))) for tsp in [TS1, TS2]]
##    for s in av_sptra2: s.normalize()
##    Spectrum.plot_all(av_sptra2, profile=True)
    
    test_lipid_names = ['PC(32:0)+K', 'PC(34:1)+K', 'PC(38:4)+K']
    test_lipid_formulas = ['C40H80NO8P', 'C42H82NO8P', 'C46H84NO8P']
    test_lipid_spectra = [Spectrum(formula=f, adduct='K', threshold=0.05) for f in test_lipid_formulas]
    test_lipid_masses = [s.confs[0][0] for s in test_lipid_spectra]
    for s in test_lipid_spectra: 
        s.normalize()
    test_images = analyze_image(cerebellum_centroided_image,
                                test_lipid_spectra,
                                MTD=0.1, MTD_th=0.15, verbose=True)
    plt.figure()
    plt.subplot(311)
    plt.imshow(test_images[...,0])
    plt.subplot(312)
    plt.imshow(test_images[...,1])
    plt.subplot(313)
    plt.imshow(test_images[...,2])
    plt.show()
