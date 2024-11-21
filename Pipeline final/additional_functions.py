"""
This file contains additional functions to simplify the notebooks.
The functions are used to visualize data.
"""

import numpy as np
from matplotlib import pyplot as plt
from collections import Counter
from masserstein import Spectrum, estimate_proportions

def plot_summary_spectrum(mass_axis, summary_spectrum_array, figsize=(8,4),
                          plot_title='Summary spectrum'):
    """
    Plots a spectrum showing the average intensity across the image
    along with the intensity's standard deviation, minimum, and maximum value.
    summary_spectrum_array is an array with four columns corresponding
    to the average, sd, min and max value.
    The number of rows in summary_spectrum_array needs to be equal to the length of mass_axis. 
    """
    assert np.all(summary_spectrum_array[:,0] >= 0), 'Intensity needs to be non-negative!'
    assert np.all(summary_spectrum_array[:,1] >= 0), 'Standard deviation needs to be non-negative!'
    assert np.all(summary_spectrum_array[:,2] <= summary_spectrum_array[:,3]), 'The minimum signal cannot be greater than the maximum!'
    av = summary_spectrum_array[:,0]
    sd = summary_spectrum_array[:,1]
    mn = summary_spectrum_array[:,2]
    mx = summary_spectrum_array[:,3]
    plt.figure(figsize=figsize)
    plt.title(plot_title)
    p1=plt.fill_between(mass_axis, av-sd, av+sd, alpha=1, color='orange')
    p2=plt.fill_between(mass_axis, av-3*sd, av+3*sd, alpha=0.2, color='orange')
    p3=plt.plot(mass_axis, av)
    p4=plt.plot(mass_axis, mx, 'k', alpha=0.4)
    p5=plt.plot(mass_axis, mn, 'k', alpha=0.4)
    plt.legend(['$\mu$', 'min', 'max', '$\mu \pm \sigma$', '$\mu \pm 3\sigma$'])
    plt.show()


def profile_ion_image(image, mz, normalize=True):
    """
    Returns an array of signal intensities at a given mz point from an MS image
    with spectra in profile mode.
    Uses a linear interpolation.
    Optionally, each pixel spectrum can be normalized to unit area under curve. 
    """
    max_coord = max(image.coordinates)[:2]
    min_coord = min(image.coordinates)[:2]
    image_shape = (max_coord[1] - min_coord[1] + 1, max_coord[0] - min_coord[0] + 1)
    #print(image_shape)
    ion_image = np.zeros(image_shape)
    for idx, (xcoord,ycoord,zcoord) in enumerate(image.coordinates):
        mz_array, intsy =  image.getspectrum(idx)
        if normalize:
            intsy = intsy / np.trapz(intsy, mz_array)
        ion_image[ycoord-1,xcoord-1] = np.interp(mz, mz_array, intsy)
    return ion_image


def centroided_ion_image(image, mz, delta=0.01, normalize=True):
    """
    Returns an array of signal intensities at a given mz point from an MS image
    with spectra in centroided mode.
    Uses a binary search to find the peak closest to mz, and then checks if it's within the mass
    accuracy given by the delta parameter (in Daltons).
    """
    max_coord = max(image.coordinates)[:2]
    min_coord = min(image.coordinates)[:2]
    image_shape = (max_coord[1] - min_coord[1] + 1, max_coord[0] - min_coord[0] + 1)
    #print(image_shape)
    ion_image = np.zeros(image_shape)
    for idx, (xcoord,ycoord,zcoord) in enumerate(image.coordinates):
        mz_array, intsy =  image.getspectrum(idx)
        if normalize:
            intsy = intsy / np.sum(intsy)
        assert all(nx>pr for nx, pr in zip(mz_array[1:], mz_array)), 'MZ array is not sorted'
        nbh = np.searchsorted(mz_array, mz)
        if nbh==0:
            if abs(mz - mz_array[0]) <= delta:
                ion_image[ycoord-1, xcoord-1] = intsy[0]
        elif nbh==len(mz_array):
            if abs(mz - mz_array[-1]) <= delta:
                ion_image[ycoord-1, xcoord-1] = intsy[-1]
        else:
            l, r = mz_array[nbh-1], mz_array[nbh]
            if mz - l < r - mz:
                # left neighbour is the closest one; update the closest neighbour index
                nbh -= 1
            if abs(mz - mz_array[nbh]) <= delta:
                ion_image[ycoord-1, xcoord-1] = intsy[nbh]
    return ion_image


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

def total_spectrum(image, mass_axis, mask=None, normalize_spectra=True):
    """
    Return a spectrum generated by summing all pixels in the image.
    The image needs to be in profile mode! Each pixel spectrum gets
    interpolated in the points of mass_axis.
    Spectra are normalized before interpolation. 
    If mask is supplied, it needs to be a binary or boolean array
    with the same shape as the image.
    Pixels corresponding to value zero or False are ignored.
    If mask is used, then the image needs to be xy-indexed starting
    from 1 (i.e. (1, 1) is the bottom-left corner; default e.g. for ImzMLParser),
    while segment_mask needs to be array-indexed (i.e. segment_mask[0, 0]
    is the bottom-left corner; default e.g. for plt.imshow).\
    
    """
    total_intensity = np.zeros(mass_axis.shape)
    skipped = 0
    for idx, (xcoord,ycoord,zcoord) in enumerate(image.coordinates):
        if mask is not None and not mask[ycoord-1, xcoord-1]:
            skipped += 1
            continue
        mz, intsy = image.getspectrum(idx)
        if normalize_spectra:
            intsy = intsy / np.trapz(intsy, mz)
        intsy = np.interp(mass_axis, mz, intsy)
        total_intensity += intsy
    #print('skipped', skipped)
    return total_intensity

def analyze_image(image, reference_spectra, MTD, MTD_th, 
                  analyzed_mass_range=None, mask=None, verbose=False, 
                  **maserstein_kwargs):
    """
    Performs a Wasserstein regression of the reference spectra against all pixel spectra of the image,
    ignoring coordinates for which mask = 0. 
    Truncates the spectra to the analyzed mass range to speed up computations. 
    Returns a 3D array of estimated proportions. 
    """
    # Get the image shape:
    i_coord_max = max(ycoord for xcoord,ycoord,zcoord in image.coordinates) 
    i_coord_min = min(ycoord for xcoord,ycoord,zcoord in image.coordinates) 
    i_range = i_coord_max - i_coord_min + 1
    j_coord_max = max(xcoord for xcoord,ycoord,zcoord in image.coordinates) 
    j_coord_min = min(xcoord for xcoord,ycoord,zcoord in image.coordinates) 
    j_range = j_coord_max - j_coord_min + 1
    proportion_images = np.zeros((i_range, j_range, len(reference_spectra)))
    for idx, (xcoord,ycoord,zcoord) in enumerate(image.coordinates):
        if verbose and not idx % 1000:
            print('Processing pixel number', idx)
        if mask is not None and mask[ycoord-1, xcoord-1] == 0: 
            continue
        # Get the pixel spectrum data and truncate
        mz, intsy = image.getspectrum(idx)
        if analyzed_mass_range is not None:
            selected_range = (analyzed_mass_range[0] <= mz)*(mz <= analyzed_mass_range[1])
            mz = mz[selected_range]
            intsy = intsy[selected_range]
        # Create a spectrum object, normalize for regression purposes
        pixel_spectrum = Spectrum(confs=list(zip(mz, intsy)))
        pixel_spectrum.normalize()
        # Regress
        regression = estimate_proportions(pixel_spectrum, reference_spectra, 
                                          MTD=MTD, MTD_th=MTD_th, 
                                          progress=False, **maserstein_kwargs)
        pr_array = np.array(regression['proportions'])
        proportion_images[ycoord-1, xcoord-1, ...] = pr_array
    return proportion_images

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
