"""
This file contains additional functions to simplify the notebooks.
The functions are used to visualize data.
"""

import numpy as np
from matplotlib import pyplot as plt

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


def profile_ion_image(image, mz):
    """
    Returns an array of signal intensities at a given mz point from an MS image
    with spectra in profile mode.
    Uses a linear interpolation.
    Optionally, the intensity can be extracted after applying a gaussian filter with a given sigma.  
    """
    max_coord = max(image.coordinates)[:2]
    min_coord = min(image.coordinates)[:2]
    image_shape = (max_coord[1] - min_coord[1] + 1, max_coord[0] - min_coord[0] + 1)
    #print(image_shape)
    ion_image = np.zeros(image_shape)
    for idx, (xcoord,ycoord,zcoord) in enumerate(image.coordinates):
        mz_array, intsy =  image.getspectrum(idx)
        ion_image[ycoord-1,xcoord-1] = np.interp(mz, mz_array, intsy)
    return ion_image


def centroided_ion_image(image, mz, delta=0.05):
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


if __name__ == '__main__':
    from pyimzml.ImzMLParser import ImzMLParser
    bladder_centroided_image = ImzMLParser('MSimages/bladder_centroided.imzML')
    centroided_ion_image(bladder_centroided_image, 798.541)


