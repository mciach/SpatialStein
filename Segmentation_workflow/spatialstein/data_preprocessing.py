def shift_coordinates(img):
    x_crd, y_crd, z_crd = zip(*img.coordinates)
    if len(set(z_crd)) > 1:
        raise ValueError('Multiple Z coordinates detected.\nSpatialStein currently works only with 2D images.')
    min_x, min_y = min(x_crd), min(y_crd)
    for i, crd in enumerate(img.coordinates):
        crd = (crd[0]-min_x, crd[1]-min_y, 1) 
        img.coordinates[i] = crd

def get_msi_shape(msi):
    x_crd, y_crd, *_ = zip(*msi.coordinates)
    return(max(y_crd)-min(y_crd)+1, max(x_crd)-min(x_crd)+1)

def centroid_and_save_msi_dataset(msi, output_path, peak_height_fraction, max_width):
    """
    Goes through each pixel, centroids it, and saves the result
    in an .imzML file in output_path.
    Memory-efficient implementation - doesn't load the whole data set.
    """
    from masserstein import Spectrum
    from pyimzml.ImzMLWriter import ImzMLWriter

    with ImzMLWriter(output_path) as writer:
        for idx, (xcoord,ycoord,zcoord) in enumerate(msi.coordinates):
            if not idx % 10000:
                print('Processing pixel number', idx)
            mz, intsy = msi.getspectrum(idx)
            S = Spectrum(confs=list(zip(mz, intsy)))
            peaks, _ = S.centroid(peak_height_fraction=peak_height_fraction, 
                                  max_width=max_width)
            mzs = [p[0] for p in peaks]
            intsys = [p[1] for p in peaks]
            writer.addSpectrum(mzs, intsys, (xcoord, ycoord))
