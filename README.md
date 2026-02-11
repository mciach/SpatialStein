# spatialstein 

This repository contains an implementation of *spatialstein*, a workflow for segmentation of mass spectrometry imaging data. Along the workflow, we provide a script for simulating MSI data sets for testing purposes.

# The contents of this repository

The directory `Segmentation workflow` contains a set of scripts for segmentation of MSI data. Most scripts are in the form of Jupyter notebooks in Python, which allows the users a high degree of control over the intermediate results of the workflow to avoid errors propagating through the analysis. The final script is an R script. 

The directory `Image simulation` contains computer codes for simulating MSI data sets. 

# Reproducing the results from the article

### Results on simulated data 

In order to reproduce the results on the simulated data, you can either create a new data set using codes in `Image stimulation`, or use the data set available in the `Data` directory of `Segmentation workflow`. 

### Results on mouse bladder and cerebellum data

The annotation results, deconvolved images, and segmentation maps described in the Spatialstein article are available on Zenodo at:
```
https://zenodo.org/records/18613414
```
In order to reproduce these results, download the mouse bladder and cerebellum data sets:

     - Mouse bladder:   https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001283
     - Mouse cerebellum:   https://www.ebi.ac.uk/metabolights/editor/MTBLS487
     
Next, extract the archives if necessary and place the imzML images in the Data directory of `Segmentation workflow`. 
The scripts in `Segmentation workflow` are prepared so that they can be used directly on these data sets. 


# Using the workflow on other data sets

In order to use the workflow on your data in imzML format, place it in the Data directory (both the .imzML and the .ibd files). Next, modify the necessary variables in the scripts: add the path to your data as well as the name of your data set in appropriate cells of the Jupyter notebooks and the final R script. 

# How to cite
If you use `spatialstein`, please cite the following article:
```
Ciach, M. A., Guo, D., Bemis, K. A., Valkenborg, D., Vitek, O., & Gambin, A. (2026). spatialstein: An Open-Source Workflow for Annotation, Deconvolution, and Spatially Aware Segmentation of Mass Spectrometry Imaging Data. Analytical Chemistry.
```
The article is available under this link:    
https://pubs.acs.org/doi/10.1021/acs.analchem.5c04737
