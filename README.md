# spatialstein 

This repository contains an implementation of *spatialstein*, a workflow for segmentation of mass spectrometry imaging data. Along the workflow, we provide a script for simulating MSI data sets for testing purposes.

# The contents of this repository

The directory `Segmentation workflow` contains a set of scripts for segmentation of MSI data. Most scripts are in the form of Jupyter notebooks in Python, which allows the users a high degree of control over the intermediate results of the workflow to avoid errors propagating through the analysis. The final script is an R script. 

The directory `Image simulation` contains computer codes for simulating MSI data sets. 

# Reproducing the results from the article

### Results on simulated data 

In order to reproduce the results on the simulated data, you can either create a new data set using codes in `Image stimulation`, or use the data set available in the `Data` directory of `Segmentation workflow`. 

### Results on mouse bladder and cerebellum data

In order to reproduce these results, download the mouse bladder and cerebellum data sets.
Next, extract the archives if necessary and place the imzML images in the Data directory of `Segmentation workflow`. 
The scripts in `Segmentation workflow` are prepared so that they can be used directly on these data sets. 

# Using the workflow on other data sets

In order to use the workflow on your data in imzML format, place it in the Data directory (both the .imzML and the .ibd files). Next, modify the necessary variables in the scripts: add the path to your data as well as the name of your data set in appropriate cells of the Jupyter notebooks and the final R script. 

