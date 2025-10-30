library(Cardinal)
set.seed(2025, kind="L'Ecuyer-CMRG")
setCardinalParallel(16)  # using 16 cores

### Input:
# Input the path to the data directory and the MSI filenames here:
data_dir <- './Deconvolved data/'
filenames <- c('lipid_MSI_deconvolved.imzML', 
               'mouse cerebellum deconvolved image.imzML', 
               'mouse bladder deconvolved image.imzML',
               'mouse brain deconvolved image.imzML')
result_dir <- './Segmentation maps/'
paths <- paste(data_dir, filenames, sep='')

### Loading the data:
msi_data_list <- lapply(paths, readMSIData)

# optional pre-processing; deconvolved images are already pre-processed
# normalized_data <- lapply(msi_data, normalize)
# normalized_data <- lapply(normalized_data, process)

# Visualizing to verify proper loading:
image(msi_data_list[[1]], i=1:3)
image(msi_data_list[[2]], mz=755.47, scale=T)
image(msi_data_list[[2]], i=1:12, scale=T)
image(msi_data_list[[3]], i=1:12, scale=T)
image(msi_data_list[[4]], i=60:66, scale=T)

### Checking the segmentation for different values of the r parameter
test_sdgmm <- spatialDGMM(msi_data_list[[1]], r=3, k=2, annealing=F, beta=6)
image(test_sdgmm, i=1:3)

image(msi_data_list[[2]], i=1:12, scale=T)
test_sdgmm <- spatialDGMM(msi_data_list[[2]], i=1:12, r=1, k=2, annealing=F, beta=4)
image(test_sdgmm, i=1:12)

image(msi_data_list[[3]], i=15:26, scale=T)
test_sdgmm <- spatialDGMM(msi_data_list[[3]], i=1:12, 
                          r=3, k=2, beta=6,
                          annealing=F, compress=F)
image(test_sdgmm, i=1:12)

image(msi_data_list[[4]], i=9:13, scale=T)
test_sdgmm <- spatialDGMM(msi_data_list[[4]], i=9:13, 
                          r=5, k=2, beta=8,
                          annealing=F, compress=F)
image(test_sdgmm, i=1:5)

# image(normalized_data[[3]], i=1:12, scale=T)
# test_sdgmm <- spatialDGMM(normalized_data[[3]], i=1:6, r=8, k=2, annealing=F, beta=20)

### Segmentation:
# spatialDGMM parameters for each image need to be specified manually here:
r_values <- c(3, 1, 3, 5)  
beta_values <- c(6, 4, 6, 8)
k_values <- c(2, 2, 2, 2)
sdgmm_list <- list()
processing_times <- list()
for(img_id in 1:length(filenames)){
  start.time <- Sys.time()
  sdgmm_list[[img_id]] <- spatialDGMM(msi_data_list[[img_id]], 
                                      r=r_values[img_id], 
                                      k=k_values[img_id],
                                      beta = beta_values[img_id],
                                      annealing = F,
                                      compress = F)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  processing_times[[img_id]] <- time.taken
}

# Visualizing to verify a proper segmentation:
image(sdgmm_list[[2]], i=1:12)
plot(sdgmm_list[[2]], i=3)
image(sdgmm_list[[3]], i=1:44)
image(sdgmm_list[[4]], i=6:14)

### Saving the results
# Saving the spatialdgmm objects for the segmentation of all MSI data sets:  
save(sdgmm_list, file=paste(result_dir, '/sdgmm_segmentation_results.RData', sep=''))

# Saving the spatialdgmm maps as arrays (optional, more portable between tools):
sdgmm_map_list <- list()
for(i in 1:length(filenames)){
  sdgmm_map <- sapply(sdgmm_list[[i]]$class, as.numeric)
  sdgmm_map <- data.frame(sdgmm_map)
  colnames(sdgmm_map) <- sdgmm_list[[i]]@featureData$mz
  coords <- data.frame('X' = sdgmm_list[[i]]@pixelData$x, 
                       'Y' = sdgmm_list[[i]]@pixelData$y)
  sdgmm_map <- cbind(coords, sdgmm_map)
  sdgmm_map_list[[i]] <- sdgmm_map
}

for(img_id in 1:length(filenames)){
  write.table(sdgmm_map_list[[img_id]], 
              file = paste(result_dir, sub('\\..*', '', filenames[img_id]), 
                           ' spatialdgmm map.csv', sep=''),
              row.names=F, col.names=F, sep='\t')
}
