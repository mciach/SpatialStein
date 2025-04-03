library(Cardinal)
library(here)
set.seed(2025, kind="L'Ecuyer-CMRG")
setCardinalParallel(12)

### Input:
# Input the path to the data directory and the MSI filenames here:
data_dir <- here() 
filenames <- c('lipid_MSI_deconvolved.imzML', 
               'cerebellum_deconvolved_image.imzML', 
               'bladder_deconvolved_image.imzML')
paths <- paste(data_dir, filenames, sep='/')

### Loading the data:
msi_images <- lapply(filenames, readMSIData)
# Visualizing to verify proper loading:
image(msi_images[[1]], i=1:3)
image(msi_images[[2]], mz=755.47, scale=T)
image(msi_images[[2]], i=15, scale=T)
image(msi_images[[3]], i=1:12, scale=T)

### Checking the segmentation for different values of the r parameter
test_sdgmm <- spatialDGMM(msi_images[[1]], r=5, k=2, annealing=T)
image(test_sdgmm, i=1:3)
test_sdgmm <- spatialDGMM(msi_images[[2]], i=10:15, r=4, k=4, annealing=F)
image(test_sdgmm, i=1:6)
test_sdgmm <- spatialDGMM(msi_images[[3]], i=1:12, r=6, k=2)
image(test_sdgmm, i=1:12)

### Segmentation:
# spatialDGMM radius parameter for each image needs to be specified manually:
r_values <- c(4, 6)  
k_values <- c(2, 2)
sdgmm_list <- list()
for(img_id in 1:length(filenames)){
  sdgmm_list[[img_id]] <- spatialDGMM(images[[img_id]], 
                                      r=r_values[img_id], 
                                      k=k_values[img_id],
                                      annealing = F)
}
# Visualizing to verify a proper segmentation:
image(sdgmm_list[[2]], i=1:12)
plot(sdgmm_list[[2]], i=3)

### Saving the results
# Saving the spatialdgmm objects:  
save(sdgmm_list, file='sdgmm_segmentation_results.RData')

# Saving the spatialdgmm maps as arrays (optional, more portable between tools):
sdgmm_array_list <- lapply(sdgmm_list, function(x) sapply(x$class, as.numeric))
for(img_id in 1:length(filenames)){
  write.table(sdgmm_array_list[[img_id]], 
              file = paste(sub('\\..*', '', filenames[img_id]), 
                           '_spatialdgmm_map.csv', sep=''),
              row.names=F, col.names=F, sep='\t')
}
  