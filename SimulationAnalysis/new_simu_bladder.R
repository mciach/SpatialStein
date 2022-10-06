library('Cardinal')
mono = readMSIData('/Users/danguo/Documents/phd/SpatialStein/NewSimulations/lipid_MSI_centroid_mode.imzML')
spec = read.table('/Users/danguo/Documents/phd/SpatialStein/NewSimulations/lipid_regression_results.tsv', sep = '\t', header = T, row.names = NULL)

deconv = mono[1:3,]
spectra(deconv) = t(spec[,3:5])

for (i in 1:3)
{
  int = matrix(as.vector(t(spec[,(i+2)])),nrow=40, byrow = TRUE)
  int<-c(int)
  spectra(deconv)[i,] = int
}

mass_list = read.table('/Users/danguo/Documents/phd/SpatialStein/NewSimulations/lipid_table.tsv', sep = '\t', header = T, row.names = NULL)
mass_list = mass_list$Monoisotopic.mass

f_list =  c()
for (i in 1:3)
{
  f_list = c(f_list, which(abs(mz(mono)-mass_list[i])==min(abs(mz(mono)-mass_list[i]))))
}

par(mfrow = c(2,2))
dgmm <- spatialDGMM(deconv, r=10, k=2, method="adaptive", init = 'gmm')
for (i in 1:3)
{
  print(image(dgmm, model=list(feature=i), values = 'class', layout=NULL, main=paste0('deconv ')))
}


source('~/Documents/phd/SpatialStein/DSDGMM/model.R')
sp_ratio=1
radius=2
w<-W_matrix(mono,sp_ratio=sp_ratio,radius=radius)
w2<-w[[2]]
w<-w[[1]]

sdgmm<-DGMM(msset=deconv, w=w, w2=w2, k=2, f=1, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")

xx<-apply(sdgmm[[10]],1, function (x) which(x==max(x)))
image(deconv, sdgmm[[11]]~x*y)
skm@resultData@listData[[1]]$cluster<-as.factor(xx)

image(skm)

load('data/Andrew_bladder.rdata')
bladder = mse_peaks
bladder <- as(bladder, "MSProcessedImagingExperiment")
spec = read.table('/Users/danguo/Documents/phd/SpatialStein/MouseBladderAnalysis/cluster3_regression_results.tsv', sep = '\t', header = T, row.names = NULL)
deconv_bladder = bladder[1:5,]

spectra(deconv_bladder) = t(spec[,3:7])
for (i in 1:5)
{
  int = matrix(as.vector(t(spec[,(i+2)])),nrow=260, byrow = TRUE)
  int<-c(int)
  spectra(deconv_bladder)[i,] = int
}


sp_ratio=1
radius=2
w<-W_matrix(deconv_bladder,sp_ratio=sp_ratio,radius=radius)
w2<-w[[2]]
w<-w[[1]]

sdgmm<-DGMM(msset=deconv_bladder, w=w, w2=w2, k=4, f=1, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")

xx<-as.factor(apply(sdgmm[[10]],1, function (x) which(x==max(x))))
image(deconv_bladder, xx~x*y)



coords <- coord(deconv_bladder)
coords <- as.data.frame(coords)

sdgmm<-DGMM(msset=deconv_bladder, w=w, w2=w2, k=3, f=1, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")
xx<-as.factor(apply(sdgmm[[10]],1, function (x) which(x==max(x))))
image(deconv_bladder, xx~x*y)

coords[,3] = xx

sdgmm<-DGMM(msset=deconv_bladder, w=w, w2=w2, k=2, f=2, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")
xx<-as.factor(apply(sdgmm[[10]],1, function (x) which(x==max(x))))
image(deconv_bladder, xx~x*y)

coords[,4] = xx


sdgmm<-DGMM(msset=deconv_bladder, w=w, w2=w2, k=4, f=3, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")
xx<-as.factor(apply(sdgmm[[10]],1, function (x) which(x==max(x))))
image(deconv_bladder, xx~x*y)
coords[,5] = xx

sdgmm<-DGMM(msset=deconv_bladder, w=w, w2=w2, k=2, f=4, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")
xx<-as.factor(apply(sdgmm[[10]],1, function (x) which(x==max(x))))
image(deconv_bladder, xx~x*y)
coords[,6] = xx



sdgmm<-DGMM(msset=deconv_bladder, w=w, w2=w2, k=3, f=5, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")
xx<-as.factor(apply(sdgmm[[10]],1, function (x) which(x==max(x))))
image(deconv_bladder, xx~x*y)
coords[,7] = xx

write.csv(coords, file = 'sdgmm_bladder_lowr.csv')

pdf('bladder_lowr_sdgmm.pdf')
par(mfrow=c(2,3))
for (i in 1:5)
{
  print(image(deconv_bladder,coords[,i+2]~x*y, layout=NULL))
}
dev.off()

spectra(deconv_bladder) <- spectra(deconv_bladder)*100000
out <- 10000
f_list = 1:5
k_list = c(3,2,4,2,3)
segments_2r <- matrix(0,nrow = ncol(deconv_bladder), ncol = length(f_list))
for (i in 1:5)
{
  sdgmm<-DGMM(msset=deconv_bladder, w=w, w2=w, k=k_list[i], f=f_list[i], sp_ratio=sp_ratio,step=1e8,iteration=200,Annl=0, sig=0.3,initialization="km")
  segments_2r[,i] <- sdgmm[[11]]
  #print(image(msset, sdgmm[[11]]~x*y,layout=NULL))
}

pdf('seg_2r_bladder.pdf')
par(mfrow=c(2,2))

for (i in 1:length(f_list))
{
  print(image(deconv_bladder, as.factor(segments_2r[,i])~x*y,layout=NULL, main = paste0('mz = ', f_list[i])))
}
dev.off()

write.csv(segments_2r, file = 'segment_2r_bladder.csv')

################k-means
k_list <- c(3,2,4,2,3)

for (i in 1:5)
{
  km = kmeans(spectra(deconv_bladder)[i,], centers = k_list[i])
  coords[,i+2] = as.factor(km$cluster)
}

write.csv(coords, file = 'km_bladder_lowr.csv')


pdf('bladder_lowr_km.pdf')
par(mfrow=c(2,3))
for (i in 1:5)
{
  print(image(deconv_bladder,coords[,i+2]~x*y, layout=NULL))
}
dev.off()

#################


pdf('bladder_image.pdf')
par(mfrow=c(2,3))
for (i in 1:5)
{
  print(image(deconv_bladder,feature = i, layout=NULL))
}
dev.off()

i=1
for (f in f_list)
{
  km <- kmeans(spectra(mono)[f,], centers = 2)
  coords[,i+2] <- ifelse(km$cluster==1,2,1)
  i = i+1
}

colnames(coords) <- c('x', 'y', 'PC(38:1)', 'PA(44:0)', 'PC(38:0)')

image(mono,coords$`PC(38:0)`~x*y)

write.csv(coords, file = 'km_simulation.csv')


coords <- coord(mono)
coords <- as.data.frame(coords)

i=1
for (f in f_list)
{
  sdgmm<-DGMM(msset=mono, w=w, w2=w, k=2, f=f, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")
  
  coords[,i+2] <- sdgmm[[11]]
  i = i+1
}


colnames(coords) <- c('x', 'y', 'PC(38:1)', 'PA(44:0)', 'PC(38:0)')
write.csv(coords, file = 'sdgmm_simulation_5r.csv')



coords <- coord(deconv)
coords <- as.data.frame(coords)

i=1
for (f in 1:3)
{
  sdgmm<-DGMM(msset=deconv, w=w, w2=w, k=2, f=f, sp_ratio=sp_ratio,step=1e10,iteration=1000,Annl=0, sig=0.3,initialization="km")
  
  coords[,i+2] <- sdgmm[[11]]
  i = i+1
}


colnames(coords) <- c('x', 'y', 'PC(38:1)', 'PA(44:0)', 'PC(38:0)')
write.csv(coords, file = 'sdgmm_simulation_deconv_5r.csv')


