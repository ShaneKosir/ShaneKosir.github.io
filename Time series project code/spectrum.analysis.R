library(stringr)
library(mice)
library(quantmod)

# Import Dataset ----------------------------------------------------------
## Import Time Series Dataset that you want to develop a power spectrum for.

# Select Column Number for Power Spectrum ----------------------------------------------------------
data<-data_preprocessed

# Scale column for power spectrum
x<-na.omit(data$carbon)
x.sd <- sd(x)
maxValue <- max(x)
minValue <- min(x)
avg <- mean(x)
  
data.col.for.power.spectrum <- as.data.frame(scale(x,center=(avg),scale = (maxValue-minValue)))
data.col.for.power.spectrum <- data.matrix(data.col.for.power.spectrum)
summary(data.col.for.power.spectrum)
spx <- rep(NA,15000)
spy <- rep(NA,15000)

power.spectrum.results <- spectrum(x,span = 2)
power.spectrum.results.strength <- t(power.spectrum.results$spec)
del<- 1 # sampling interval, 1 year--> rad/year
  
spx <- power.spectrum.results$freq/del
spy <- 2*power.spectrum.results$spec
summary(spx)
summary(spy)
options(scipen=5) # turn scientific notation off
plot(spx,spy,xlab="Frequency, 1/year",ylab="Spectral Density",type="l",
     xlim=c(0,.5),ylim=c(0,2.5e17))

# Find Peaks ----
peak.indices<-findPeaks(spy)
peak.indices<-peak.indices-1
peak.frequencies<-spx[peak.indices]
peak.periods<-1/peak.frequencies
peak.spectral.densities<-spy[peak.indices]
length.peaks=length(peak.indices)

# Compile Power Spectrum ----
spectral.results<-matrix(0,ncol=6,nrow=length(spx))
spectral.results<-data.frame(spectral.results)
names(spectral.results)<-c('Frequency, 1/yr','Spectral Density','Peak Frequencies, 1/yr',
                           'Peak Periods, yr','Peak Spectral Densities','Peak Indices')
spectral.results[1]<-spx
spectral.results[2]<-spy
spectral.results[1:length.peaks,3]<-peak.frequencies
spectral.results[1:length.peaks,4]<-peak.periods
spectral.results[1:length.peaks,5]<-peak.spectral.densities
spectral.results[1:length.peaks,6]<-peak.indices

write.csv(spectral.results,"spectral.results.csv",row.names=FALSE)