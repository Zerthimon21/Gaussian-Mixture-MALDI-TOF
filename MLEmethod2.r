rm(list = ls())
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/mzxml" 
setwd(FilePath)
source("IsotpoicFitLoad11.R") 

#Generate a model with no baseline noise or Ar(1) error
PepMass=MassCalc('DRVYIHPF')
PepIsotope <- IsotopicCountFFT('DRVYIHPF')  ##this must match peaks in figure
length(PepIsotope)
PepIsotope =PepIsotope/sum(PepIsotope)
sum(PepIsotope) 
mzError = .09
Sigma = .1
Range = seq(PepMass-(2*1.008664916),PepMass+(9*1.008664916),.02)
#This acts as if already standardized
PD <- MultiNormCalc(Range,'DRVYIHPF',PepMass,Sigma,mzError)
SpectraAverage <- PD[PD[,1]>=(PepMass-1.008664916) & PD[,1]<=(PepMass+(length(PepIsotope)*1.008664916)),]
plot(PD[,1],PD[,2],type='l')
plot(SpectraAverage[,1],SpectraAverage[,2],type='l') ##No normalization
range(SpectraAverage[,1])
##SpectraAverage[,2]=SpectraAverage[,2]/sum(SpectraAverage[,2])

#DeltaHat/mzError Estimator
#calculate expected value of kappa sum(point distrobution * isotope expected mass) This is for (K-1)N
ExK = 0
for(i in 1:length(PepIsotope)){
	ExK <- ExK + PepIsotope[i]*i 
}
#Sum of h,labda,x, including sum of hi
DoubleSum=0
for (i in 1:length(SpectraAverage[,2])){
	for(k in 1:length(PepIsotope)){
		DoubleSum <- DoubleSum+as.numeric((SpectraAverage[i,2]/sum(SpectraAverage[,2]))*PepIsotope[k]*SpectraAverage[i,1])
}
}
#mzError Estimator, delta hat
mzErrorHat = DoubleSum - PepMass - ExK*1.008664916 + 1.008664916 ##
mzErrorHat
(mzErrorHat-mzError)/mzError

########################################################################################################
##One peak of simulation
SpectraAverage <- PD[PD[,1]>=(PepMass+.5*1.008664916) & PD[,1]<=(PepMass+1.5*1.008664916),]
SSigma = 0
for (i in 1:length(SpectraAverage[,2])){
	bump <- SpectraAverage[i,2]*((SpectraAverage[i,1]-(PepMass+1.008664916+mzErrorHat))^2)
	bump <- bump/sum(SpectraAverage[,2])
	SSigma <- SSigma+bump
}
SSigma
Sigma=sqrt(SSigma)
Sigma

#Series of 7 peaks of simulation, note both of these miss *PepIsotope[k] in bump line
for(k in 1:length(PepIsotope)){
SpectraAverage <- PD[PD[,1]>=(PepMass+(k-1.5)*1.008664916) & PD[,1]<=(PepMass+(k-.5)*1.008664916),]
SSigma = 0
StdSum=sum(SpectraAverage[,2])
for (i in 1:length(SpectraAverage[,2])){
		bump <- SpectraAverage[i,2]*((SpectraAverage[i,1]-(PepMass+((k-1)*1.008664916)+mzErrorHat))^2)
		bump <- bump/StdSum
		SSigma <- SSigma+bump
}
##print(sum(SpectraAverage[,2]))
SSigma
Sigma=sqrt(SSigma)
print(Sigma)
}
##########################################################################################################

#SigmaSquared/variance,standard deviation/resolution/SigmaHat Estimator
SSigma = 0
SpectraAverage <- PD[PD[,1]>=(PepMass-(.5*1.008664916)) & PD[,1]<=(PepMass+((length(PepIsotope)-.5)*1.008664916)),]

for (i in 1:length(SpectraAverage[,2])){
	for(k in 1:length(PepIsotope)){
bump <- SpectraAverage[i,2]*PepIsotope[k]*((SpectraAverage[i,1]-(PepMass+((k-1)*1.008664916)+mzErrorHat))^2)
bump <- bump/sum(SpectraAverage[,2])
SSigma <- SSigma+bump
}}
SSigma/sum(SpectraAverage[,2])
sqrt(SSigma/sum(SpectraAverage[,2]))