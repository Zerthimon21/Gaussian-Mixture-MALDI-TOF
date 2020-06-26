rm(list = ls())
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/mzxml" 
setwd(FilePath)
source("IsotpoicFitLoad11.R") ##now includes library(xcms), AA

FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/SpectraSim" ##create folder 
setwd(FilePath)

SpectraSimulationParameters=as.data.frame(matrix(0, ncol = 8, nrow = 1))
colnames(SpectraSimulationParameters)=c('Name','Peptide','PeptideMass','Clusters','mzError','Sigma','BaselineNoise','PeakAreas')

PeptideMin=4
PeptideMax=15
SimNumber=100
PeakArea=c(1500,1000,750,500,250,100)

Peptide=NA
PeptideMass=NA

for(z in 1:SimNumber){
x = ceiling(runif(1,PeptideMin,PeptideMax))
pep=rep(0,x)
for (i in 1:x){
	pep[i] = AA[ceiling(runif(1,0,20))]
}
peptype=runif(1,0,1)
if(peptype>=.9){Cluster = 4} 
else if(peptype>=.8){Cluster = 2}
else if(peptype>=.7){Cluster = 2.5} ##strange peak offset, mimic 2nd compund
else {Cluster = 1}
SpectraSimulationParameters[z,1] = paste("FakePeptide",z)
SpectraSimulationParameters[z,2] = paste(pep,collapse = "")
SpectraSimulationParameters[z,3] = MassCalc(SpectraSimulationParameters[z,2])
SpectraSimulationParameters[z,4] = Cluster
SpectraSimulationParameters[z,5] = round(runif(1,-.4,.4),digits = 2)
SpectraSimulationParameters[z,6] = round(runif(1,.01,.3),digits = 2)
SpectraSimulationParameters[z,7] = round(runif(1,9,55),digits = 2)
SpectraSimulationParameters[z,8] = PeakArea[ceiling(runif(1,min=0,max=length(PeakArea)))]
}
SpectraSimulationParameters
write.csv(SpectraSimulationParameters,"SpectraSimulationParameters.csv")

SpectraSimulations=as.data.frame(matrix(0, ncol = 1, nrow = 505))
MZDense=50
##Noise Addition
rho=.7
UnifRange=2
Base=10
##Slice of simulated spectra that covers noise
SDback=2.5
Peaks=3		##N+Peaks of spectra to be shown
SDforward=2.5
for(i in 1:length(SpectraSimulationParameters[,1])){

PepRange = seq(SpectraSimulationParameters[i,3]-3*1.008664916,SpectraSimulationParameters[i,3]+7*1.008664916,1/MZDense)##~50 data points per Neutron mass

Base=rep(Base,(length(PepRange)+200)))
for(l in 1:(length(PepRange)+200)){
Base[l+1]=rho*Base[l]+runif(1,0,UnifBase)
}
NoiseAddition=(Base[100:length(PepRange)]) - mean(Base[100:length(PepRange)])+SpectraSimulationParameters[i,7]

PeakSample=MultiNormCalc(PepRange,SpectraSimulationParameters[i,2],SpectraSimulationParameters[i,3],SpectraSimulationParameters[i,6],SpectraSimulationParameters[i,5])
PeakSample[,2]=(PeakSample[,2]*SpectraSimulationParameters[i,8])+SpectraSimulationParameters[i,7]  ##+runif(length(PepRange),1/SpectraSimulationParameters[i,8],1/SpectraSimulationParameters[i,8]) ##peak noise...?

SpectraPeptideRange <- PeakSample[PeakSample[,1]>=(SpectraSimulationParameters[i,3]+SpectraSimulationParameters[i,5]-SDback*SpectraSimulationParameters[i,6]) & PeakSample[,1]<=(SpectraSimulationParameters[i,3]+SpectraSimulationParameters[i,5]+SDforward*SpectraSimulationParameters[i,6]+Peaks*1.008664916),]

for(m in 1:(length(PepRange))){
if(PeakSample[m,1]<(SpectraSimulationParameters[i,3]+SpectraSimulationParameters[i,5]-SDback*SpectraSimulationParameters[i,6]))
{PeakSample[m,2]=NoiseAddition[m]}
if(PeakSample[m,1]>(SpectraSimulationParameters[i,3]+SpectraSimulationParameters[i,5]+SDforward*SpectraSimulationParameters[i,6]+Peaks*1.008664916)))
{PeakSample[m,2]=NoiseAddition[m]}
}

colnames(PeakSample)=c('MZ','Intensity')
FileTitle=paste(SpectraSimulationParameters[i,1],'Spectra.csv',collapse = "")
write.csv(PeakSample,FileTitle)
##print(length(PeakSample[,2]))
SpectraSimulations=cbind(SpectraSimulations,PeakSample)
}
SpectraSimulations=SpectraSimulations[,-1]
write.csv(SpectraSimulations,"SpectraSimulations.csv")

##FileList <- list.files(FilePath,pattern='.*csv')

#####################################################################################################################################################
par(mfrow = c(3,1))

SpectraSimulationParameters[1,6]=.1
SpectraSimulationParameters[1,6]=.2
SpectraSimulationParameters[1,6]=.3

SpectraSimulationParameters[1,]
PepRange = seq(SpectraSimulationParameters[1,3]-3*1.008664916,SpectraSimulationParameters[1,3]+7*1.008664916,.02)##~50 data points per Neutron mass
PeakSample=MultiNormCalc(PepRange,SpectraSimulationParameters[1,2],SpectraSimulationParameters[1,3],SpectraSimulationParameters[1,6],SpectraSimulationParameters[1,5])
PeakSample[,2]=(PeakSample[,2]*SpectraSimulationParameters[1,8])+SpectraSimulationParameters[1,7]
plot(PeakSample[,1],PeakSample[,2],type='l')

SpectraSimulationParameters[1,6]=.2
##PeakSample = SpectraSimulations[,1:2]
Peaks = PeakSample[PeakSample[,1]>=(SpectraSimulationParameters[1,3]-(SpectraSimulationParameters[1,6])^2*2) & PeakSample[,1]<=(SpectraSimulationParameters[1,3]+4*1.008664916+(SpectraSimulationParameters[1,6])^2*4),]
plot(Peaks[,1],Peaks[,2],type='l')

Peaks = PeakSample[PeakSample[,1]>=(SpectraSimulationParameters[1,3]-(SpectraSimulationParameters[1,6])^2*2),]
Peaks = PeakSample[PeakSample[,1]<=(SpectraSimulationParameters[1,3]+3*1.008664916+(SpectraSimulationParameters[1,6])^2*4),]
##()^2 vs 
ErrorRange1 <- PeakSample[PeakSample[,1]<=(SpectraSimulationParameters[1,3]-SpectraSimulationParameters[1,6]*2) & PeakSample[,1]<=(SpectraSimulationParameters[1,3]+4*1.008664916+SpectraSimulationParameters[1,6]*4),]
runif(ReplacedLength,-1,1)
runif(ReplacedLength,-2,2)
runif(ReplacedLength,-(SpectraSimulationParameters[1,7]/10),(SpectraSimulationParameters[1,7]/10))

ErrorRange2 <- PeakSample[PeakSample[,1]>=(SpectraSimulationParameters[1,3]-SpectraSimulationParameters[1,6]*3) & PeakSample[,1]<=(SpectraSimulationParameters[1,3]+4*1.008664916+SpectraSimulationParameters[1,6]*4),]

ErrorRange3 <- PeakSample[PeakSample[,1]>=(SpectraSimulationParameters[1,3]-SpectraSimulationParameters[1,6]*4) & PeakSample[,1]<=(SpectraSimulationParameters[1,3]+4*1.008664916+SpectraSimulationParameters[1,6]*4),]

ErrorRange4 <- PeakSample[PeakSample[,1]>=(SpectraSimulationParameters[1,3]-SpectraSimulationParameters[1,6]*2) & PeakSample[,1]<=(SpectraSimulationParameters[1,3]+4*1.008664916+SpectraSimulationParameters[1,6]*2),]

x=seq(0,2,.02)
A=2 		##Amplitude
W=2*pi 	##Wavelength
y=A*sin(x*W)
plot(x,y)

x=seq(0,2,.02)
A=2 		##Amplitude
W=1*pi 	##Wavelength
y=A*sin(x*W)
plot(x,y,type='l')

x=seq(0,25,.02)
z=y+runif(length(y),-1,1)
plot(x,z)
###########################################################################################################
rm(list = ls())
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/mzxml" 
setwd(FilePath)
source("IsotpoicFitLoad11.R") ##now includes library(xcms), AA

FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/HEKAT1TP0"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')

SpectraImput <- xcmsRaw(FileList[1])
SpectraAverage <- getSpec(SpectraImput,mzrange=c(1045,1046.4))
length(SpectraAverage[,2])
plot(SpectraAverage[,1],SpectraAverage[,2],type='l')

SpectraAverage <- getSpec(SpectraImput,mzrange=c(1045,1046.6))
length(SpectraAverage[,2])
plot(SpectraAverage[,1],SpectraAverage[,2])

##AR(1)##
BL=10
BaselineSample1=BL
BaselineSample2=BL
BaselineSample3=BL
BaselineSample4=BL
rho=c(.7,.5,.2,-.3,-.51)

for(i in 1:700){
BaselineSample1[i+1]=rho*BaselineSample1[i]+runif(1,0,1)
BaselineSample2[i+1]=rho*BaselineSample2[i]+runif(1,0,BL/10)
BaselineSample3[i+1]=rho*BaselineSample3[i]+runif(1,0,BL/100)
BaselineSample4[i+1]=rho*BaselineSample4[i]+runif(1,0,2)

}
plot(seq(0,700,1),BaselineSample1,type='l',col='red',ylim=c(0,16))
abline(h=13.5,col='black')
lines(seq(0,700,1),BaselineSample2,type='l',col='green')
lines(seq(0,700,1),BaselineSample3,type='l',col='blue')
lines(seq(0,700,1),BaselineSample4,type='l',col='yellow')
##################################################################
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/HEKAT1TP0"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')

SpectraImput <- xcmsRaw(FileList[1])
SpectraAverage <- getSpec(SpectraImput,mzrange=c(1284,1294))

SmallAnalysis1=data.frame()
xyz=arima(SpectraAverage[,2],order=c(1,0,0))
xxx=unname(coef(xyz)[1])
is(xxx)
SmallAnalysis1[1,1]=coef(xyz)[1]
##################################################################
##Analysis of fitting
rm(list = ls())
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/mzxml" 
setwd(FilePath)
source("IsotpoicFitLoad11.R") ##now includes library(xcms), AA

FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/HEKAT1TP0"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')

SmallAnalysis1=c(0,0,0,0,0,0,0,0,0,0)
SmallAnalysis=c(0,0,0,0,0,0,0,0,0,0)
LargeAnalysis1=rep(0,length(FileList))

for (i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	LargeSpectraAverage <- getSpec(SpectraImput,mzrange=c(1194,1294))
	xyz1=arima(LargeSpectraAverage[,2],order=c(1,0,0))
	LargeAnalysis1[i]=as.numeric(as.character(coef(xyz1)[1]))
	for(j in 1:10){
		SmallSpectraAverage <- getSpec(SpectraImput,mzrange=c(1184+10*j,1194+10*j))
		xyz2=arima(SmallSpectraAverage[,2],order=c(1,0,0))
		SmallAnalysis[j]=as.numeric(as.character(coef(xyz2)[1]))
	}
	SmallAnalysis1=rbind(SmallAnalysis1,SmallAnalysis)
}

FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/Blanks"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')

SmallAnalysis2=c(0,0,0,0,0,0,0,0,0,0)
SmallAnalysis=c(0,0,0,0,0,0,0,0,0,0)
LargeAnalysis2=rep(0,length(FileList))

for (i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	LargeSpectraAverage <- getSpec(SpectraImput,mzrange=c(1194,1294))
	xyz1=arima(LargeSpectraAverage[,2],order=c(1,0,0))
	LargeAnalysis2[i]=coef(xyz1)[1]
	for(j in 1:10){
		SmallSpectraAverage <- getSpec(SpectraImput,mzrange=c(1184+10*j,1194+10*j))
		xyz2=arima(SmallSpectraAverage[,2],order=c(1,0,0))
		SmallAnalysis[j]=as.numeric(as.character(coef(xyz2)[1]))
	}
	SmallAnalysis2=rbind(SmallAnalysis2,SmallAnalysis)
}

FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/HEKAT1TP0"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')

SmallAnalysis3=c(0,0,0,0,0,0,0,0,0,0)
SmallAnalysis=c(0,0,0,0,0,0,0,0,0,0)
LargeAnalysis3=rep(0,length(FileList))

for (i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	LargeSpectraAverage <- getSpec(SpectraImput,mzrange=c(940,1040))
	xyz1=arima(LargeSpectraAverage[,2],order=c(1,0,0))
	LargeAnalysis3[i]=as.numeric(as.character(coef(xyz1)[1]))
	for(j in 1:10){
		SmallSpectraAverage <- getSpec(SpectraImput,mzrange=c(930+10*j,940+10*j))
		xyz2=arima(SmallSpectraAverage[,2],order=c(1,0,0))
		SmallAnalysis[j]=as.numeric(as.character(coef(xyz2)[1]))
	}
	SmallAnalysis3=rbind(SmallAnalysis3,SmallAnalysis)
}

LargeAnalysis3=rep(0,length(FileList))
for (i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	LargeSpectraAverage <- getSpec(SpectraImput,mzrange=c(946,1046))
	xyz1=arima(LargeSpectraAverage[,2],order=c(1,0,0))
	LargeAnalysis3[i]=as.numeric(as.character(coef(xyz1)[1]))
}


##HEKAT1TP0-Blank
SmallAnalysis1
LargeAnalysis1
##Blanks
SmallAnalysis2
LargeAnalysis2
##HEKAT1TP0-Near ANGII
SmallAnalysis3
LargeAnalysis3

x=rbind(SmallAnalysis1,LargeAnalysis1,SmallAnalysis2,LargeAnalysis2,SmallAnalysis3,LargeAnalysis3)
for (i in 1:length(x[,1])){
y=mean(x[i,])
print(y)
}
##############################################################################
##############################################################################
##############################################################################
##Printing
rm(list = ls())
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/mzxml" 
setwd(FilePath)
source("IsotpoicFitLoad11.R") 
##########
pdf(file='GraphicalNoiseModeling1.pdf')
par(mfrow=c(4,1))
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/HEKAT1TP0"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')
SmallAnalysis1=as.data.frame(matrix(0, ncol = 10, nrow = length(FileList)))
SmallAnalysis=c(0,0,0,0,0,0,0,0,0,0)
LargeAnalysis1=rep(0,length(FileList))
for (i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	LargeSpectraAverage <- getSpec(SpectraImput,mzrange=c(1194,1294))
	plot(LargeSpectraAverage[,1],LargeSpectraAverage[,2],type='l',main=FileList[i])
	for(j in 1:10){
		SmallSpectraAverage <- getSpec(SpectraImput,mzrange=c(1184+10*j,1194+10*j))
		plot(SmallSpectraAverage[,1],SmallSpectraAverage[,2],type='l',main=FileList[i])
	}
}
dev.off()

##########
pdf(file='GraphicalNoiseModeling2.pdf')
par(mfrow=c(4,1))
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/Blanks"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')
for (i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	LargeSpectraAverage <- getSpec(SpectraImput,mzrange=c(1194,1294))
	plot(LargeSpectraAverage[,1],LargeSpectraAverage[,2],type='l',main=FileList[i])
	for(j in 1:10){
		SmallSpectraAverage <- getSpec(SpectraImput,mzrange=c(1184+10*j,1194+10*j))
		plot(SmallSpectraAverage[,1],SmallSpectraAverage[,2],type='l',main=FileList[i])
	}
}
dev.off()

##########
pdf(file='GraphicalNoiseModeling3.pdf')
par(mfrow=c(4,1))
FilePath <- "C:/Documents and Settings/Zerthimon21/Desktop/spectra/HEKAT1TP0"
setwd(FilePath)
FileList <- list.files(FilePath,pattern='.*mzXML')
for (i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	
	SpectraAverage <- getSpec(SpectraImput,mzrange=c(1045,1051))
	plot(SpectraAverage[,1],SpectraAverage[,2],type='l',main=FileList[i])

	LargeSpectraAverage <- getSpec(SpectraImput,mzrange=c(946,1046))
	plot(LargeSpectraAverage[,1],LargeSpectraAverage[,2],type='l',main=FileList[i])
	for(j in 1:10){
		SmallSpectraAverage <- getSpec(SpectraImput,mzrange=c(936+10*j,946+10*j))
		plot(SmallSpectraAverage[,1],SmallSpectraAverage[,2],type='l',main=FileList[i])
		}
}
dev.off()