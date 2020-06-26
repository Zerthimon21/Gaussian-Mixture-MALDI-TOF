##reworked Pepnfo for easier table reconstruction from v8, more user friendly 
##and allows for easier change of peptides outside angiotensin searches of mass spectra

library(xcms)
###########################################################
UnlabeledSequence = NA
UnlabeledPeptides = NA
LabeledSequence = NA
LabeledPeptides = NA

###########################################################
#Construct AAcomp globally used as a lookup tabel for amino acid compositions, avoid repeated reconstruction in loops
	AAComp=rbind(c(3,2,1,1,0,0),c(5,3,1,1,0,0),c(5,3,1,2,0,0),c(7,5,1,1,0,0),c(9,5,1,1,0,0),c(7,4,1,2,0,0),c(5,3,1,1,0,1),c(11,6,1,1,0,0),c(11,6,1,1,0,0),c(6,4,2,2,0,0),c(5,4,1,3,0,0),c(8,5,2,2,0,0),c(12,6,2,1,0,0),c(7,5,1,3,0,0),c(9,5,1,1,0,1),c(7,6,3,1,0,0),c(9,9,1,1,0,0),c(12,6,4,1,0,0),c(9,9,1,2,0,0),c(10,11,2,1,0,0))
	rownames(AAComp)=c('G','A','S','P','V','T','C','L','I','N','D','Q','K','E','M','H','F','R','Y','W')
	colnames(AAComp)=c('H','C','N','O','P','S')
	AA=c('G','A','S','P','V','T','C','L','I','N','D','Q','K','E','M','H','F','R','Y','W')
###########################################################
#Construct matrix of atom totals for a given peptide, includes a water molecule
AtomMatrix=function(pep, Water=TRUE){
	AtomMatrix = matrix(0,nchar(pep)+1,6)  ##+1 for extra water at ends of peptide/ single amino acid

	for (i in 1:nchar(pep)){
		AtomMatrix[i,]= AAComp[substring(pep,i,i),] ##Extract substrings in a character vector 
		}
		if(Water){
	AtomMatrix[nchar(pep)+1,]=c(2,0,0,1,0,0)  ##Add water to count
			}
	AtomMatrix
}
###########################################################
##Gives monotopic mass (isotope 0) of a given peptide summed as H,C,N,O,P,S
MassCalc <- function(pep, AAComp, proton=1, MH=TRUE){
	AtomMatrix <- AtomMatrix(pep)
 
	Mass <- sum(AtomMatrix[,1])*1.007825032+sum(AtomMatrix[,2])*12+sum(AtomMatrix[,3])*14.003074004+sum(AtomMatrix[,4])*15.994914619+sum(AtomMatrix[,5])*30.973761632+sum(AtomMatrix[,6])*31.972071002##CHECK CRC
	if(MH){
		Mass <- (Mass+proton*1.007276467)/proton##Note peaks become(MassNuetron)/#Protons apart
	}
	Mass
}
###########################################################
##Generates point distrobutins for a given combination of atoms (peptide)
IsotopicCountFFT = function (pep,cutoff=.9999){

	AtomMatrix=AtomMatrix(pep)

	Molecule=c(H=sum(AtomMatrix[,1]),C=sum(AtomMatrix[,2]),N=sum(AtomMatrix[,3]),O=sum(AtomMatrix[,4]),P=sum(AtomMatrix[,5]),S=sum(AtomMatrix[,6]))

	##check for all atoms being present using molecule to check for 0's
	#ScanMolecule = Molecule[Molecule!=0]

	##CRC values for isotopic point distrobutions for each element
	hd = c(0.99985, 0.00015) ##1,2  
	cd = c(0.98900, 0.01100) ##12,13
	nd = c(0.99634, 0.00366) ##14,15
	od = c(0.99760, 0.00039, 0.00201) ##16,17,18
	pd = c(1.00000) ##31
	sd = c(0.95020, 0.00750, 0.04210, 0.00000, 0.00020) ##32,33,34,35,36  CHANGETHIS!!!!

	## Counting each occurance of a given atom in the selected molecule
	hcount = Molecule['H'] 
	ccount = Molecule['C']
	ncount = Molecule['N']
	ocount = Molecule['O']
	pcount = Molecule['P']
	scount = Molecule['S']

	##Creating vector for Fast Fourier Transform based on number of atoms of a given type present
	##FFT work faster on vectors that have a length of 2^x, use 256 here since we look at low dalton (mass) proteins
	##look at 128 innstead of 256 for speed
 
	MoleculeFFT=1
	if (hcount>0){
		hd2 = append(hd,rep(0,256-length(hd)))
		hd2=fft(hd2)^hcount
		MoleculeFFT=MoleculeFFT*hd2
	}
	if (ccount>0){
		cd2 = append(cd,rep(0,256-length(cd)))
		cd2 = fft(cd2)^ccount
		MoleculeFFT=MoleculeFFT*cd2
	}
	if (ncount>0){
		nd2 = append(nd,rep(0,256-length(nd)))
		nd2 = fft(nd2)^ncount
		MoleculeFFT=MoleculeFFT*nd2
	}
	if (ocount>0){
		od2 = append(od,rep(0,256-length(od)))
		od2 = fft(od2)^ocount
		MoleculeFFT=MoleculeFFT*od2
	}
	if (pcount>0){
		pd2 = append(pd,rep(0,256-length(pd)))
		pd2 = fft(pd2)^pcount
		MoleculeFFT=MoleculeFFT*pd2
	}
	if (scount>0){
		sd2 = append(sd,rep(0,256-length(sd)))
		sd2 = fft(sd2)^scount
		MoleculeFFT=MoleculeFFT*sd2
	}

	distri = fft(MoleculeFFT,inverse=TRUE)/length(MoleculeFFT)

	##Remove the complex part of each entry in dristi
	distri = Re(distri)

	# Remove those peaks (higher masses) that contribute less than 1/100% of the total mass
	stop = min(which(cumsum(distri)>cutoff))
	distri[1:stop]

	###distri
	##This produces slight differences from old convolution method due rounding errors to transformation, power raise, inverse transformation
	##still produces - negnumbers use rnd if this becomes a problem
}
#############################################
#############################################
##Generate a Mixture of Normals model, Mean separated by one Neutron standars deviations 
##are equal this curve is used to model the peaks of a mass spectrum
MultiNormCalc=function(SpectraRange,peptide,mass,sigma,mzError=0){

	Mass=mass ##This calulated previously
	h=IsotopicCountFFT(peptide)
	mz=SpectraRange#x-axis, the mass charge
	peaks=rep(0,length(mz))
	for(i in 1:length(h)){
		peaks=peaks+dnorm(mz-Mass,((i-1)*1.008664916)+mzError,sigma)*h[i]
	}
	Predicted=matrix(0,length(mz),2)
	Predicted[,1]=mz
	Predicted[,2]=peaks
	Predicted
}
#############################################
PeptideInfo2=function(UnlabeledPeptides,UnlabeledSequence,LabeledPeptides,LabeledSequence,AQUAState="Both",Neutron=0,ClusterCut=5){
switch(AQUAState,
	Both = {##Makes list for both labeled and unlabeled peptides
	Peptides=cbind(UnlabeledSequence)
	rownames(Peptides)=UnlabeledPeptides
	HeavyPeptides=cbind(LabeledSequence)
	rownames(HeavyPeptides)=LabeledPeptides
	
PeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjustSlope=NA,NoiseAdjustIntercept=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

	for(i in 1:length(Peptides)){ 
	PeptideInformation[i,1] <- rownames(Peptides)[i]
	PeptideInformation[i,2] <- Peptides[i]
	PeptideInformation[i,3] <- MassCalc(Peptides[i]) 
	}
for(i in 1:length(HeavyPeptides)){ 	
	PeptideInformation[i+length(Peptides),1] <- paste('AQUA',rownames(HeavyPeptides)[i])
	PeptideInformation[i+length(Peptides),2] <- HeavyPeptides[i]
	PeptideInformation[i+length(Peptides),3] <- MassCalc(HeavyPeptides[i])+Neutron*1.008664916
}

##Order peptides based on mass
PeptideInformation <- PeptideInformation[order(-PeptideInformation$Mass),]

##Assign clusters to each group of peaks less than ClusterCut mz apart
counter=0
for(i in 1:length(PeptideInformation$Mass)){
	if(i<=length(PeptideInformation$Mass)-1){
		if ((PeptideInformation$Mass[i]-PeptideInformation$Mass[i+1])<ClusterCut){
			PeptideInformation$Cluster[i]=counter+1
			PeptideInformation$Cluster[i+1]=counter+1
		}
		else{
			PeptideInformation$Cluster[i]=counter+1
			PeptideInformation$Cluster[i+1]=counter+2
			counter=counter+1}
	}
	else
		if ((PeptideInformation$Mass[i-1]-PeptideInformation$Mass[i])<ClusterCut){
			PeptideInformation$Cluster[i]=counter+1
		}
		else{
			PeptideInformation$Cluster[i]=counter+1
		}
}
PeptideInformation		
		}, 
	AQUA = {##Makes list for just labeled peptides
	HeavyPeptides=cbind(LabeledSequence)
	rownames(HeavyPeptides)=LabeledPeptides
	
PeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjustSlope=NA,NoiseAdjustIntercept=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)


for(i in 1:length(HeavyPeptides)){
			PeptideInformation[i,1] <- paste('AQUA',rownames(HeavyPeptides)[i])
			PeptideInformation[i,2] <- HeavyPeptides[i]
			PeptideInformation[i,3] <- MassCalc(HeavyPeptides[i])+Neutron*1.008664916
	}


##Order peptides based on mass
PeptideInformation <- PeptideInformation[order(-PeptideInformation$Mass),]

##Assign clusters to each group of peaks less than ClusterCut mz apart
counter=0
if(length(PeptideInformation$Mass) < 2){
PeptideInformation$Cluster[1]=1
}
else{
for(i in 1:length(PeptideInformation$Mass)){
	if(i<=length(PeptideInformation$Mass)-1){
		if ((PeptideInformation$Mass[i]-PeptideInformation$Mass[i+1])<ClusterCut){
			PeptideInformation$Cluster[i]=counter+1
			PeptideInformation$Cluster[i+1]=counter+1
		}
		else{
			PeptideInformation$Cluster[i]=counter+1
			PeptideInformation$Cluster[i+1]=counter+2
			counter=counter+1}
	}

	else
		if ((PeptideInformation$Mass[i-1]-PeptideInformation$Mass[i])<ClusterCut){
			PeptideInformation$Cluster[i]=counter+1
		}
		else{
			PeptideInformation$Cluster[i]=counter+1
		}
}}
PeptideInformation
		}, 
	Normal = {##Makes list for just unlabeled peptides
	Peptides=cbind(UnlabeledSequence)
	rownames(Peptides)=UnlabeledPeptides
	
PeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjustSlope=NA,NoiseAdjustIntercept=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

for(i in 1:length(Peptides)){ 
	PeptideInformation[i,1] <- rownames(Peptides)[i]
	PeptideInformation[i,2] <- Peptides[i]
	PeptideInformation[i,3] <- MassCalc(Peptides[i]) 
}

##Order peptides based on mass
PeptideInformation <- PeptideInformation[order(-PeptideInformation$Mass),]

##Assign clusters to each group of peaks less than ClusterCut mz apart
counter=0
if(length(PeptideInformation$Mass) < 2){
PeptideInformation$Cluster[1]=1
}
else{
for(i in 1:length(PeptideInformation$Mass)){
	if(i<=length(PeptideInformation$Mass)-1){
		if ((PeptideInformation$Mass[i]-PeptideInformation$Mass[i+1])<ClusterCut){
			PeptideInformation$Cluster[i]=counter+1
			PeptideInformation$Cluster[i+1]=counter+1
		}
		else{
			PeptideInformation$Cluster[i]=counter+1
			PeptideInformation$Cluster[i+1]=counter+2
			counter=counter+1}
	}

	else
		if ((PeptideInformation$Mass[i-1]-PeptideInformation$Mass[i])<ClusterCut){
			PeptideInformation$Cluster[i]=counter+1
		}
		else{
			PeptideInformation$Cluster[i]=counter+1
		}
}}
PeptideInformation
		})
}

#############################################
##Calculate R squared of a given fit from the data
RsqCalc <- function(SpectraPeptideRange,yhat){

	r1=SpectraPeptideRange[,2]-yhat
	r2=SpectraPeptideRange[,2]-mean(SpectraPeptideRange[,2])
	R=1-(sum(r1^2)/sum(r2^2))
		
	R
}
#############################################
##Brute force search of parameter space consisting of sigma and mzError
##For looking at flat baseline
MultiPeakPeptideScan=function(FileList,PeptideInformation,SigmaRange,mzErrorRange,SpectraMin=500,SpectraMax=1350,RsqMin=0,NormalizeMethod="NoNormalization"){
##Modifying Peptide Information from PeptideInfo2 to fit this function
PeptideInformation=PeptideInformation[,-11]
colnames(PeptideInformation)[10] = "NoiseAdjust"

##Construct data frame for inforamtion on peptides looked for
TotalPeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjust=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

##Run through each of the spectra in the file, comparingthe data to the constructed PeptideInformation data frame and filling in the blanks 
Sigma=SigmaRange
mzError=mzErrorRange

for(x in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[x])
	TIC <- mean(SpectraImput@tic)
	SpectraAverage <- getSpec(SpectraImput,mzrange=c(SpectraMin,SpectraMax))

	cat("Loaded",x,"Spectra of",length(FileList),"\n")

	for(i in 1:max(PeptideInformation$Cluster)){
		Cluster <- PeptideInformation[PeptideInformation$Cluster==(i),]
		PeptideInformation[PeptideInformation$Cluster==(i),]$Spectra <- FileList[x]
		SpectraPeptideRange <- SpectraAverage[SpectraAverage[,1]>=(min(Cluster$Mass)-1.008664916) & SpectraAverage[,1]<=(max(Cluster$Mass)+5*1.008664916),]
		
		Rsq <- RsqMin
		for(j in 1:length(Sigma)){
			for(k in 1:length(mzError)){
				PredictedData <- rep(1,length(SpectraPeptideRange[,1]))  ##generate new PredictedData
				for(l in 1:length(Cluster$Mass)){
					PD <- MultiNormCalc(SpectraPeptideRange[,1],Cluster$Sequence[l],Cluster$Mass[l],Sigma[j],mzError[k])
					PredictedData=cbind(PD[,2],PredictedData)
				}
				
				Soln=qr.solve(PredictedData,SpectraPeptideRange[,2])
				yhat=PredictedData %*% Soln  ##generate fit from model ###Correct by using model...??

				R=RsqCalc(SpectraPeptideRange,yhat)
				
				if(R>Rsq){
					Rsq=R
					PeptideInformation[PeptideInformation$Cluster==(i),]$mzError=mzError[k]
					PeptideInformation[PeptideInformation$Cluster==(i),]$Sigma=Sigma[j]
						MatrixAreas=Soln[-length(Soln)]
						ReturnedAreas=rev(MatrixAreas)
					PeptideInformation[PeptideInformation$Cluster==(i),]$Area=ReturnedAreas
					PeptideInformation[PeptideInformation$Cluster==(i),]$NoiseAdjust=Soln[length(Soln)]
					PeptideInformation[PeptideInformation$Cluster==(i),]$Rsq=Rsq
					PeptideInformation[PeptideInformation$Cluster==(i),]$SpectraCount=x
						Residuals=SpectraPeptideRange[,2]-yhat
						mR=mean(Residuals)
						lR=length(Residuals)
					PeptideInformation[PeptideInformation$Cluster==(i),]$ResidualMean=mean(Residuals)
						variance=(1/(lR - 1))*(sum((Residuals - mR)^2))
					PeptideInformation[PeptideInformation$Cluster==(i),]$ResidualVariance=variance
					PeptideInformation[PeptideInformation$Cluster==(i),]$TIC=TIC
				}
			}
		}
	}
##Use Peak or Peptide for comparisons with in a spectra, Ion to compare across spectra
	switch(NormalizeMethod[1],
		Peak = {PeptideInformation$NormalizedArea=PeptideInformation$Area/max(PeptideInformation$Area)}, ##Normalizes by largest peak looked for
		Peptide = {PeptideInformation$NormalizedArea=PeptideInformation$Area/ PeptideInformation[PeptideInformation$Sequence==NormalizeMethod[2],]$Area}, ##Normalizes bu selected peptide
		Ion = {PeptideInformation$NormalizedArea=PeptideInformation$Area/TIC}, ##Normalizes by Total Ion Count/Current
		NoNormalization = {PeptideInformation$NormalizedArea=NA}) ##No normalization asked for

	TotalPeptideInformation=rbind(TotalPeptideInformation,PeptideInformation) ##Link all generated PeptideInformation data frames together
}
TotalPeptideInformation=TotalPeptideInformation[-1,]
TotalPeptideInformation
}
#############################################
##Psudo Newton-Raphson fitting method used for fitting model, looking at sigma and mzError parameters
##quandrant gradient assnet method fo crossing search space implemented
##Root bysection in two dimensions
##add. wif,.t2d,.mgf formats with text,csv,tsd
MultiPeakPeptideScan2=function(FileList,PeptideInformation,SigmaRange,mzErrorRange,RangeCutoff=c(1,5),SpectraMin=500,SpectraMax=1350,RsqMin=0,NormalizeMethod="NoNormalization",InterationLimit=100,StartPoint=c(0,.1,0)){

TotalPeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjustSlope=NA,NoiseAdjustIntercept=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

IntChoice = 1:InterationLimit %% 2 == 0

for(i in 1:length(FileList)){
	SpectraImput <- xcmsRaw(FileList[i])
	TIC <- mean(SpectraImput@tic)
	SpectraAverage <- getSpec(SpectraImput,mzrange=c(SpectraMin,SpectraMax))

	cat("Loaded",i,"Spectra of",length(FileList),"\n")
	for(j in 1:max(PeptideInformation$Cluster)){
		SavePnt=StartPoint
		Cluster <- PeptideInformation[PeptideInformation$Cluster==(j),]
		PeptideInformation[PeptideInformation$Cluster==(j),]$Spectra <- FileList[i]
		SpectraPeptideRange <- SpectraAverage[SpectraAverage[,1]>=(min(Cluster$Mass)-RangeCutoff[1]*1.008664916) & SpectraAverage[,1]<=(max(Cluster$Mass)+RangeCutoff[2]**1.008664916),]
		
		Rsq <- RsqMin
		RCheck=rep(0,InterationLimit)
		for(k in 1:InterationLimit){ 
			##Here we check a range of mzError
			if(IntChoice[k] == TRUE){
				for(l in 1:length(mzError)){ 
					PredictedData <- cbind(SpectraPeptideRange[,1],rep(1,length(SpectraPeptideRange[,1])))
					for(m in 1:length(Cluster$Mass)){  
						PD <- MultiNormCalc(SpectraPeptideRange[,1],Cluster$Sequence[m],Cluster$Mass[m],sigma=SavePnt[2],mzError=mzError[l])
						PredictedData=cbind(PD[,2],PredictedData)
					}
				
					Soln=qr.solve(PredictedData,SpectraPeptideRange[,2])
					yhat=PredictedData %*% Soln  ##generate fit from model 

					R=RsqCalc(SpectraPeptideRange,yhat)
				
					if(R>Rsq){
						Rsq=R
						PeptideInformation[PeptideInformation$Cluster==(j),]$mzError=mzError[l]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Sigma=SavePnt[2]
							MatrixAreas=Soln[-c(length(Soln)-1,length(Soln))]
							ReturnedAreas=rev(MatrixAreas)
						PeptideInformation[PeptideInformation$Cluster==(j),]$Area=ReturnedAreas
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustIntercept=Soln[length(Soln)]
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustSlope=Soln[length(Soln)-1]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Rsq=Rsq
						PeptideInformation[PeptideInformation$Cluster==(j),]$SpectraCount=i
							Residuals=SpectraPeptideRange[,2]-yhat
							mR=mean(Residuals)
							lR=length(Residuals)
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualMean=mean(Residuals)
							variance=(1/(lR - 1))*(sum((Residuals - mR)^2))
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualVariance=variance
						PeptideInformation[PeptideInformation$Cluster==(j),]$TIC=TIC
				
							SavePnt[1]=mzError[l]
							SavePnt[3]=Rsq
							RCheck[k]=Rsq
					}
					}	
				}

			##Start with mzError, here we check a range of Sigma
			if(IntChoice[k] == FALSE){
				for(n in 1:length(Sigma)){ 
					PredictedData <- cbind(SpectraPeptideRange[,1],rep(1,length(SpectraPeptideRange[,1])))
					for(o in 1:length(Cluster$Mass)){  
						PD <- MultiNormCalc(SpectraPeptideRange[,1],Cluster$Sequence[o],Cluster$Mass[o],sigma=Sigma[n],mzError=SavePnt[1])
						PredictedData=cbind(PD[,2],PredictedData)
					}
				
					Soln=qr.solve(PredictedData,SpectraPeptideRange[,2])
					yhat=PredictedData %*% Soln  ##generate fit from model 

					R=RsqCalc(SpectraPeptideRange,yhat)
				
					if(R>Rsq){
						Rsq=R
						PeptideInformation[PeptideInformation$Cluster==(j),]$mzError=SavePnt[1]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Sigma=Sigma[n]
							MatrixAreas=Soln[-c(length(Soln)-1,length(Soln))]
							ReturnedAreas=rev(MatrixAreas)
						PeptideInformation[PeptideInformation$Cluster==(j),]$Area=ReturnedAreas
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustIntercept=Soln[length(Soln)]
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustSlope=Soln[length(Soln)-1]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Rsq=Rsq
						PeptideInformation[PeptideInformation$Cluster==(j),]$SpectraCount=i
							Residuals=SpectraPeptideRange[,2]-yhat
							mR=mean(Residuals)
							lR=length(Residuals)
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualMean=mean(Residuals)
							variance=(1/(lR - 1))*(sum((Residuals - mR)^2))
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualVariance=variance
						PeptideInformation[PeptideInformation$Cluster==(j),]$TIC=TIC
							SavePnt[2]=Sigma[n]
							SavePnt[3]=Rsq
							RCheck[k]=Rsq
					}
				}
			}

			if(k>1){
				if(RCheck[k-1]==RCheck[k]){
					break
				}
			}
		}
	}
	##Use Peak or Peptide for comparisons with in a spectra, Ion to compare across spectra
	switch(NormalizeMethod[1],
		Peak = {PeptideInformation$NormalizedArea=PeptideInformation$Area/max(PeptideInformation$Area)}, ##Normalizes by largest peak looked for
		Peptide = {PeptideInformation$NormalizedArea=PeptideInformation$Area/ PeptideInformation[PeptideInformation$Sequence==NormalizeMethod[2],]$Area}, ##Normalizes bu selected peptide
		Ion = {PeptideInformation$NormalizedArea=PeptideInformation$Area/TIC}, ##Normalizes by Total Ion Count/Current
		NoNormalization = {PeptideInformation$NormalizedArea=NA}) ##No normalization asked for

	TotalPeptideInformation=rbind(TotalPeptideInformation,PeptideInformation) ##Link all generated PeptideInformation data frames together
	}
TotalPeptideInformation=TotalPeptideInformation[-1,]
TotalPeptideInformation
}
#############################################
##For removing given nth row of a data 
RowRemove <- function(Record, n){
Record[-(seq(n,to=nrow(Record),by=n)),]
}
#############################################	
##For editing results of peptide search, looking at rsq of fit and predicted area
QuickEditResults <- function(Record,RsqCutoff=.7,AreaCutoff=10){
	OutPut <- Record[Record$Rsq>=(RsqCutoff),]##Filters by R squared values
	OutPut <- OutPut[OutPut$Area>=(AreaCutoff=10),]##Filters for negative/low areas, nonsensical information in this context
	OutPut
}
#############################################
##Looks at .mzXML files
##Must set file destination, use .pdf in FileTitle and Record for fit from model, 
##will dump pdf into file containing spectra
ViewRawSpectra=function(FileTitle,FileList,PageRow=3,PageCol=2,SpectraStart=650,SpectraFinish=1350){
	pdf(file = FileTitle)
	par(mfrow = c(PageRow,PageCol))
	
	for(i in 1:length(FileList)){
		SpectraImput <- xcmsRaw(as.character(FileList[i]))
		SpectraAverage <- getSpec(SpectraImput,mzrange=c(SpectraStart,SpectraFinish))
		plot(SpectraAverage[,1],SpectraAverage[,2],type='l',xlab="MZ",ylab="Intensity",main=paste(FileList[i]))
	}
	dev.off()
	cat("Please remember to shutdown R before opening any large .PDF, thank you.","\n")
}
##For Spectra with multiple peptide clusters
##Producing NA titles??##NEEDS REFINEMENT AND CORRECDTION
ViewFitResults=function(SpectraFitRecord,FileTitle,PageRow=3,PageCol=2,SpectraStart=650,SpectraFinish=1350){
	pdf(file = FileTitle)
	par(mfrow = c(PageRow,PageCol))

	for(i in 1:max(SpectraFitRecord$SpectraCount)){
		SampleSpectra <- SpectraFitRecord[SpectraFitRecord$SpectraCount==(i),]
		SpectraImput <- xcmsRaw(as.character(SampleSpectra$Spectra[1]))

		SpectraAverage <- getSpec(SpectraImput,mzrange=c(SpectraStart,SpectraFinish))##Get spectra data

		for(j in 1:max(SampleSpectra$Cluster)){ ##NUmber of clusters is = to number of pics per spectra we want
			Cluster <- SampleSpectra[SampleSpectra$Cluster==(j),]
			SpectraPeptideRange <- SpectraAverage[SpectraAverage[,1]>=(min(Cluster$Mass)-1.008664916) & SpectraAverage[,1]<=(max(Cluster$Mass)+5*1.008664916),]	
			if(length(SpectraFitRecord[1,])<16){##for use with flat baseline, change to if there is only NoiseAdjust
				ModelPrediction <- rep(0,length(SpectraPeptideRange[,1]))
				for(k in 1:length(Cluster$Mass)){
					ModelCurve <- MultiNormCalc(SpectraPeptideRange[,1],as.character(Cluster$Sequence[k]),Cluster$Mass[k],Cluster$Sigma[k],Cluster$mzError[k])
					Curve <- (ModelCurve[,2]*Cluster$Area[k])
					ModelPrediction <- ModelPrediction+Curve
				}
				ModelPrediction <- ModelPrediction+Cluster$NoiseAdjust[j]			
			}
			else{ 
				ModelPrediction <- rep(0,length(SpectraPeptideRange[,1]))
				for(k in 1:length(Cluster$Mass)){
					ModelCurve <- MultiNormCalc(SpectraPeptideRange[,1],as.character(Cluster$Sequence[k]),Cluster$Mass[k],Cluster$Sigma[k],Cluster$mzError[k])
					Curve <- (ModelCurve[,2]*Cluster$Area[k])
					ModelPrediction <- ModelPrediction+Curve
				}
				BL=Cluster$NoiseAdjustSlope[k]*SpectraPeptideRange[,1]+Cluster$NoiseAdjustIntercept[k]
				ModelPrediction <- ModelPrediction+BL
			}
			plot(SpectraPeptideRange[,1],SpectraPeptideRange[,2],type='l',xlab="MZ",ylab="Intensity")
			title(substr(as.character(Cluster$Spectra[k]),1, 25),line=3)
			title(as.character(Cluster$Peptide[k]),line=2)			
			lines(SpectraPeptideRange[,1],ModelPrediction,col='red')
		}

	}
	dev.off()
	cat("Please remember to shutdown R before opening any large .PDF, thank you.","\n")
}
#############################################
##Looking at the first 3 peak heights (intensies) and
##calculates Riemann sums for the peaks, Must set file destination
##for .mzXML files
PeakAreaHeightScan=function(SpectraFitRecord){
HeightStore=data.frame(SpectraFile=NA,Peptide=NA,Peak1=NA,Peak2=NA,Peak3=NA,Sum=NA,AUCMinusBaseline1=NA,AUCMinusBaseline2=NA,AUCMinusBaseline3=NA,AUCTotal=NA)

	for(i in 1:length(SpectraFitRecord[,1])){
		SpectraImput <- xcmsRaw(as.character(SpectraFitRecord$Spectra[i]))
		
			cat("Loaded",i,"Spectra of",length(SpectraFitRecord[,1]),"\n")
		
		SpectraAverage <- getSpec(SpectraImput,mzrange=c(650,1350))
						
		HeightStore[i,1]=as.character(SpectraFitRecord$Spectra[i])
		HeightStore[i,2]=as.character(SpectraFitRecord$Peptide[i])
		PeakEstimate=SpectraFitRecord$Mass[i]+SpectraFitRecord$mzError[i]-1.008664916  ##Fixed twice

		if(length(SpectraFitRecord[1,])<16){ ##for use with linear baseline estimate
			for(j in 1:3){
				counter=j*1.008664916
				PeakRange=PeakEstimate+counter
				Peak=SpectraAverage[SpectraAverage[,1]>=(PeakRange-0.5043325) & SpectraAverage[,1]<=(PeakRange+0.5043325),]
			
				HeightStore[i,2+j]= max(Peak[,2]) - SpectraFitRecord$NoiseAdjust[i]
				HeightStore[i,6+j]=CurveAreaCalc(Peak,SpectraFitRecord$NoiseAdjust[i])
			}
		}
		else{
			for(j in 1:3){
				counter=j*1.008664916
				PeakRange=PeakEstimate+counter
				Peak=SpectraAverage[SpectraAverage[,1]>=(PeakRange-0.5043325) & SpectraAverage[,1]<=(PeakRange+0.5043325),]
			
				BL=SpectraFitRecord$NoiseAdjustSlope[1]*Peak[,1]+SpectraFitRecord$NoiseAdjustIntercept[1]
			
				HeightStore[i,2+j]= max(Peak[,2]) - BL[which.max(Peak[,2])]
				HeightStore[i,6+j]=CurveAreaCalc(Peak,BL)
			}
		}
		HeightStore[i,6]=sum(HeightStore[i,3:5])
		HeightStore[i,10]=sum(HeightStore[i,7:9])
	}
	HeightStore
}
############################################
##Trapazoidal AUC calc
CurveAreaCalc = function(test,baseline){
	PArea=0

	if(length(baseline)==1){
		baseline=rep(baseline,length(test[,1]))
	}
	if(length(baseline)!=length(test[,1])){
		print('CurveAreaCalc has stopped because baseline and vector with peak data are of unequal lenght, 
		\nbaseline must describe whole area under peak of be one element to describe flat baseline')
		break
	}
	for(i in 1:(length(test[,1])-1)){
		area=.5*(test[i+1,2]+test[i,2])*(test[i+1,1]-test[i,1])
		area=area-(.5*(test[i+1,1]-test[i,1])*(baseline[i+1]+baseline[i]))
		PArea=PArea+area
	}
	PArea
}


#############################################
#############################################
#############################################
#############################################
##Spectra Simulation Work
set.seed(31851)
PeakSim=function(SimCount,PickCount,Pep,mzError,Sigma,Slope,Intercept,Area,StepBinSize=.02,rho=.7,UnifRange=c(0,2),titles,SIS=0,mzR=c(1,5)){

InfoStore=list()

Mass=rep(0,length(Pep))
for(n in 1:length(Pep)){
Mass[n]=MassCalc(Pep[n])+SIS*1.008664916
}
Mass=sort(Mass)

SpectraData=data.frame(x=seq(min(Mass)-(mzR[1]*1.008664916),max(Mass)+(mzR[2]*1.008664916),StepBinSize))  
mzRange=seq(min(Mass)-(mzR[1]*1.008664916),max(Mass)+(mzR[2]*1.008664916),StepBinSize)+.5*StepBinSize

for(i in 1:SimCount){
SpectraInfo=data.frame(Simulation=NA,Peptide=NA,Mass=NA,mzError=NA,Sigma=NA,Area=NA,Slope=NA,Intercep=NA,PeakCount=NA,RemovedCount=NA)
SpectraInfo[1:length(Mass),]=0

AreaCurrent=list()
for(m in 1:length(Pep)){
	lambda=IsotopicCountFFT(Pep[m])
	StorePeakDraws=list()
	Peak=rmultinom(1,PickCount,lambda)
		for(j in 1:length(lambda)){
			PeakDraws=rnorm(Peak[j],Mass[m]+((j-1)*1.008664916)+mzError[m],Sigma[m])
			StorePeakDraws[[j]]=PeakDraws
}
Points=unlist(StorePeakDraws)

CurrentCount=sum(table(cut(Points,breaks=seq(min(Mass)-(mzR[1]*1.008664916),max(Mass)+(mzR[2]*1.008664916),StepBinSize))))

RawCurrent=as.vector(table(cut(Points,breaks=seq(min(Mass)-(mzR[1]*1.008664916),max(Mass)+(mzR[2]*1.008664916)+StepBinSize,StepBinSize))))
NormalizedCurrent=(RawCurrent/CurrentCount)*(1/StepBinSize)

AreaCurrent[[m]]=NormalizedCurrent*Area[m]
}
TotalAreaCurrent=Reduce('+',AreaCurrent)

##AR(1) Noise
if(rho>0){
Base=rep(10,(length(mzRange)+200))
for(k in 1:(length(Base)-1)){
	Base[k+1]=rho*Base[k]+runif(1,UnifRange[1],UnifRange[2])
}
AR1Noise=Base[101:(length(Base)-100)]-mean(Base[101:(length(Base)-100)])
##adjust to mean of noise to 0

##Combine AR(1) and peaks
for(l in 1:length(TotalAreaCurrent)){
if(TotalAreaCurrent[l]==0){
	TotalAreaCurrent[l]=AR1Noise[l] ##remove l+100
}
if(TotalAreaCurrent[l]<AR1Noise[l]){
	TotalAreaCurrent[l]=AR1Noise[l]
}
}
}

ModelCurrent = TotalAreaCurrent

if(abs(Slope)>0){
##ABS() since slope can be +/-
Lift = Slope*mzRange+((min(Mass)*-1*Slope)+Intercept)
ModelCurrent = TotalAreaCurrent+Lift
}

SpectraData[2*i-1] <- mzRange 
SpectraData[2*i] <- ModelCurrent

SpectraInfo[,1]=i
SpectraInfo[,2]=Pep
SpectraInfo[,3]=Mass
SpectraInfo[,4]=mzError
SpectraInfo[,5]=Sigma
SpectraInfo[,6]=Area
SpectraInfo[,7]=Slope
SpectraInfo[,8]=((Mass*-1*Slope)+Intercept)
SpectraInfo[,9]=PickCount
SpectraInfo[,10]=PickCount-CurrentCount

InfoStore[[i]]=SpectraInfo
}

SumInfo=data.frame(Simulation=NA,Peptide=NA,Mass=NA,mzError=NA,Sigma=NA,Area=NA,Slope=NA,Intercep=NA,PeakCount=NA,RemovedCount=NA)
for(i in 1:length(InfoStore)){
SumInfo=rbind(SumInfo,InfoStore[[i]])
}
SumInfo=SumInfo[-1,]

write.csv(SumInfo,titles[1])
write.csv(SpectraData,titles[2])
}
##For csv files
MultiPeakPeptideScan3=function(FileList,PeptideInformation,SigmaRange,mzErrorRange,RangeCutoff=c(1,5),RsqMin=0,NormalizeMethod="NoNormalization",InterationLimit=100,StartPoint=c(0,.1,0)){

TotalPeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjustSlope=NA,NoiseAdjustIntercept=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

IntChoice = 1:InterationLimit %% 2 == 0

SpectraData <- read.csv(FileList[1])
SpectraNumber=length(SpectraData[1,])/2

if(SpectraNumber <1){
SpectraData <- read.csv(FileList,sep='\t')
print('File in tab delimited form')
}
if(SpectraNumber %% 2 == .5){
SpectraData <- SpectraData[,-1]
SpectraNumber=length(SpectraData[1,])/2
print('File has uneven number of columns, removing first column')
}

for(i in 1:SpectraNumber){

	cat("Loaded",i,"Spectra of",SpectraNumber,"\n")
	TIC <- mean(SpectraData[,2*i])
	SpectraAverage <- cbind(SpectraData[,2*i-1],SpectraData[,2*i])

	for(j in 1:max(PeptideInformation$Cluster)){
		SavePnt=StartPoint
		Cluster <- PeptideInformation[PeptideInformation$Cluster==(j),]
		PeptideInformation[PeptideInformation$Cluster==(j),]$Spectra <- FileList
		SpectraPeptideRange <- SpectraAverage[SpectraAverage[,1]>=(min(Cluster$Mass)-RangeCutoff[1]*1.008664916) & SpectraAverage[,1]<=(max(Cluster$Mass)+RangeCutoff[2]**1.008664916),]
		
		Rsq <- RsqMin
		RCheck=rep(0,InterationLimit)
		for(k in 1:InterationLimit){ 
			##Here we check a range of mzError
			if(IntChoice[k] == TRUE){
				for(l in 1:length(mzError)){ 
					PredictedData <- cbind(SpectraPeptideRange[,1],rep(1,length(SpectraPeptideRange[,1])))
					for(m in 1:length(Cluster$Mass)){  
						PD <- MultiNormCalc(SpectraPeptideRange[,1],Cluster$Sequence[m],Cluster$Mass[m],sigma=SavePnt[2],mzError=mzError[l])
						PredictedData=cbind(PD[,2],PredictedData)
					}
				
					Soln=qr.solve(PredictedData,SpectraPeptideRange[,2])
					yhat=PredictedData %*% Soln  ##generate fit from model 

					R=RsqCalc(SpectraPeptideRange,yhat)
				
					if(R>Rsq){
						Rsq=R
						PeptideInformation[PeptideInformation$Cluster==(j),]$mzError=mzError[l]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Sigma=SavePnt[2]
							MatrixAreas=Soln[-c(length(Soln)-1,length(Soln))]
							ReturnedAreas=rev(MatrixAreas)
						PeptideInformation[PeptideInformation$Cluster==(j),]$Area=ReturnedAreas
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustIntercept=Soln[length(Soln)]
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustSlope=Soln[length(Soln)-1]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Rsq=Rsq
						PeptideInformation[PeptideInformation$Cluster==(j),]$SpectraCount=i
							Residuals=SpectraPeptideRange[,2]-yhat
							mR=mean(Residuals)
							lR=length(Residuals)
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualMean=mean(Residuals)
							variance=(1/(lR - 1))*(sum((Residuals - mR)^2))
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualVariance=variance
						PeptideInformation[PeptideInformation$Cluster==(j),]$TIC=TIC
				
							SavePnt[1]=mzError[l]
							SavePnt[3]=Rsq
							RCheck[k]=Rsq
					}
					}	
				}

			##Start with mzError, here we check a range of Sigma
			if(IntChoice[k] == FALSE){
				for(n in 1:length(Sigma)){ 
					PredictedData <- cbind(SpectraPeptideRange[,1],rep(1,length(SpectraPeptideRange[,1])))
					for(o in 1:length(Cluster$Mass)){  
						PD <- MultiNormCalc(SpectraPeptideRange[,1],Cluster$Sequence[o],Cluster$Mass[o],sigma=Sigma[n],mzError=SavePnt[1])
						PredictedData=cbind(PD[,2],PredictedData)
					}
				
					Soln=qr.solve(PredictedData,SpectraPeptideRange[,2])
					yhat=PredictedData %*% Soln  ##generate fit from model 

					R=RsqCalc(SpectraPeptideRange,yhat)
				
					if(R>Rsq){
						Rsq=R
						PeptideInformation[PeptideInformation$Cluster==(j),]$mzError=SavePnt[1]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Sigma=Sigma[n]
							MatrixAreas=Soln[-c(length(Soln)-1,length(Soln))]
							ReturnedAreas=rev(MatrixAreas)
						PeptideInformation[PeptideInformation$Cluster==(j),]$Area=ReturnedAreas
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustIntercept=Soln[length(Soln)]
						PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjustSlope=Soln[length(Soln)-1]
						PeptideInformation[PeptideInformation$Cluster==(j),]$Rsq=Rsq
						PeptideInformation[PeptideInformation$Cluster==(j),]$SpectraCount=i
							Residuals=SpectraPeptideRange[,2]-yhat
							mR=mean(Residuals)
							lR=length(Residuals)
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualMean=mean(Residuals)
							variance=(1/(lR - 1))*(sum((Residuals - mR)^2))
						PeptideInformation[PeptideInformation$Cluster==(j),]$ResidualVariance=variance
						PeptideInformation[PeptideInformation$Cluster==(j),]$TIC=TIC
							SavePnt[2]=Sigma[n]
							SavePnt[3]=Rsq
							RCheck[k]=Rsq
					}
				}
			}

			if(k>1){
				if(RCheck[k-1]==RCheck[k]){
					break
				}
			}
		}
	}
	##Use Peak or Peptide for comparisons with in a spectra, Ion to compare across spectra
	switch(NormalizeMethod[1],
		Peak = {PeptideInformation$NormalizedArea=PeptideInformation$Area/max(PeptideInformation$Area)}, ##Normalizes by largest peak looked for
		Peptide = {PeptideInformation$NormalizedArea=PeptideInformation$Area/ PeptideInformation[PeptideInformation$Sequence==NormalizeMethod[2],]$Area}, ##Normalizes bu selected peptide
		Ion = {PeptideInformation$NormalizedArea=PeptideInformation$Area/TIC}, ##Normalizes by Total Ion Count/Current
		NoNormalization = {PeptideInformation$NormalizedArea=NA}) ##No normalization asked for

	TotalPeptideInformation=rbind(TotalPeptideInformation,PeptideInformation) ##Link all generated PeptideInformation data frames together
	}
TotalPeptideInformation=TotalPeptideInformation[-1,]
TotalPeptideInformation
}
##Random Peptide Generation and Filtering
