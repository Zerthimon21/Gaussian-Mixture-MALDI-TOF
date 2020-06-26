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
##Just one of either peptide?? gives error for one peptide
PeptideInfo=function(UnlabeledPeptides,UnlabeledSequence,LabeledPeptides,LabeledSequence,AQUA=FALSE,Neutron=0,ClusterCut=5){

	Peptides=cbind(UnlabeledSequence)
	rownames(Peptides)=UnlabeledPeptides
	HeavyPeptides=cbind(LabeledSequence)
	rownames(HeavyPeptides)=LabeledPeptides
	
	PeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjust=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

for(i in 1:length(Peptides)){ 
	PeptideInformation[i,1] <- rownames(Peptides)[i]
	PeptideInformation[i,2] <- Peptides[i]
	PeptideInformation[i,3] <- MassCalc(Peptides[i]) 
	if(AQUA==TRUE){
		for(i in 1:length(HeavyPeptides)){
			PeptideInformation[i+length(Peptides),1] <- paste('AQUA',rownames(HeavyPeptides)[i])
			PeptideInformation[i+length(Peptides),2] <- HeavyPeptides[i]
			PeptideInformation[i+length(Peptides),3] <- MassCalc(HeavyPeptides[i])+Neutron*1.008664916
		}
	}
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
}

##Check if peptide or heavypeptide have anything in them
##return a error otherwise
##if any names empty say unknown #?? fill in blanks
##if unlabeled() DONE
##if labeled() DONE
#order by mass DONE
#assign cluster DONE
#return pepinfo data frame DONE

##Updated for either/or peptides and single pepides
PeptideInfo2=function(UnlabeledPeptides,UnlabeledSequence,LabeledPeptides,LabeledSequence,AQUAState="Both",Neutron=0,ClusterCut=5){
switch(AQUAState,
	Both = {##Makes list for both labeled and unlabeled peptides
	Peptides=cbind(UnlabeledSequence)
	rownames(Peptides)=UnlabeledPeptides
	HeavyPeptides=cbind(LabeledSequence)
	rownames(HeavyPeptides)=LabeledPeptides
	
	PeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjust=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

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
	
	PeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjust=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)


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
	
	PeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjust=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

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
MultiPeakPeptideScan=function(FileList,PeptideInformation,SigmaRange,mzErrorRange,SpectraMin=500,SpectraMax=1350,RsqMin=0,NormalizeMethod="NoNormalization"){

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
##Root bysection in two dimensions
##add. wif,.t2d,.mgf formats with text,csv,tsd
MultiPeakPeptideScan2=function(FileList,PeptideInformation,SigmaRange,mzErrorRange,SpectraMin=500,SpectraMax=1350,RsqMin=0,NormalizeMethod="NoNormalization",InterationLimit=10,StartPoint=c(0,.1,0)){

TotalPeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjust=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA,TIC=NA,NormalizedArea=NA)

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
		SpectraPeptideRange <- SpectraAverage[SpectraAverage[,1]>=(min(Cluster$Mass)-1.008664916) & SpectraAverage[,1]<=(max(Cluster$Mass)+5*1.008664916),]
		
		Rsq <- RsqMin
		RCheck=rep(0,InterationLimit)
		for(k in 1:InterationLimit){ 
			##Here we check a range of mzError
			if(IntChoice[k] == TRUE){
				for(l in 1:length(mzError)){ 
				PredictedData <- rep(1,length(SpectraPeptideRange[,1]))  ##generate new PredictedData
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
						MatrixAreas=Soln[-length(Soln)]
						ReturnedAreas=rev(MatrixAreas)
					PeptideInformation[PeptideInformation$Cluster==(j),]$Area=ReturnedAreas
					PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjust=Soln[length(Soln)]
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
					PredictedData <- rep(1,length(SpectraPeptideRange[,1]))  ##generate new PredictedData
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
						MatrixAreas=Soln[-length(Soln)]
						ReturnedAreas=rev(MatrixAreas)
					PeptideInformation[PeptideInformation$Cluster==(j),]$Area=ReturnedAreas
					PeptideInformation[PeptideInformation$Cluster==(j),]$NoiseAdjust=Soln[length(Soln)]
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
##Producing NA titles??
ViewFitResults=function(SpectraFitRecord,FileTitle,PageRow=3,PageCol=2,SpectraStart=650,SpectraFinish=1350){
	pdf(file = FileTitle)
	par(mfrow = c(PageRow,PageCol))
	
	for(i in 1:max(SpectraFitRecord$SpectraCount)){
		SampleSpectra <- SpectraFitRecord[SpectraFitRecord$SpectraCount==(i),]
		SpectraImput <- xcmsRaw(as.character(SampleSpectra[1,4]))
		SpectraAverage <- getSpec(SpectraImput,mzrange=c(SpectraStart,SpectraFinish))##set check for lenght fo spectra
		for(j in 1:max(SampleSpectra$Cluster)){
			Cluster <- SampleSpectra[SampleSpectra$Cluster==(j),]
			SpectraPeptideRange <- SpectraAverage[SpectraAverage[,1]>=(min(Cluster$Mass)-1.008664916) & SpectraAverage[,1]<=(max(Cluster$Mass)+5*1.008664916),]	
			plot(SpectraPeptideRange[,1],SpectraPeptideRange[,2],type='l',xlab="MZ",ylab="Intensity",main=paste(SpectraFitRecord[j,4],SpectraFitRecord[j,1]))
			ModelPrediction <- rep(0,length(SpectraPeptideRange[,1]))
			for(k in 1:length(Cluster$Mass)){
				ModelCurve <- MultiNormCalc(SpectraPeptideRange[,1],Cluster$Sequence[k],Cluster$Mass[k],Cluster$Sigma[k],Cluster$mzError[k])
				Curve <- (ModelCurve[,2]*Cluster$Area[k])
				ModelPrediction <- ModelPrediction+Curve
			}
			ModelPrediction <- ModelPrediction+Cluster$NoiseAdjust[1]
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

##change all colum references to column title references
PeakAreaHeightScan=function(SpectraFitRecord){
HeightStore=data.frame(SpectraFile=NA,Peptide=NA,Peak1=NA,Peak2=NA,Peak3=NA,Sum=NA,AUCMinusBaseline1=NA,AUCMinusBaseline2=NA,AUCMinusBaseline3=NA,AUCTotal=NA)

	for(i in 1:length(SpectraFitRecord[,1])){
		SpectraImput <- xcmsRaw(SpectraFitRecord[i,4])
		SpectraAverage <- getSpec(SpectraImput,mzrange=c(650,1350))
						
		HeightStore[i,1]=SpectraFitRecord[i,4]
		HeightStore[i,2]=SpectraFitRecord[i,1]
		PeakEstimate=SpectraFitRecord[i,3]+SpectraFitRecord[i,7]-1.008664916  ##Fixed twice
		
		for(j in 1:3){
			counter=j*1.008664916
			PeakRange=PeakEstimate+counter
			Peak=SpectraAverage[SpectraAverage[,1]>=(PeakRange-0.5043325) & SpectraAverage[,1]<=(PeakRange+0.5043325),]
			
			HeightStore[i,2+j]= max(Peak[,2]) - SpectraFitRecord[i,10]
			HeightStore[i,6+j]=CurveAreaCalc(Peak,SpectraFitRecord[i,10])
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
for(i in 1:(length(test[,1])-1)){
	area=.5*(test[i+1,2]+test[i,2])*(test[i+1,1]-test[i,1])
	area=area-((test[i+1,1]-test[i,1])*baseline)
	PArea=PArea+area
}
PArea
}


#############################################
#############################################
#############################################
#############################################
##Spaectra Simulation Work
#############################################
##Simulates spectra from a given set of peptides
SpectraSimulation=function(SimNumber,SpectraRange,PeakIntensity,SamplePeptides,PeakWidth,Noise=c(1,1,1,1),AquaPep=NA,AquaInt=NA,N=6){

#Generate data storage for output
SpectraInfo=data.frame(Simulation=NA,Peptide=NA,TrueMass=NA,EstMass=NA,mzError=NA,Cluster=NA,TrueArea=NA,TrueSigma=NA,NoiseDist=NA,BaselineShift=NA)
SpectraInfo=SpectraInfo[-1,]
SampledSpectra=matrix(SpectraRange,length(SpectraRange),1)

for(x in 1:SimNumber){
print(x)
#data storage for each spectra we look at
PepInfo=data.frame(Simulation=NA,Peptide=NA,TrueMass=NA,EstMass=NA,mzError=NA,Cluster=NA,TrueArea=NA,TrueSigma=NA,NoiseDist=NA,BaselineShift=NA)
PepInfo=PepInfo[-1,]

##Generate a list of peptides to look at for a given Spectra Simulation
NumPeaks=ceiling(runif(1,0,(length(SamplePeptides))))
SelectPeaks=SamplePeptides
SpectraSample="Null"

##Select the peptides for each spectra at random
if(Noise[1]>0){
	for(a in 1:NumPeaks){
		Choose=ceiling(runif(1,0,(length(SelectPeaks))))
		SpectraSample[a]=SelectPeaks[Choose]
		SelectPeaks=SelectPeaks[-Choose]
	}
}
else{
SpectraSample=SamplePeptides
}

##Calc Mass of peptides
SpectraSampleMass=data.frame(Simulation=rep(x,length(SpectraSample)),Peptide=SpectraSample,Mass=rep(0,length(SpectraSample)),Cluster=rep(0,length(SpectraSample)),Intensity=rep(0,length(SpectraSample)))
for(q in 1:length(SpectraSample)){
	SpectraSampleMass$Mass[q]=MassCalc(SpectraSample[q])
	IntIntensity=ceiling(runif(1,0,(length(PeakIntensity))))
	SpectraSampleMass$Intensity[q]=PeakIntensity[IntIntensity]
}

##Add Aqua Peptides if wanted
if(is.na(AquaPep[1])==FALSE){	
	AquaMass=data.frame(Simulation=rep(x,length(AquaPep)),Peptide=AquaPep,Mass=rep(0,length(AquaPep)),Cluster=rep(0,length(AquaPep)),Intensity=AquaInt)	
	for (r in 1:length(AquaPep)){
		AquaMass$Mass[r]=MassCalc(AquaPep[r])+(N*1.008664916)
	}
	SpectraSampleMass=rbind(AquaMass,SpectraSampleMass)
}


##Order peptides based on mass
SpectraSampleMass <- SpectraSampleMass[order(-SpectraSampleMass$Mass),]

##Assign clusters to each group of peaks less than 7 mz apart
counter=0
for(i in 1:length(SpectraSampleMass$Mass)){
	if(i<=length(SpectraSampleMass$Mass)-1){
		if ((SpectraSampleMass$Mass[i]-SpectraSampleMass$Mass[i+1])<7){	##Change from 5 to make sure labeled peaks have same sigma as unlabeled, hardcoded here
			SpectraSampleMass$Cluster[i]=counter+1
			SpectraSampleMass$Cluster[i+1]=counter+1
		}
		else{
			SpectraSampleMass$Cluster[i]=counter+1
			SpectraSampleMass$Cluster[i+1]=counter+2
			counter=counter+1}
	}
	else
		if ((SpectraSampleMass$Mass[i-1]-SpectraSampleMass$Mass[i])<7){##5->7
			SpectraSampleMass$Cluster[i]=counter+1
		}
		else{
			SpectraSampleMass$Cluster[i]=counter+1
		}
}

##keeping strings as characters for later use
SpectraSampleMass[,2] = as.character(SpectraSampleMass[,2])
##Record keeping
p=length(SpectraSampleMass[,1])
PepInfo[1:p,1]=SpectraSampleMass$Simulation
PepInfo[1:p,2]=SpectraSampleMass$Peptide
PepInfo[1:p,3]=SpectraSampleMass$Mass
PepInfo[1:p,6]=SpectraSampleMass$Cluster
PepInfo[1:p,7]=SpectraSampleMass$Intensity


##generate the Combinded spectra for the randomly selected peptides
SimStore=rep(0,length(SpectraRange))

##Generate the Spectra for each peak/cluster and paste them together
InfoCount=0
for(i in 1:max(SpectraSampleMass$Cluster)){
	PeptideGroup = SpectraSampleMass[SpectraSampleMass$Cluster==(i),]

	PeptideGroup[,2] = as.character(PeptideGroup[,2])

	##Generate rnd var
	IntSigma=ceiling(runif(1,0,(length(PeakWidth))))	

	for(j in 1:length(PeptideGroup[,1])){
		InfoCount=InfoCount+1
	
		PepMass=PeptideGroup[j,3]
		if(Noise[2]>0){##mzError based on evidence
		SpectraMass=(5e-7)*(PepMass)^2+.9993*(PepMass)
		}
		else{##can be exchanged with a rnd,unf,-.5,.5
		SpectraMass=PepMass-round(runif(1,-.5,.5),2)
		}
		mzError=SpectraMass-PepMass	
		
		PepInfo[InfoCount,4]=SpectraMass
		PepInfo[InfoCount,5]=mzError
		PepInfo[InfoCount,8]=PeakWidth[IntSigma]

		##Calc Peak with adjustments
		PeakSample=MultiNormCalc(SpectraRange,PeptideGroup[j,2],PeptideGroup[j,3],PeakWidth[IntSigma],mzError)
		##Give AUC of Peak
		PeakSample[,2]=PeakSample[,2]*PeptideGroup$Intensity[j]

		SimStore=PeakSample[,2]+SimStore
		}
}

##Add error to generated spectra
SimStoreError=SimStore
if(Noise[3]>0){##per point error range
	SelectDir=rbinom(1,1,0.4288417)
	if(SelectDir>0){##small error range
		LogVar=rnorm(1,0.3620344,1.237027)
		Var=exp(LogVar)
		sdev=sqrt(Var)
		name='dist1'
		}
	else{##large error range
		LogVar=rnorm(1,5.3593586,1.547924)
		Var=exp(LogVar)
		sdev=sqrt(Var)
		name='dist2'
		}
	PepInfo[1:p,9]=name
}
else{
	PepInfo[1:p,9]=NA
	sdev=0
	}

if(Noise[4]>0){##baseline lift
	draw=rnorm(1,3.279316,0.6534015)##fitting a normal to the dirtbution of the log of the baselineAdj	
	baseline=exp(draw)
	PepInfo[1:p,10]=baseline
}
else{
PepInfo[1:p,10]=0
baseline=0
}

for(i in 1:length(SimStore)){##adding of error
		adj=rnorm(1,0,sdev)
		SimStoreError[i]=SimStore[i]+adj+baseline
		if(SimStoreError[i]<0){SimStoreError[i]=0}
}

SpectraInfo=rbind(SpectraInfo,PepInfo)
SampledSpectra=cbind(SampledSpectra,SimStoreError)

}
Output=list(SpectraInfo,SampledSpectra)
Output
}
#############################################
##Find the three highest peak intensities and AUC in simulated spectra
PeakHeightSimulation=function(SData,spectra,peaks=3){
FullStore=data.frame(Simulation=NA,EstMass=NA,Peak1=NA,Peak2=NA,Peak3=NA,Sum=NA,AUCMinusBaseline1=NA,AUCMinusBaseline2=NA,AUCMinusBaseline3=NA,AUCTotal=NA)
	##redo naming convention to reflect mutable number of peaks
	for (i in 1:max(SData$Simulation)){
	HeightStore=data.frame(Simulation=NA,EstMass=NA,Peak1=NA,Peak2=NA,Peak3=NA,Sum=NA,AUCMinusBaseline1=NA,AUCMinusBaseline2=NA,AUCMinusBaseline3=NA,AUCTotal=NA)
	SpectraData = SData[SData$Simulation==(i),]
		for(j in 1:length(SpectraData$EstMass)){	
		HeightStore[j,1]=i
		HeightStore[j,2]=SpectraData$EstMass[j]
		height=rep(0,length(peaks))
		AUCnoBaseline=rep(0,length(peaks))
		for(k in 1:peaks){
			adj=(k-1)*1.008664916
			point=round((SpectraData$EstMass[j]-spectra[1,1]+(adj))*100,0)
			a=point-50
			b=point+50
			height[k]=max(spectra[a:b,1+i])
			cut=cbind(spectra[a:b,1],spectra[a:b,1+i])
			AUCnoBaseline[k]=CurveAreaCalc(cut,SData$BaselineShift[j])
		}
		HeightStore[j,3]=height[1]
		HeightStore[j,4]=height[2]
		HeightStore[j,5]=height[3]
		HeightStore[j,6]=sum(height)
		HeightStore[j,7]=AUCnoBaseline[1]
		HeightStore[j,8]=AUCnoBaseline[2]
		HeightStore[j,9]=AUCnoBaseline[3]
		HeightStore[j,10]=sum(AUCnoBaseline)
	}
	FullStore=rbind(FullStore,HeightStore)

	}
FullStore=FullStore[-1,]	
FullStore
}
#############################################
##Look at generated data with my method,
SpectrumScan=function(spectra,PeptideInformation,SigmaRange,mzErrorRange,ClusterCut=5,RsqMin=0){

##Construct data frame for inforamtion on peptides looked for
TotalPeptideInformation=data.frame(Peptide=NA,Sequence=NA,Mass=NA,Spectra=NA,SpectraCount=NA,Cluster=NA,mzError=NA,Sigma=NA,Area=NA,NoiseAdjust=NA,Rsq=NA,ResidualMean=NA,ResidualVariance=NA)

##Run through each of the spectra in the file, comparingthe data to the constructed PeptideInformation data frame and filling in the blanks 
Sigma=SigmaRange
mzError=mzErrorRange

for(x in 1:(length(spectra[1,])-1)){
	SpectraAverage <- cbind(spectra[,1],spectra[,x+1])

	cat("Loaded",x,"Spectra of",length(spectra[1,])-1,"\n")

	for(i in 1:max(PeptideInformation$Cluster)){
		Cluster <- PeptideInformation[PeptideInformation$Cluster==(i),]
		PeptideInformation[PeptideInformation$Cluster==(i),]$Spectra <- x
		SpectraPeptideRange <- SpectraAverage[SpectraAverage[,1]>=(min(Cluster$Mass)-1.008664916) & SpectraAverage[,1]<=(max(Cluster$Mass)+ClusterCut*1.008664916),]##changed 5->7
		
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
					PeptideInformation[PeptideInformation$Cluster==(i),]$Area=rev(Soln[-length(Soln)])
					PeptideInformation[PeptideInformation$Cluster==(i),]$NoiseAdjust=Soln[length(Soln)]
					PeptideInformation[PeptideInformation$Cluster==(i),]$Rsq=Rsq
					PeptideInformation[PeptideInformation$Cluster==(i),]$SpectraCount=x
						Residuals=SpectraPeptideRange[,2]-yhat
						mR=mean(Residuals)
						lR=length(Residuals)
					PeptideInformation[PeptideInformation$Cluster==(i),]$ResidualMean=mean(Residuals)
						variance=(1/(lR - 1))*(sum((Residuals - mR)^2))
					PeptideInformation[PeptideInformation$Cluster==(i),]$ResidualVariance=variance
					}
			}
		}
	}
	##Link all generated PeptideInformation data frames together...
	TotalPeptideInformation=rbind(TotalPeptideInformation,PeptideInformation)
}
TotalPeptideInformation=TotalPeptideInformation[-1,]
TotalPeptideInformation
}
#############################################

