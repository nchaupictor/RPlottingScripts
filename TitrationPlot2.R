#library(ggplot2)

#Run this block first - prompts for number of slides and filenames
rm(list=ls(all=TRUE))

#Set working directory
FILE <- file.choose()
DIR  <- dirname(FILE)
setwd(DIR)

n <- -1
while(n < 0){ 
 	 n <- readline("Enter number of slides: ")
	 test <- readline("Enter type of slide (1 = Biotin/Strep, 2 = Other):")
	 test <- as.integer(test)
 	 n <- as.integer(n)
	 tempF <- mat.or.vec(n,1)
	 temp <- mat.or.vec(n,1)

 	 if (is.na(n)){
  		  n <- readline("Enter number of slides: ")
  	}
 	 for(i in 1:n){
		tempF[i] <- readline("Enter file name (excluding .csv): ")
	      temp[i] <- paste(tempF[i] , ".csv",sep ="")
	}
}



#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

#Run the rest of this code 
global <- mat.or.vec(320,n)

for(m in 1:n){

	data <- read.csv(temp[m],header=TRUE)
	data1 <- as.matrix(data)

	count <- 0
	count2 <- 0
	count3 <- 1

	# Convert raw data into a single 1x320 vector sorting from Spot 1, Spot 6, Spot 2, Spot 7, etc
	a <- mat.or.vec(320, 1)

	for(i in 1:16){
		for(j in 1:10){
			for(k in 1:2){
				if(j <= 5){
					a[count3] <- data1[j+count,k+count2]
					count3 <- count3 + 1
				}
				else{
					a[count3] <- data1[j+count-5,k+2+count2]
					count3 <- count3 + 1
				}
			}
		}
			count <- count+7
			if(i == 8){
				count <- 0
				count2 <- count2 + 6
			}

		if(n ==1){
			global <- a
		}else{
			global[,m] <- a
		}
	}
}
if(test ==2){
	aMat = array(0, dim = c(16,2,n))
for(m in 1:n){
	
	#Grab 4th and 13th position spots 
	for(i in 1:16){
		for(j in 1:2){
			if(j == 1){
			aMat[i,j,m] <- global[4+(20*(i-1)),m]
			}else{
			aMat[i,j,m] <- global[13+(20*(i-1)),m]
			}
		}
	}
}
}



#--------------------------------------------------------------------------
for(l in 1:n){

	#Create 8x2 plot 
	#fileP <- paste(tempF[l],"Raw Intensity.pdf")
	#pdf(fileP)
	x11()
	par(mfrow=c(8,2))
	par(mar=c(1,2,1,1))
	if(test ==1){
	for(i in 1:16){
		plot(as.numeric(global[1+(20*(i-1)+1):(20*i)-1,l]),ylim=c(0,4e6))
		title1 <- paste(temp[l],"Raw Intensity")
		mtext(title1, side = 3, line = -1.1, outer = TRUE)
		mtext(title1, side = 3, line = -1.1, outer = TRUE)
			
		}
	}else{
	for(i in 1:16){
		plot(as.numeric(aMat[i,,l]), xaxt="n",col = "red",ylim =c(0,2e6))
		title1 <- paste(temp[l], "Raw Intensity")
		mtext(title1, side = 3, line = -1.1, outer = TRUE)
		mtext(title1, side = 3, line = -1.1, outer = TRUE)
		axis(1, xaxp=c(1, 2, 1), las=2)
		mtext(paste("Well",i),side=1,line=0.25)

	}
	}
	
	#dev.off();	
}
#--------------------------------------------------------------------------
#Calculate Mean Intensity /STD per Slide
if(test ==1){
#Split vector into 16x20
#aMat <- mat.or.vec(16,20)
aMat = array(0, dim = c(16,20,n))

for(l in 1:n){
	count <- 1
	for(i in 1:16){
		for(j in 1:20){
			aMat[i,j,l] <- global[count,l]
			count <- count + 1
		}
	}
}

#Average columns
aMean <- array(0,dim=c(1,20,n))
aSd <- array(0,dim=c(1,20,n))
aCv <- array(0,dim=c(1,20,n))

for(l in 1:n){
	for(j in 1:20){
		aMean[1,j,l] <- mean(as.numeric(aMat[,j,l]))
		aSd[1,j,l] <- sd(as.numeric(aMat[,j,l]))
	}
}
#aCv[1,j,l] <- aSd / aMean * 100
aCv <- aSd / aMean * 100
}else{
aMean <- array(0,dim=c(16,1,n))
aSd <- array(0,dim=c(16,1,n))
aCv <- array(0,dim=c(16,1,n))

for(l in 1:n){
	for(i in 1:16){
		aMean[i,1,l] <- mean(as.numeric(aMat[i,,l]))
		aSd[i,1,l] <- sd(as.numeric(aMat[i,,l]))
	}
}
aCv <- aSd / aMean * 100
}

#--------------------------------------------------------------------------
#Plot Intra-slide Intensity and CV across Wells of same slide on two axis
for(l in 1:n){
	x11() 
	par(mfrow=c(1,1))
	par(mar=c(4,4,2,4))

	plot(as.numeric(aMean[,,l]),col="red",ylab="",xlab="",axes=FALSE, yaxt="n")
	if(test ==1){
	axis(2,ylim=c(0,4e6),col="red",col.axis="red",las=1, at = c(0,1e6,2e6,3e6,4e6))
	}else{
	axis(2,ylim=c(0,2e6),col="red",col.axis="red",las=1, at = c(0,1e6,2e6))
	}
	mtext(side = 2, line = 3, 'Intensity')

	par(new=TRUE)
	if(test ==1){
	plot(as.numeric(aCv[,,l]),col="blue",ylab="",xlab="Spots",axes=FALSE, yaxt="n")
	}else{
	plot(as.numeric(aCv[,,l]),col="blue",ylab="",xlab="Wells",axes=FALSE, yaxt="n")
	}
	axis(1, xaxp=c(1, 20, 19), las=2)
	axis(4,ylim=c(0,100),col="blue",col.axis="blue",las=1,at = c(0,20,40,60,80,100))

	mtext(side = 4, line = 3, 'CV (%)')
	if(test ==1){
	title2 <- paste(temp[l],"Raw Intensity vs CV")
	}else{
	title2 <- paste(temp[l],"Average Duplicate Intensity vs CV")
	}
	mtext(title2, side = 3, line = -1.1, outer = TRUE)
	mtext(title2, side = 3, line = -1.1, outer = TRUE)
	box()

	#dev.copy(jpeg,filename="IntensityvsCV.jpg");
	#dev.off ();
}




#--------------------------------------------------------------------------
#Plot of inter-slide spot concentrations 
if(test ==1){
spotAv <- mat.or.vec(10,n)
spotAvCv <- mat.or.vec(10,n)
aMatCombine <- array(0, dim = c(32,10,n))

for (k in 1:n){
	count <- 1
	for (i in seq(2,20,by=2)){
		for (j in 1:2){
			aMatCombine[1+(16*(j-1)):(16*j-1),i/2,k] <- aMat[,count,k]
			count <- count + 1
		}
	}
}


for (j in 1:4){
	count <- 1
	for(i in seq(2,20,by=2)){
		spotAv[count,j] <- mean(as.numeric(aMean[1,(i-1):i,j]))
		spotAvCv[count,j] <- mean(as.numeric(aCv[1,(i-1):i,j]))
		count <- count + 1
	}
}

x11()
par(mfrow=c(5,2))
par(mar=c(2,2,1,1))

for (i in seq(2,20,by=2)){
	for(j in 1:n){
		#boxplot(x = as.list(as.numeric(spotAv[i,])))
		boxplot(as.numeric(aMat[,(i-1):i,]))
		axis(side = 1, xlab="", xaxt="n", axes = FALSE)
	}
}
mtext("Mean Inter-Slide Spot Intensity", side = 3, line = -1.1, outer = TRUE)

x11()
par(mfrow=c(5,2))
par(mar=c(2,2,1,1))
#Plot of inter-slide spot CV
for (i in 1:10){
	plot(as.numeric(spotAvCv[i,]),xaxt="n",col = "red")
	#axis(side = 1, xlab="", xaxt="n", axes = FALSE)
	axis(side = 1, xaxp=c(1, 4, 3), las=2)
	mtext(paste("Conc",i),side=1,line=0.5)
}

mtext("Mean Inter-Slide Spot Intensity CV", side = 3, line = -1.1, outer = TRUE)
}

#--------------------------------------------------------------------------
#Plot of inter-slide intensity
x11()
par(mfrow=c(8,2))
par(mar=c(1,2,1,1))

for(i in 1:16){
	if(test ==1){
	boxplot(x = as.list(as.numeric(aMat[i,,])),ylim=c(0,2e6))
 	mtext(paste("Well",i),side=1,line=0.25)
	}
	else{
	boxplot(x = as.list(as.numeric(aMean[i,,])),ylim=c(0,2e6))
      mtext(paste("Well",i),side=1,line=0.25)
	}
	axis(side = 1,xlab="",xaxt="n", axes =FALSE)
}
mtext("Inter-Slide Intensity", side = 3, line = -1.1, outer = TRUE)

#--------------------------------------------------------------------------
#Final mean inter-slide plot 
if(test ==1){
aMeaninter <- mat.or.vec(1,20)
aSdinter <- mat.or.vec(1,20)

for(i in 1:20){
	aMeaninter[i] <- mean(as.numeric(aMean[,i,]))
	aSdinter[i] <- sd(as.numeric(aMean[,i,]))
}
aCvinter <- aSdinter / aMeaninter *100
}else{
aMeaninter <- mat.or.vec(1,16)
aSdinter <- mat.or.vec(1,16)

for(i in 1:16){
	aMeaninter[i] <- mean(as.numeric(aMean[i,,]))
	aSdinter[i] <- sd(as.numeric(aMean[i,,]))
}
aCvinter <- aSdinter / aMeaninter * 100
}

#Plot Inter-slide Intensity and CV across Wells of same slide on two axis
x11() 
par(mfrow=c(1,1))
par(mar=c(4,4,2,4))

plot(as.numeric(aMeaninter),col="red",ylab="",xlab="",axes=FALSE, yaxt="n")
axis(2,ylim=c(0,4e6),col="red",col.axis="red",las=1, at = c(0,1e6,2e6,3e6,4e6))
mtext(side = 2, line = 3, 'Intensity')

par(new=TRUE)

plot(as.numeric(aCvinter),col="blue",ylab="",xlab="Spots",axes=FALSE, yaxt="n")
axis(1, xaxp=c(1, 20, 19), las=2)
axis(4,ylim=c(0,100),col="blue",col.axis="blue",las=1,at = c(0,20,40,60,80,100))

mtext(side = 4, line = 3, 'CV (%)')
title3 <- "Inter-Slide Mean Raw Intensity vs CV"
mtext(title3, side = 3, line = -1.1, outer = TRUE)
mtext(title3, side = 3, line = -1.1, outer = TRUE)
box()

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

