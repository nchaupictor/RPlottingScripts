#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Generic Template for Membrane Intensity Variations
#Date: 27/04/16

#Ensure raw intensity results from results.xls has been copied and pasted into a new excel
#document and saved as a .csv file. The csv should have the following format:

#Experiment Number / Vendor / Membrane .csv

# e.g. 68E08_260616_8 M biodyne.csv


#Clear current plots
rm(list=ls(all=TRUE))

#Set working directory

FILE <- file.choose()
DIR  <- dirname(FILE)
setwd(DIR)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#File input 
n <- -1
while(n < 0){ 
 	 n <- readline("Enter number of slides: ")
 	 n <- as.integer(n)
 	 if (is.na(n)){
  		  n <- readline("Enter number of slides: ")
 		  n <- as.integer(n)
  	}
	tempF <- mat.or.vec(n,1)
	temp <- mat.or.vec(n,1)
 	 for(i in 1:n){
		tempF[i] <- readline("Enter file name (excluding .csv): ")
	      temp[i] <- paste(tempF[i] , ".csv",sep ="")
	}
}

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


strLabels <- strsplit(tempF, " ")
n <- as.integer(n)

#Parse in results.xlsx and sort  

global <- mat.or.vec(400,n) #Take all 25 spots in 16 wells 

for(m in 1:n){

	data <- read.csv(temp[m],header=TRUE)
	data1 <- as.matrix(data)

	count <- 0
	count2 <- 0
	count3 <- 1

	# Convert raw data into a single 1x320 vector sorting from Spot 1, Spot 6, Spot 2, Spot 7, etc
	a <- mat.or.vec(400, 1)

	for(i in 1:16){
		for(j in 1:5){
			for(k in 1:5){
					a[count3] <- data1[k+count,j+count2]
					count3 <- count3 + 1	
			}
		}
			count <- count+7
			if(i == 8){
				count <- 0
				count2 <- count2 + 6
			}
		if(n != 1){
		global[,m] <- a
		}else{
		global <- a
		}
	}
}
#Split global into 16 x 25 matrices into wells

well <- array(0, dim = c(25,16,n))
count <- 25

if(n != 1){
	for(m in 1:n){	
		for(i in 1:16){
			well[,i,m] <- global[(1+(count*(i-1))):(25+(count*(i-1))),m]
		}
	}
}else{
	for(i in 1:16){
		well[,i,1] <- global[(1+(count*(i-1))):(25+(count*(i-1)))]
	}
}

#Input desired columns
prompt <- "enter the columns of interest (column 1 = spots 6 7 8 9, column 2 = 11 12 13 14)"
columns <- as.integer(strsplit(readline(prompt), " ")[[1]])

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Input desired wells 
prompt <- "enter the well number of interest: "
strWell <- as.integer(strsplit(readline(prompt), " ")[[1]])

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Input number of membranes  
memNum <- -1
while(memNum < 0){
memNum <- readline("Enter number of membranes: ")
memNum <- as.integer(memNum)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
venNum <- -1
while(venNum < 0){
venNum <- readline("Enter number of vendors: ")
venNum <- as.integer(venNum)	
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Plot venNum by memNum plot showing a single sample intensity values for A1,A2,A3,A4
count <- 0
tempA <- mat.or.vec(4,4)
tempA[,1] <- c(6,7,8,9)
tempA[,2] <- c(11,12,13,14)
tempA[,3] <- c(16,17,18,19)
tempA[,4] <- c(21,22,23,24)

str <- (tempA[,columns])



for(j in 1 : (length(str)/4)){
	x11()
	par(mfrow=c(venNum,memNum))
	par(mar=c(2,4,3,1))
	par(oma=c(0,0,2,0))

	for(i in 1:(venNum*memNum)){
		maxBar <- max(well[str[,j],strWell,i])
		cols <- c("blue", "red")[(well[str[,j],strWell,i] == maxBar) +1] 
		barplot(as.numeric(well[str[,j],strWell,i]),col = cols,axes=FALSE,ylab="",xlab="",ylim = c(0,as.numeric(max(as.numeric(well[str[,j],strWell,])))*1.05))

		#barplot(as.numeric(well[spotD,strWell,i]),col = cols,axes=FALSE,ylab="",xlab="",ylim = c(0,as.numeric(max(as.numeric(well[spotD,strWell,])))*1.05))

		axis(2,ylim = c(0,max(as.numeric(well[str[,j],strWell,]))))
		#axis(2, ylim = c(0,max(as.numeric(well[spotD,strWell,]))))
		mtext(side = 1, line = 0.75,'Conc',cex = 0.75)
		mtext(side = 3, line = 0.75, strLabels[[i]][3])
		if(i == 1){
			mtext(side = 2, line = 2, strLabels[[1]][2],col="black")
		}
		if(i == (venNum*memNum/2)+1){
			mtext(side = 2, line = 2, strLabels[[venNum*memNum/2+1]][2],col ="steelblue1")
		}
		if(i <= (venNum + memNum)/2){
		box(col="black")	
		}else{
		box(col="steelblue1")
		}
	}
	#mtext(paste(strLabels[[1]][1],'NS-1','Sample',strWell),outer =TRUE, cex = 1.5)
	#file <- paste(strLabels[[1]][1],"NS-1","Sample",strWell,".pdf",sep="")
	mtext(paste(strLabels[[1]][1],'Antigen',j,'Sample', strWell), outer = TRUE, cex = 1.5)
	file <- paste(strLabels[[1]][1],"Antigen",j,"Sample",strWell,".pdf",sep="")
	dev.copy(pdf,file=paste(DIR,file,sep="/"))
	dev.off()
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
for(k in 1:2){
	#Plot box whisker of 16 wells intra slide 
	x11()
	par(mfrow=c(venNum,memNum))
	par(mar=c(4,4,3,1))
	par(oma=c(0,0,2,0))

	if(k == 1){
		str <- c(1,5,25)
	}else{
		str <- c(2,4,15,20)
	}
	strWell2 = strWell
	#strWell2 <- seq(1,16)

	transpose <- aperm(well,c(2,1,3))
	wellSub <- transpose[strWell2,str,]


	#Convert to n data frames

	df.names <- paste("dat",1:n,sep="")
	yMax <- 0
	for(i in 1:n){
		#d.frame<-data.frame(wellSub[,,i])
		d.frame<-data.frame(wellSub[,])
		d.frame<-as.data.frame(lapply(d.frame,function(x){as.numeric(as.character(x))}))
		colnames(d.frame) <- c(paste('Spot',str))
		assign(df.names[i],d.frame)
		tempMax <- max(d.frame)
		if(tempMax > yMax){
			yMax = tempMax
		}
	}

	for(m in 1:n){
	
		a <- eval(as.name(df.names[m]))

		boxdat <- boxplot(a,pch=19,ylab="",xlab="",axes=FALSE,ylim = c(0,yMax))	
      	if(length(boxdat$out)!=0){
  			out.rows<-sapply(1:length(boxdat$out),function(i) which(a[,boxdat$group[i]]==boxdat$out[i]))
 			text(boxdat$group,boxdat$out,
   	      	rownames(a)[out.rows],
   	      	pos=3,
			offset = 0.75)
		}

		axis(2,ylim = c(0,max(eval(as.name(df.names[])))))
		axis(1, at = 1:length(str), labels = c(str), las=2)
		mtext(side = 1, line = 3, 'Spots')
		mtext(side = 3, line = 0.75, strLabels[[m]][3])
		if(m == 1){
			mtext(side = 2, line = 2, strLabels[[1]][2])
		}
		if(m == (venNum*memNum/2)+1){
			mtext(side = 2, line = 2, strLabels[[venNum*memNum/2+1]][2],col="steelblue1")
		}
		if(m <= 4){
			box()
		}else{
			box(col = "steelblue1")
		}

	}
		if(k == 1){
			mtext(paste(strLabels[[1]][1],'Average Controls'), side = 3, line = 0.5, outer =TRUE)
			file <- paste(strLabels[[1]][1],"Average Controls",".pdf",sep="")
			dev.copy(pdf,file=paste(DIR,file,sep="/"))
			dev.off()
		}else{
			mtext(paste(strLabels[[1]][1],'Average Print Buffer'), side = 3, line = 0.5, outer =TRUE)
			file <- paste(strLabels[[1]][1],"Average Print Buffer",".pdf",sep="")
			dev.copy(pdf,file=paste(DIR,file,sep="/"))
			dev.off()
		}
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

	spotD <- 15
	x11()
	par(mfrow=c(venNum,memNum))
	par(mar=c(2,4,3,1))
	par(oma=c(0,0,2,0))

	for(i in 1:(venNum*memNum)){
		#maxBar <- max(well[str[,j],strWell,i])
		#cols <- c("blue", "red")[(well[str[,j],strWell,i] == maxBar) +1] 
		#barplot(as.numeric(well[str[,j],strWell,i]),col = cols,axes=FALSE,ylab="",xlab="",ylim = c(0,as.numeric(max(as.numeric(well[str[,j],strWell,])))*1.05))

		barplot(as.numeric(well[spotD,strWell,i]),col = cols,axes=FALSE,ylab="",xlab="",ylim = c(0,as.numeric(max(as.numeric(well[spotD,strWell,])))*1.05))

		#axis(2,ylim = c(0,max(as.numeric(well[str[,j],strWell,]))))
		axis(2, ylim = c(0,max(as.numeric(well[spotD,strWell,]))))
		mtext(side = 1, line = 0.75,'Conc',cex = 0.75)
		mtext(side = 3, line = 0.75, strLabels[[i]][3])
		if(i == 1){
			mtext(side = 2, line = 2, strLabels[[1]][2],col="black")
		}
		if(i == (venNum*memNum/2)+1){
			mtext(side = 2, line = 2, strLabels[[venNum*memNum/2+1]][2],col ="steelblue1")
		}
		if(i <= (venNum+memNum)/2){
		box(col="black")	
		}else{
		box(col="steelblue1")
		}
	}
	mtext(paste(strLabels[[1]][1],'NS-1','Sample',strWell),outer =TRUE, cex = 1.5)
	file <- paste(strLabels[[1]][1],"NS-1","Sample",strWell,".pdf",sep="")
	#mtext(paste(strLabels[[1]][1],'Antigen',j,'Sample', strWell), outer = TRUE, cex = 1.5)
	#file <- paste(strLabels[[1]][1],"Antigen",j,"Sample",strWell,".pdf",sep="")
	dev.copy(pdf,file=paste(DIR,file,sep="/"))
	dev.off()




