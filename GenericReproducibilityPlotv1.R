#Run this block first - prompts for number of slides and filenames
rm(list=ls(all=TRUE))

#Set working directory
FILE <- file.choose()
DIR  <- dirname(FILE)
setwd(DIR)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Run the code in chunks  

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
		global[,m] <- a
	}
}
#Calculate Mean Intensity /STD per spot across one slide
#Split global into 16 x 25 matrices into wells

well <- array(0, dim = c(25,16,n))
count <- 25

for(m in 1:n){	
	for(i in 1:16){
		well[,i,m] <- global[(1+(count*(i-1))):(25+(count*(i-1))),m]
	}
}

aMean <- mat.or.vec(25,n)
aSd <- mat.or.vec(25,n)
aCv <- mat.or.vec(25,n)

for(m in 1:n){
	for(j in 1:25){
		aMean[j,m] <- mean(as.numeric(well[j,,m]))
		aSd[j,m] <- sd(as.numeric(well[j,,m]))
		}
}

aCv <- aSd / aMean * 100

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#Plot particular spots 
prompt <- "enter the spot numbers of interest (spaced between numbers or 25 for all)"
str <- as.integer(strsplit(readline(prompt), " ")[[1]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

prompt <- "enter the well numbers of interest (spaced between numbers or 16 for all)"
strWell <- as.integer(strsplit(readline(prompt), " ")[[1]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#If 25 , plot all 
if (str == 25){
	str <- 1:25
}

if (strWell == 16){
	strWell <- 1:16
}

transpose <- aperm(well,c(2,1,3))
wellSub <- transpose[strWell,str,]

#Convert to n data frames

df.names <- paste("dat",1:n,sep="")
yMax <- 0
for(i in 1:n){
	d.frame<-data.frame(wellSub[,,i])
	d.frame<-as.data.frame(lapply(d.frame,function(x){as.numeric(as.character(x))}))
	colnames(d.frame) <- c(paste('Spot',str))
	assign(df.names[i],d.frame)
	tempMax <- max(d.frame)
	if(tempMax > yMax){
		yMax = tempMax
	}
}

#Plot box whisker of 16 wells intra slide 
for(m in 1:n){
	x11()
	par(mfrow=c(1,1))
	par(mar=c(4,4,2,4))
	a <- eval(as.name(df.names[m]))

	boxdat <- boxplot(a,pch=19,ylab="",xlab="",axes=FALSE,ylim = c(0,yMax))	
      if(length(boxdat$out)!=0){
  		out.rows<-sapply(1:length(boxdat$out),function(i) which(a[,boxdat$group[i]]==boxdat$out[i]))
 		text(boxdat$group,boxdat$out,
   	      rownames(a)[out.rows],
   	      pos=3,
		offset = 0.75
  	)
}

	axis(2,ylim = c(0,max(eval(as.name(df.names[])))))
	mtext(side = 2, line = 3, 'Intensity')
	axis(1, at = 1:length(str), labels = c(str), las=2)
	mtext(side = 1, line = 3, 'Spots')
	if(length(strWell) == 16){
	title3 <- paste(tempF[m],"Intra-Slide Variation","Wells 1-16")
	}else{
	title3 <- paste(tempF[m],"Intra-Slide Variation","Wells",paste(strWell,collapse=" "))
	}
	mtext(title3, side = 3, line = -1.1, outer =TRUE)
	box()
	file <- paste(tempF[m],"Box",".pdf",sep="")
	dev.copy(pdf,file=paste(DIR,file,sep="/"))
	dev.off()
}


for(m in 1:n){
	x11()
	par(mfrow=c(1,1))
	par(mar=c(4,4,2,4))

	#plot(as.numeric(aMean[str,m]),pch = 19,col="red",ylab="",xlab="",axes=FALSE, yaxt="n")
	plot(as.numeric(aMean[str,m]),pch = 19,col="red",ylab="",xlab="",axes=FALSE,ylim = c(0,max(aMean[str,])))
	lines(1:length(str),aMean[str,m],col="red")
	#axis(2,ylim=c(0,4e6),col="red",col.axis="red",las=1, at = c(0,1e6,2e6,3e6,4e6))
	axis(2,col="red",col.axis="red",ylim = c(min(aMean[str,]),max(aMean[str,])))
	mtext(side = 2, line = 3, 'Intensity',col = "red")

	par(new=TRUE)
	#plot(as.numeric(aCv[str,m]),pch=19,col="blue",ylab="",xlab="Spots",axes=FALSE, yaxt="n")
	plot(as.numeric(aCv[str,m]),pch=19,col="blue",ylab="",xlab="Spots",axes=FALSE,ylim = c(0 ,max(aCv[str,])))
	lines(1:length(str),aCv[str,m],col="blue")
	axis(1, at = 1:length(str), labels = c(str), las=2)
	#axis(4,ylim=c(0,100),col="blue",col.axis="blue",las=1,at = c(0,20,40,60,80,100))
	axis(4,col="blue",col.axis="blue",ylim = c(min(aMean[str,]),max(aMean[str,])))
	abline(h = 20, lty = 2, col = "black")
	mtext(side = 4, line = 3, 'CV (%)', col = "blue")
	if(length(strWell) == 16){
	title2 <- paste(tempF[m],"Average Intra-Slide Intensity vs CV","Wells 1-16")
	}else{
	title2 <- paste(tempF[m],"Average Intra-Slide Intensity vs CV","Wells",paste(strWell,collapse=" "))
	}
	mtext(title2, side = 3, line = -1.1, outer = TRUE)
	box()
	file <- paste(tempF[m],"IntraSlide",".pdf",sep="")
	dev.copy(pdf,file=paste(DIR,file,sep="/"))
	dev.off()
}


#Plot inter-slide scatter differences 
x11()
par(mfrow=c(1,1))
par(mar=c(4,4,2,4))
color <- rainbow(n)

b<- 26:50
for(m in 1:n){
	#plot(as.numeric(aMean[str,m]),pch = 19, col = palette()[m], ylab="",xlab="",axes=FALSE,ylim = c (0,max(aMean[str,])))
	#lines(1:length(str),aMean[str,m],col=palette()[m])
	#axis(2,ylim = c(min(aMean[str,]),max(aMean[str,])))
	plot(as.numeric(global[b,m]),pch = 19, col = color[m], ylab="",xlab="",axes=FALSE,ylim = c (0,max(as.numeric(global[b,]))))
	lines(1:length(str),global[b,m],col=color[m])
	axis(2,ylim = c(min(global[b,]),max(global[b,])))
	mtext(side = 2, line = 3, 'Intensity')
	par(new=TRUE)
}
axis(1, at = 1:length(str), labels = c(str), las=2)
mtext(side = 1, line = 3, 'Spots')
legend("topright",c(tempF[1:m]),lty=c(1,1),lwd=c(2.5,2.5),col=color[1:m])
if(length(strWell) == 16){
	mtext(side = 3, line = -1.1, outer = TRUE, 'Inter-Slide Average Intensity Comparison Wells 1-16')
}else{
	title2 <- paste(tempF[m],"Inter-Slide Average Intensity Comparison","Wells",paste(strWell,collapse=" ")) 
	mtext(title2, side = 3, line = -1.1, outer = TRUE)
}
box()

file <- paste(tempF[m],"InterSlide",".pdf",sep="")
dev.copy(pdf,file=paste(DIR,file,sep="/"))
dev.off()

#Plot inter-slide CV differences
x11()
par(mfrow=c(1,1))
par(mar=c(4,4,2,4))

for(m in 1:n){
	plot(as.numeric(aCv[str,m]),pch = 19, col = palette()[m], ylab="",xlab="",axes=FALSE,ylim = c (0,max(aCv[str,])))
	lines(1:length(str),aCv[str,m],col=palette()[m])
	axis(2,ylim = c(min(aCv[str,]),max(aCv[str,])))
	mtext(side = 2, line = 3, 'CV (%)')
	par(new=TRUE)
}

axis(1, at = 1:length(str), labels = c(str), las=2)
mtext(side = 1, line = 3, 'Spots')
legend("topleft",c(tempF[1:m]),lty=c(1,1),lwd=c(2.5,2.5),col=palette()[1:m])
if(length(strWell) == 16){
	mtext(side = 3, line = -1.1, outer = TRUE, 'Inter-Slide Average CV Comparison Wells 1-16')
}else{
	title2 <- paste(tempF[m],"Inter-Slide Average CV Comparison","Wells",paste(strWell,collapse=" "))
	mtext(title2, side = 3, line = -1.1, outer = TRUE)
}
abline(h = 20, lty = 2, col = "black")
box()
file <- paste(tempF[m],"InterSlideCV",".pdf",sep="")
dev.copy(pdf,file=paste(DIR,file,sep="/"))
dev.off()




