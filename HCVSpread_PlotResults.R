######################################################################
# PROJECT: SPATIAL SIMULATION OF VIRAL SPREAD WITHIN A 2D LAYER OF CELLS
#
# R-file to plot the spatial distributions
# @Authors: Peter Kumberger & Frederik Graw
######################################################################

rm(list=ls())

# APPROACH: Use a matrix of size lines x columns in combination with levelplot function
library(lattice)
library(plotrix)
#setwd(path.figure)


#######
# (1.) Function to fill the layer
fill.layer <- function(vec,status,matrix.old,length_side){
	
	if (length(vec)>0){
		for (i in 1:length(vec)){
      # calculate the position in axial coordinates for the id
		  xpos <- floor((vec[i]-1)/(2*length_side-1))+1
		  ypos <- (vec[i]-1)%%(2*length_side-1)+1
      
      ## Convert from axial to offset coordinates
      ypos <- ypos + (xpos+xpos%%2)/2 - ceiling(length_side/2)

			matrix.old[xpos,ypos] <- status
		}
	}
	matrix.old
}
	

#######
# (2.) Determine the development for one particular simulation and save them into a list of matrices
# path.data <- path on your directory where the data are stored
# path.figure <- path where you want to store your figures
# name - name of the scenario that was simulated and used to save the data
# timespan - vector with the time points at which the plots should be generated
# simu - number of the simulation that should be plotted (if several repeats were run for one scenario)
# length_side - dimension of the hexagonal grid
# factor - factor to scale the timespan (minutes, days, seconds, ..)
# unit - additional factor to scale the time span (default set to 1)


plot.timeseries <- function(path.data,path.figure,name,simu,timespan,length_side,factor,unit=1){
    
  ## Calculate total number of cells an,d maximal ID
	cell.no <- 3*length_side^2-3*length_side+1
    max.id <- (2*length_side-1)^2-length_side+1
	out <- list(NULL)
# read empty cells and get data set
	Emp.data <- read.table(paste(path.data,"Empty_",name,".csv",sep=""),sep=",",fill=TRUE,
    header=FALSE,col.names=c("sim","time",1:cell.no))
    Empty <- subset(Emp.data,Emp.data$sim==simu)
# read infected cells and get data set
	Inf.data <- read.table(paste0(path.data,"HEPinf_",name,".csv"),sep=",",fill=TRUE,
						   header=FALSE,col.names=c("sim","time",1:cell.no))
	Infected <- subset(Inf.data,Inf.data$sim==simu)
# read infectious cells and get data set
	InfT.data <- read.table(paste0(path.data,"HEPinfectious_",name,".csv"),sep=",",fill=TRUE,
						   header=FALSE,col.names=c("sim","time",1:cell.no))
	Infectious <- subset(InfT.data,InfT.data$sim==simu)
# read the type of infection for each cell
	InfectType.data <- read.table(paste0(path.data,"HEPInfectClusterType_",name,".csv"),sep=",",header=FALSE,col.names=c("ID","infType","clustID"))
  	InfectType <- InfectType.data[cell.no*simu+(1:cell.no),]


# create the matrix with the different cells
  layer <- matrix(NA,nrow=2*length_side-1,ncol=2*length_side-1)
  layer <- fill.layer(vec = InfectType$ID,status = 0,matrix.old = layer,length_side)
	
# (I.) fill the start with the first line of timespan
	J <- subset(Infected,Infected$time==timespan[1])[,-c(1,2)]
	J <- J[!is.na(J)]
  InfectTypeCC <- subset(InfectType,InfectType[,2]==2 | (InfectType[,2]==4)) # Infection by free diffusion and cell to cell transmission is shown as cell to cell transmission
  JCC <- subset(J,is.element(J,InfectTypeCC[,1]))
  InfectTypeFD <- subset(InfectType,(InfectType[,2]==3)) # Infection by free diffusion
  JFD <- subset(J,is.element(J,InfectTypeFD[,1]))
  InfectTypeB <- subset(InfectType,InfectType[,2]==1) # Infectious cells at the beginning of the simulation
  JB <- subset(J,is.element(J,InfectTypeB[,1]))

	I <- subset(Infectious,Infectious$time==timespan[1])[,-c(1,2)]
	I <- I[!is.na(I)]
  ICC <- subset(I,is.element(I,InfectTypeCC[,1]))
  IFD <- subset(I,is.element(I,InfectTypeFD[,1]))
  IB <- subset(I,is.element(I,InfectTypeB[,1]))
  E <- subset(Empty,Empty$time==timespan[1])[,-c(1,2)]
  E <- E[!is.na(E)]

    layer <- fill.layer(JB,status=7,layer,length_side)
    layer <- fill.layer(IB,status=1,layer,length_side)
	layer <- fill.layer(JCC,status=2,layer,length_side)
    layer <- fill.layer(ICC,status=3,layer,length_side)
    layer <- fill.layer(JFD,status=4,layer,length_side)
    layer <- fill.layer(IFD,status=5,layer,length_side)
    layer <- fill.layer(E,status=6,layer,length_side)

	out[[1]] <- layer
	hour <- floor(timespan[1]/factor)
	minute <- timespan[1]%%factor*unit
	if (minute==0){
		minute <- "00"
	}
	title <- paste("Time ",hour,":",minute,sep="")
	
	folder <- paste(name,"/",sep="")
	path.figure <- paste0(path.figure,name,"/figures/")
  system(paste("mkdir ",path.figure,sep=""))
  path.figure <- paste0(path.figure,folder)
  system(paste("mkdir ",path.figure,sep=""))	

## Color code for plots later on
cellcol <- matrix(NA,nrow=nrow(layer),ncol=ncol(layer))
    cellcol[layer==0] <- "#a7a7a7"
    cellcol[layer==1] <- "#ff00ff"
    cellcol[layer==2] <- "#ffd891"
    cellcol[layer==3] <- "#ffa500"
    cellcol[layer==4] <- "#9191ff"
    cellcol[layer==5] <- "#00008b"
    cellcol[layer==6] <- "white"
    cellcol[layer==7] <- "#ff91ff"

  # plot the matrix and save the plot
par(oma=c(4,0,0,0))
pdf(paste(path.figure,"fig_",name,"_",str_pad(timespan[1],nchar(max(timespan)),pad=0),".pdf",sep=""),width=5,height=5,bg="transparent")
print(color2D.matplot(layer,cellcolors=cellcol,show.legend=F,xlim=c(0,2*length_side-1),ylim=c(0,2*length_side-1),xlab="",ylab="",do.hex=T,show.values=F,na.color="white",border=NA))
legend.col <- c("#a7a7a7","#ff00ff","#ffd891", "#ffa500","#9191ff","#00008b","white","#ff91ff")
legend.vec <- c("Uninfected","Seed_Infected","Seed_Infectious","Infected_CC","Infectious_CC","Infected_CF","Infectious_CF","Empty")
color.legend(-10,-25,110,-20,legend=legend.vec,align="rb",gradient="x",legend.col,cex=.6)
dev.off()

# (II.) change the "new" ones
for (i in 2:length(timespan)){
    print(timespan[i])
    # new infected
    Jn <- subset(Infected,Infected$time==timespan[i])[,-c(1,2)]
    Jn <- Jn[!is.na(Jn)]
    J.next <- setdiff(Jn,J)
    JBn <- subset(Jn,is.element(Jn,InfectTypeB[,1]))
    JB.next <- setdiff(JBn,JB)
    JCCn <- subset(Jn,is.element(Jn,InfectTypeCC[,1]))
    JCC.next <- setdiff(JCCn,JCC)
    JFDn <- subset(Jn,is.element(Jn,InfectTypeFD[,1]))
    JFD.next <- setdiff(JFDn,JFD)
    
    layer <- fill.layer(JB.next,status=7,layer,length_side)
    layer <- fill.layer(JCC.next,status=2,layer,length_side)
    layer <- fill.layer(JFD.next,status=4,layer,length_side)
    J <- Jn
    JB <- JBn
    JCC <- JCCn
    JFD <- JFDn
    
    # new infectious and targets
    In <- subset(Infectious,Infectious$time==timespan[i])[,-c(1,2)]
    In <- In[!is.na(In)]
    I.next <- setdiff(In,I)
    IBn <- subset(In,is.element(In,InfectTypeB[,1]))
  	IB.next <- setdiff(IBn,IB)
  	ICCn <- subset(In,is.element(In,InfectTypeCC[,1]))
  	ICC.next <- setdiff(ICCn,ICC)
  	IFDn <- subset(In,is.element(In,InfectTypeFD[,1]))
  	IFD.next <- setdiff(IFDn,IFD)
    En <- subset(Empty,Empty$time==timespan[i])[,-c(1,2)]
    E.next <- as.numeric(setdiff(E,En))
    
    layer <- fill.layer(E.next,status=0,layer,length_side)
  	layer <- fill.layer(IB.next,status=1,layer,length_side)
    layer <- fill.layer(ICC.next,status=3,layer,length_side)
  	layer <- fill.layer(IFD.next,status=5,layer,length_side)
    
    
    I <- In
    IB <- IBn
  	ICC <- ICCn
  	IFD <- IFDn
    E <- En

		
		out[[i]] <- layer
		
		hour <- floor(timespan[i]/factor)
		minute <- timespan[i]%%factor*unit
		if (minute==0){
			minute <- "00"
		}
		title <- paste("Time ",hour,":",minute,sep="")
    
		## Color code for plots later on
		cellcol <- matrix(NA,nrow=nrow(layer),ncol=ncol(layer))
        cellcol[layer==0] <- "#a7a7a7"
        cellcol[layer==1] <- "#ff00ff"
        cellcol[layer==2] <- "#ffd891"
        cellcol[layer==3] <- "#ffa500"
        cellcol[layer==4] <- "#9191ff"
        cellcol[layer==5] <- "#00008b"
        cellcol[layer==6] <- "white"
        cellcol[layer==7] <- "#ff91ff"
    
  # plot the new matrix
    par(oma=c(4,0,0,0))
  	pdf(paste(path.figure,"fig_",name,"_",str_pad(timespan[i],nchar(max(timespan)),pad=0),".pdf",sep=""),width=5,height=5,bg="transparent")
  	print(color2D.matplot(layer,cellcolors=cellcol,show.legend=F,xlim=c(0,2*length_side-1),ylim=c(0,2*length_side-1),xlab="",ylab="",do.hex=T,show.values=F,na.color="white",border=NA))
    #   legend.col <- c("#a7a7a7","#ff00ff","#ffd891", "#ffa500","#9191ff","#00008b","white","#ff91ff")
    #legend.vec <- c("Uninfected","Seed_Infected","Seed_Infectious","Infected_CC","Infectious_CC","Infected_CF","Infectious_CF","Empty")
    # color.legend(-10,-25,110,-20,legend=c("Uninfected","Seed","Infected_CC","Infectious_CC","Infected_FD","Infectious_FD"),align="rb",gradient="x",testcol,cex=.6)
  	dev.off()
	}
	out
}


######################################################################
# CALCULATE THE MATRICES AND PLOT THEM INTO A PDF-FILE


# save every 4 hours a picture
hour.step <- seq(60,5760,60)
#hour.step <- 14400
name <- "name of scenario"
LENGTH_SIDE <- 60
simulation.plot <- plot.timeseries(path.data,path.figure,name,0,hour.step,LENGTH_SIDE,factor=60,unit=1)




