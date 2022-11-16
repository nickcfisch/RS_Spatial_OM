
#######################################################################################
#Getting continuous (effectively) depth for the GOM (Florida) and making a grid
#######################################################################################

#GOM Depth
#   Depth<-read.delim("C:/Users/nfisch/Documents/GOM Data/GOM_Depth.txt", header=FALSE, col.names = c("Longitude", "Latitude", "Depth"))

#Subset for Florida
#   Depth<-Depth[Depth$Longitude > -87.5 & Depth$Longitude < -81 & Depth$Latitude > 24.5 & Depth$Latitude < 30.5,]
#   Depth<-Depth[-c(which(Depth$Longitude> -82 & Depth$Latitude > 28)),] #This takes out St Johns and Atlantic Depths

#Subset for negative depths 
#   Depth<-Depth[Depth$Depth<0,]
#   Depth$Depth<--Depth$Depth  #Making Depths positive

#write.table(Depth,"C:/Users/nfisch/Documents/GOM_Data/Depth_FLA.txt", row.names=FALSE)

Depth<-read.delim("C:/Users/nfisch/Documents/GOM_Data/Depth_FLA.txt", sep=" ",header=TRUE)

#Looking at it in a map
par(mar=c(1,4,0,1))
plotMap(NWApolys_i, xlim=c(271.5,280), ylim=c(24, 32), bg=clr$sea, col=clr$land, tck=-0.014, las=1, plt=NULL)
D_pts<-as.EventData(data.frame("X"=Depth$Longitude+360, "Y"=Depth$Latitude, "EID"=seq(1,dim(Depth)[1]), "Z"=Depth$Depth), projection="LL")
addBubbles(D_pts, symbol.bg=colfunc_abun(ceiling(max(Depth$Depth))), min.size=0.0001, max.size=0.0001,legend.breaks=seq(0,ceiling(max(Depth$Depth)),1), symbol.zero=1, legend.pos=NULL)
par(mar=c(1,4,0,1))
addPolys(NWApolys_i,col=clr$land)

abline(v=seq(272.5,279,0.1), h=seq(24.5, 30.5, 0.1))

x_coords<-seq(272.5,279,0.1)
y_coords<-seq(24.5,30.5,0.1)

for (i in 1:length(x_coords)){
  lines(rep(x_coords[i],length(y_coords)),y_coords)
}
for (j in 1:length(y_coords)){
  lines(x_coords,rep(y_coords[j],length(x_coords)))
}

par(mar=c(1,4,0,1))
addPolys(NWApolys_i,col=clr$land)

#Creating data frame to hold data for grids 
Cells<-data.frame("Longitude"=rep(round(seq(272.55,278.95,0.1),2),each=length(seq(24.55,30.45,0.1))), "Latitude"=rep(round(seq(24.55,30.45,0.1),2),length(seq(272.55,278.95,0.1))))
Cells$Longneg<-Cells$Longitude-360
#for loop that gets average depth within the grids
for (i in 1:dim(Cells)[1]){
Cells[i,"Depth"]<-mean(Depth[Depth$Longitude > Cells$Longneg[i]-0.05 & Depth$Longitude < Cells$Longneg[i]+0.05 & Depth$Latitude < Cells$Latitude[i] + 0.05 & Depth$Latitude > Cells$Latitude[i] - 0.05,"Depth"])
}

#Subsetting for just GOM
Cells<-Cells[-which(is.na(Cells[,4])),] #Taking out cells that have NA depth
Cells<-Cells[-which(Cells$Depth < 10),] #Keeping cells With depths greater than 10 meters
Cells<-Cells[-which(Cells$Longitude >= 278.45 & Cells$Latitude <= 24.65),] #Taking out Atlantic cells
Cells<-Cells[-which(Cells$Depth>500),]
#write.table(Cells,"C:/Users/nfisch/Documents/GOM_Data/Depth_GOMFLA.txt", row.names=FALSE)

par(mar=c(1,4,0,1))
plotMap(NWApolys_i, xlim=c(271.5,280), ylim=c(24, 32), bg=clr$sea, col=clr$land, tck=-0.014, las=1, plt=NULL, main="Depth")
Cells_pts<-as.EventData(data.frame("X"=Cells$Longitude, "Y"=Cells$Latitude, "EID"=seq(1,dim(Cells)[1]), "Z"=Cells$Depth), projection="LL")
addBubbles(Cells_pts, symbol.bg=colfunc_abun(ceiling(max(Cells$Depth))), min.size=0.09, max.size=0.09,legend.breaks=seq(0,ceiling(max(Cells$Depth)),1), symbol.zero=1, legend.pos=NULL)

#Logged for visualization
par(mar=c(1,4,0,1))
plotMap(NWApolys_i, xlim=c(271.5,280), ylim=c(24, 32), bg=clr$sea, col=clr$land, tck=-0.014, las=1, plt=NULL, main="Depth")
Cells_log_pts<-as.EventData(data.frame("X"=Cells$Longitude, "Y"=Cells$Latitude, "EID"=seq(1,dim(Cells)[1]), "Z"=log(Cells$Depth)), projection="LL")
addBubbles(Cells_log_pts, symbol.bg=colfunc_abun(ceiling(max(log(Cells$Depth)))), min.size=0.09, max.size=0.09,legend.breaks=seq(0,ceiling(max(log(Cells$Depth))),1), symbol.zero=1, legend.pos=NULL)
addBubbles(Cells_log_pts, symbol.bg=heat.colors(ceiling(max(log(Cells$Depth))), alpha=0.5), min.size=0.09, max.size=0.09,legend.breaks=seq(0,ceiling(max(log(Cells$Depth))),1), symbol.zero=1, legend.type="horiz", legend.title="Depth (Logged)")
###########################################################


#Checking Zach Depth vs my depth 
x<-matrix(NA, nrow=dim(Cells)[1], ncol=2)
for (i in 1:dim(Cells)[1]){
  if(length(ssdf[which(ssdf$Lon > Cells$Longneg[i]-0.05 & ssdf$Lon < Cells$Longneg[i]+0.05 & ssdf$Lat < Cells$Latitude[i] + 0.05 & ssdf$Lat > Cells$Latitude[i] - 0.05),"Depth"])!=0){
  x[i,1:length(ssdf[which(ssdf$Lon > Cells$Longneg[i]-0.05 & ssdf$Lon < Cells$Longneg[i]+0.05 & ssdf$Lat < Cells$Latitude[i] + 0.05 & ssdf$Lat > Cells$Latitude[i] - 0.05),"Depth"])]<-ssdf[which(ssdf$Lon > Cells$Longneg[i]-0.05 & ssdf$Lon < Cells$Longneg[i]+0.05 & ssdf$Lat < Cells$Latitude[i] + 0.05 & ssdf$Lat > Cells$Latitude[i] - 0.05),"Depth"]
  }
}

#################################################################
##########################################
#Reading in hardbottom data 
##########################################
#################################################################
setwd("C:/Users/nfisch/Documents/Mapping")
library(PBSmapping)
library(maptools)
#Downloading high resolution world coastlines
data(worldLLhigh)
clr <- .PBSclr()
gshhg = "./gshhg-bin-2.3.0/"
#Map parameters
limits <- list(x = c(265, 290), y = c(20, 35))
NWApolys_i <- importGSHHS("Gshhg/gshhs_i.b", xlim = limits$x, limits$y, maxLevel = 4)  # automatically converts to °W

#Reading in shapefile and subsetting it to my area
hardbottom<-importShapefile(fn="C:/Users/nfisch/Documents/GOM_Data/GOM_Sediments/usSEABED_GOM_Sediments",readDBF=TRUE, projection = "LL")
hardbottom$X<-hardbottom$X+360
hardbottom<-hardbottom[-which(hardbottom$X < 272.5),] 
hardbottom<-hardbottom[-which(hardbottom$Y < 24.5),] 
hardbottom<-hardbottom[-which(hardbottom$X >= 278.45),] #Taking out Atlantic cells
hardbottom<-hardbottom[-which(hardbottom$X <= 275 & hardbottom$Y <= 27.5),] #Cutting out data that likely won't be included in model

#Full map
par(mar=c(1,4,0,1))
plotMap(NWApolys_i, xlim=c(271.5,280), ylim=c(24, 32), bg=clr$sea, col=clr$land, tck=-0.014, las=1, plt=NULL, main="Substrate Types")
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==3000)],], col = 2, border=2)
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==2000)],], col = 3, border=3)
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==300)],], col = 4, border=4)
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==200)],], col = 5, border=5)
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==30)],], col = 6, border=6)
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==20)],], col = 7, border=7)
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==3)],], col = 8, border=8)
addPolys(hardbottom[hardbottom$PID %in% attributes(hardbottom)$PolyData$PID[which(attributes(hardbottom)$PolyData$gom_domnc==2)],], col = 9, border=9)
legend("topright", c("Dominant Rock","Subdominant Rock","Dominant Gravel","Subdominant Gravel","Dominant Sand","Subdominant Sand","Dominant Mud","Subdominant Mud"), col=c(2,3,4,5,6,7,8,9), pch=16)

load("C:/Users/nfisch/Documents/RSnapper_Data/ssdf_drop_2.RDATA", verbose=T) 
head(ssdf)
ssdf_FLAsubset<-ssdf[ssdf$Lon > -87.5, ]
ssdf_FLAsubset$Lon<-ssdf_FLAsubset$Lon+360
points<-as.EventData(x=data.frame(EID=c(1:nrow(ssdf_FLAsubset)),X=ssdf_FLAsubset$Lon, Y=ssdf_FLAsubset$Lat), projection="LL")
addPoints(points)

#Adding a substrate column by extracting from attributes of the polydata
for (i in 1:nrow(hardbottom)){
 hardbottom$substrate[i]<-attributes(hardbottom)$PolyData[which(attributes(hardbottom)$PolyData$PID == hardbottom$PID[i]),12]
}
hardbottom$substrate<-ifelse(hardbottom$substrate==-99,0,hardbottom$substrate) #If a value is -99, then switch to zero code

#Equivalent Way of doing it 
par(mfrow=c(1,1),mar=c(1,4,0,1))
plotMap(NWApolys_i, xlim=c(271.5,280), ylim=c(24, 32), bg=clr$sea, col=clr$land, tck=-0.014, las=1, plt=NULL, main="Substrate Types")
addPolys(hardbottom[hardbottom$substrate==3000,], col = 2, border=2)
addPolys(hardbottom[hardbottom$substrate==2000,], col = 3, border=3)
addPolys(hardbottom[hardbottom$substrate==300,], col = 4, border=4)
addPolys(hardbottom[hardbottom$substrate==200,], col = 5, border=5)
addPolys(hardbottom[hardbottom$substrate==30,], col = 6, border=6)
addPolys(hardbottom[hardbottom$substrate==20,], col = 7, border=7)
addPolys(hardbottom[hardbottom$substrate==3,], col = 8, border=8)
addPolys(hardbottom[hardbottom$substrate==2,], col = 9, border=9)
legend("topright", c("Dominant Rock","Subdominant Rock","Dominant Gravel","Subdominant Gravel","Dominant Sand","Subdominant Sand","Dominant Mud","Subdominant Mud"), col=c(2,3,4,5,6,7,8,9), pch=16)

#Getting the points that are in polygons
Poly_points<-findPolys(points, hardbottom)
Poly_points$substrate<-rep(0,nrow(Poly_points))
for(i in 1:nrow(Poly_points)){  #Adding a substrate column for Polypoints
  Poly_points$substrate[i]<-hardbottom[which(hardbottom$PID==Poly_points$PID[i])[1],6]    
}

#Now adding substrate column to points data
points$substrate<-rep(0,nrow(points))
for(i in 1:nrow(points)){
  if(sum(Poly_points$EID==points$EID[i])!=0){
  points$substrate[i]<-Poly_points[which(Poly_points$EID==points$EID[i])[1],5]    
  }
}

#Making color vector
points$col_vec<-ifelse(points$substrate==3000,2,0)
points$col_vec<-ifelse(points$substrate==2000,3,points$col_vec)
points$col_vec<-ifelse(points$substrate==300,4,points$col_vec)
points$col_vec<-ifelse(points$substrate==200,5,points$col_vec)
points$col_vec<-ifelse(points$substrate==30,6,points$col_vec)
points$col_vec<-ifelse(points$substrate==20,7,points$col_vec)
points$col_vec<-ifelse(points$substrate==3,8,points$col_vec)
points$col_vec<-ifelse(points$substrate==2,9,points$col_vec)
#Adding substrate column to points dataframe
points_wsub<-as.EventData(x=data.frame(EID=points$EID,X=points$X, Y=points$Y, Z=points$substrate), projection="LL")
addPoints(points_wsub,col=points$col_vec, cex=0.5, pch=16)

#Adding substrate class to red snapper grid & catch database
ssdf_FLAsubset$substrate<-points_wsub$Z
#Making a column that has the actual names of the substrates
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==3000,"Dominant Rock","Unclassified")
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==2000,"Subdominant Rock",ssdf_FLAsubset$substrate_lab)
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==300,"Dominant Gravel",ssdf_FLAsubset$substrate_lab)
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==200,"Subdominant Gravel",ssdf_FLAsubset$substrate_lab)
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==30,"Dominant Sand",ssdf_FLAsubset$substrate_lab)
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==20,"Subdominant Sand",ssdf_FLAsubset$substrate_lab)
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==3,"Dominant Mud",ssdf_FLAsubset$substrate_lab)
ssdf_FLAsubset$substrate_lab<-ifelse(ssdf_FLAsubset$substrate==2,"Subdominant Mud",ssdf_FLAsubset$substrate_lab)

Agg_list_substrate<-list()
#par(mfrow=c(4,3), mar=c(5,2.5,4,2))
for (i in 1:11){
#  Agg_list_substrate[[i]]<-do.call(data.frame, aggregate(ssdf_FLAsubset[,50+i],by=list(ssdf_FLAsubset$substrate_lab), FUN=sum))
  Agg_list_substrate[[i]]<-do.call(data.frame, aggregate(ssdf_FLAsubset[,50+i],by=list(ssdf_FLAsubset$substrate), FUN=sum))
  Agg_list_substrate[[i]]$percent<-c(NA,Agg_list_substrate[[i]]$x[2:9]/sum(Agg_list_substrate[[i]]$x[2:9]))
  Agg_list_substrate[[i]]$label<-as.factor(c("Unclassified","Sub Mud", "Dom Mud", "Sub Sand", "Dom Sand", "Sub Gravel", "Dom Gravel", "Sub Rock" ,"Dom Rock"))
#  plot(Agg_list_substrate[[i]]$label, Agg_list_substrate[[i]]$percent, las=2, main=paste0("Age ",i-1), ylim=c(0,0.5))
}

par(mfrow=c(1,1), mar=c(6,3,1,1))
plot(Agg_list_substrate[[1]]$label, rep(NA,9), las=2, ylim=c(0,0.5), type="b", pch=16, xlab="", xaxt="n", xlim=c(0.5,9.5))
axis(side=1, at=1:9,labels = c("Dom Gravel", "Dom Mud","Dom Rock","Dom Sand","Sub Gravel","Sub Mud","Sub Rock","Sub Sand", "Unclassified"), las=2)
for(i in 1:11){
  points(Agg_list_substrate[[i]]$label, Agg_list_substrate[[i]]$percent, pch=16, col=i)
}
legend("topright", c("Age 0", "Age 1","Age 2", "Age 3","Age 4", "Age 5","Age 6", "Age 7","Age 8", "Age 9","Age 10", "Age 11"), col=1:11, pch=16)

#Writing Preference dataframe

Pref_substrate<-data.frame(Age0=Agg_list_substrate[[1]]$percent[2:9],Age1=Agg_list_substrate[[2]]$percent[2:9],
                 Age2=Agg_list_substrate[[3]]$percent[2:9],Age3=Agg_list_substrate[[4]]$percent[2:9],
                 Age4=Agg_list_substrate[[5]]$percent[2:9],Age5=Agg_list_substrate[[6]]$percent[2:9],
                 Age6=Agg_list_substrate[[7]]$percent[2:9],Age7=Agg_list_substrate[[8]]$percent[2:9],
                 Age8=Agg_list_substrate[[9]]$percent[2:9],Age9=Agg_list_substrate[[10]]$percent[2:9],
                 Age10=Agg_list_substrate[[11]]$percent[2:9])
#write.table(Pref_substrate,"C:/Users/nfisch/Documents/Snapper Simulation Assessments/substrate_pref.txt", row.names=FALSE)
Pref_sub<-read.delim("C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/substrate_pref.txt", sep=" ",header=TRUE)

############################################################################
#Ok now getting the substrate type in the cell matrix used in spatial model 
############################################################################
Cells<-read.delim("C:/Users/nfisch/Documents/GOM_Data/Depth_GOMFLA.txt", sep=" ",header=TRUE)

Cells_points<-as.EventData(x=data.frame(EID=c(1:nrow(Cells)),X=Cells$Longitude, Y=Cells$Latitude), projection="LL")

#Finding the substrate type of each cell points
Cells_Poly_points<-findPolys(Cells_points, hardbottom)
Cells_Poly_points$substrate<-rep(0,nrow(Cells_Poly_points))
for(i in 1:nrow(Cells_Poly_points)){
  Cells_Poly_points$substrate[i]<-hardbottom[which(hardbottom$PID==Cells_Poly_points$PID[i])[1],6]    
}

#Finding the substrate type of each cell points
Cells_points$substrate<-rep(0,nrow(Cells_points))
for(i in 1:nrow(Cells_points)){
  if(sum(Cells_Poly_points$EID==Cells_points$EID[i])!=0){
    Cells_points$substrate[i]<-Cells_Poly_points[which(Cells_Poly_points$EID==Cells_points$EID[i])[1],5]    
  }
}
  
Cells$substrate<-Cells_points$substrate

col_vec<-ifelse(Cells$substrate==3000,2,0)
col_vec<-ifelse(Cells$substrate==2000,3,col_vec)
col_vec<-ifelse(Cells$substrate==300,4,col_vec)
col_vec<-ifelse(Cells$substrate==200,5,col_vec)
col_vec<-ifelse(Cells$substrate==30,6,col_vec)
col_vec<-ifelse(Cells$substrate==20,7,col_vec)
col_vec<-ifelse(Cells$substrate==3,8,col_vec)
col_vec<-ifelse(Cells$substrate==2,9,col_vec)

addPoints(Cells_points,col=col_vec, cex=0.5, pch=16)

#Manually Filling in zeroes with closest cell, because I am lazy
Empty_substrate<-Cells[which(col_vec==0),]
Full_substrate<-Cells[which(col_vec!=0),]

Cells[46,5]<-30
Cells[173,5]<-200
Cells[193,5]<-200
Cells[233,5]<-200
Cells[273,5]<-20
Cells[313:315,5]<-20
Cells[331,5]<-3
Cells[347,5]<-3
Cells[383,5]<-3
Cells[408:409,5]<-c(3000,30)
Cells[435:436,5]<-c(3000,20)
Cells[752:755,5]<-30
Cells[801:804,5]<-30
Cells[853,5]<-2
Cells[854,5]<-20
Cells[1559,5]<-2000

#Coding for substrate
substrate_code<-ifelse(Cells$substrate==3000,8,0)
substrate_code<-ifelse(Cells$substrate==2000,7,substrate_code)
substrate_code<-ifelse(Cells$substrate==300,6,substrate_code)
substrate_code<-ifelse(Cells$substrate==200,5,substrate_code)
substrate_code<-ifelse(Cells$substrate==30,4,substrate_code)
substrate_code<-ifelse(Cells$substrate==20,3,substrate_code)
substrate_code<-ifelse(Cells$substrate==3,2,substrate_code)
substrate_code<-ifelse(Cells$substrate==2,1,substrate_code)

Cells$substrate_code<-substrate_code #Getting actual code from aggregate function (matches it with preference function)
#write.table(Cells,"C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/GOMFLA_Depth_n_substrate.txt", row.names=FALSE)

Cells<-read.delim("C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/GOMFLA_Depth_n_substrate.txt", sep=" ",header=TRUE)

par(mfrow=c(1,1),mar=c(4,4,1,1))
plotMap(NWApolys_i, xlim=c(271.5,280), ylim=c(24, 32), bg=clr$sea, col=clr$land, tck=-0.014, las=1, plt=NULL, main="Substrate Types")
addPoints(as.EventData(data.frame(EID=seq(1,nrow(Cells)),X=Cells$Longitude, Y=Cells$Latitude), projection="LL"),col=substrate_code, cex=1.2, pch=16)
legend(x=276,y=32, c("Dominant Rock","Subdominant Rock","Dominant Gravel","Subdominant Gravel","Dominant Sand","Subdominant Sand","Dominant Mud","Subdominant Mud"), col=c(8,7,6,5,4,3,2,1), pch=16)

#############
#CVs
#############
#Sample size vs ESS
par(mar=c(4,5,2,2))
plot(rowSums(Output_SG$Bio_sample_yr), xlab="Year", ylab="Composition Sample Size",las=1, ylim=c(0,1.2e5))

Comp_yr<-list()
for (i in 1:nyear){
  Comp_yr[[i]]<-Output_SG$Bio_samples[i,1,,]
  for (p in 2:num_ports){
    Comp_yr[[i]]<-rbind(Comp_yr[[i]],Output_SG$Bio_samples[i,p,,])
  }
  Comp_yr[[i]]<-Comp_yr[[i]][-which(is.na(Comp_yr[[i]][,1])),] #Taking out NAs
}

CV_SS<-matrix(NA, nrow=100, ncol=21)
for (i in 1:100){
  for (j in 1:21){
   CV_SS[i,j]<-sd(Comp_yr[[i]][,j])/mean(Comp_yr[[i]][,j])
  }
}

Expected_CV<-matrix(NA, nrow=100, ncol=21)
for (i in 1:100){
  for (j in 1:21){
    Expected_CV[i,j]<-((1-Output_SG$Age_Comp_Data[i,j])/(Output_SG$Age_Comp_Data[i,j]*rowSums(Output_SG$Bio_sample_yr)[i]))^0.5
  }
}

plot(0:20,CV_SS[6,], type="l", ylim=c(0,6), xlab="Age", ylab="CV of Composition", pch=16, col=3)
lines(0:20,Expected_CV[6,], lwd=2)
for (i in 7:100){
  lines(0:20,CV_SS[i,], col=3)
  lines(0:20,Expected_CV[i,], lwd=2)
}
legend("topleft", c("Expected_CV", "Actual CVs"), col=c(1,3), lwd=1, bty="n")

#Based on expected proportion
plot(Output_SG$Age_Comp_Data[6,],Expected_CV[6,], type="p", xlim=c(0,0.25),ylim=c(0,2), xlab="Expected Proportion", ylab="CV of Composition (across samples)", pch=16, col=1)
points(Output_SG$Age_Comp_Data[6,],CV_SS[6,], col=3, pch=16)
for (i in 7:100){
  points(Output_SG$Age_Comp_Data[i,],CV_SS[i,], col=3, pch=16)
  points(Output_SG$Age_Comp_Data[i,],Expected_CV[i,], lwd=2, pch=16)
}
legend("topleft", c("Expected CV", "Actual CVs"), col=c(1,3), lwd=1, bty="n")


##################################
#Red Snapper Tag Data
##################################
RSnap_Tag<-read.csv("C:/Users/nfisch/Documents/RSnapper_Data/Tag_snapper.csv")
RSnap_Tag<-RSnap_Tag[!is.na(RSnap_Tag[,1]),]

RSnap_Tag$Dist_Traveled<-earth.dist(RSnap_Tag[,"Lon1"], RSnap_Tag[,"Lat1"], RSnap_Tag[,"Lon2"], RSnap_Tag[,"Lat2"])
hist(RSnap_Tag$Dist_Traveled, xlab="Distance Traveled (km)", breaks=seq(0,400,5), main="Red Snapper Distance Traveled", las=1)

RSnap_Tag$Dist_Traveled_peryr<-RSnap_Tag$Dist_Traveled*(365/RSnap_Tag$days_at_large)
hist(RSnap_Tag$Dist_Traveled_peryr, xlab="Distance Traveled (km)", breaks=seq(0,1010,10), main="Red Snapper Distance Traveled Per Year", las=1)

func_move<-function(theta){
  lam_dist<-exp(theta[1])
  #p_fit<-exp(-lam_dist*RSnap_Tag$Dist_Traveled_peryr)
  #NLL<--1*sum(p_fit*(1-p_fit))
  
  NLL<--1*sum(dexp(RSnap_Tag$Dist_Traveled_peryr,lam_dist,log=T))
  #NLL<--1*sum(log(lam_dist*exp(-lam_dist*RSnap_Tag$Dist_Traveled*(RSnap_Tag$days_at_large/365))))
  return(NLL)
}
fit_move<-optim(log(0.05), fn=func_move, method="BFGS")
exp(fit_move$par)

dexp(RSnap_Tag$Dist_Traveled[1], rate=0.02)
0.02*exp(-0.02*RSnap_Tag$Dist_Traveled[1]*(RSnap_Tag$days_at_large[1]/RSnap_Tag$Dist_Traveled[1]))    #Density function in code

plot(0:500, exp(-exp(fit_move$par)*0:500), type="l", las=1, xlab="Distance (km)", ylab="Probability", lwd=2)
points(RS_Dist_Traveled_peryr_hist[RS_Dist_Traveled_peryr_hist$Freq!=0,"Dist"], RS_Dist_Traveled_peryr_hist[RS_Dist_Traveled_peryr_hist$Freq!=0,"Prop"])

hist(rexp(length(RSnap_Tag$Dist_Traveled_peryr),exp(fit_move$par)), breaks=seq(0,2000,1), xlim=c(0,1000), ylim=c(0,10), las=1)
hist(RSnap_Tag$Dist_Traveled_peryr, xlab="Distance Traveled (km)", breaks=seq(0,1010,1), xlim=c(0,1000), main="Red Snapper Distance Traveled Per Year", las=1, ylim=c(0,10), las=1)

RS_Dist_Traveled_peryr_hist<-as.data.frame(table(cut(RSnap_Tag$Dist_Traveled_peryr, breaks=seq(0,1010,1))))
RS_Dist_Traveled_peryr_hist$Dist<-seq(1,1010)

RS_Dist_Traveled_peryr_hist$Prop<-RS_Dist_Traveled_peryr_hist$Freq/max(RS_Dist_Traveled_peryr_hist$Freq)
plot(RS_Dist_Traveled_peryr_hist$Dist, RS_Dist_Traveled_peryr_hist$Prop, ylim=c(0,1), xlab="Distance Traveled (km)", pch=16, las=1, ylab="Frequency")

plot(RS_Dist_Traveled_peryr_hist[RS_Dist_Traveled_peryr_hist$Freq!=0,"Dist"], RS_Dist_Traveled_peryr_hist[RS_Dist_Traveled_peryr_hist$Freq!=0,"Prop"], ylim=c(0,1), xlab="Distance Traveled (km)", pch=16, las=1, ylab="Frequency")
move_distcost<-nls(Prop ~ exp(-lam_dist*Dist),data=RS_Dist_Traveled_peryr_hist[RS_Dist_Traveled_peryr_hist$Freq!=0,], start=list(lam_dist=0.03))
points(RS_Dist_Traveled_peryr_hist[RS_Dist_Traveled_peryr_hist$Freq!=0,"Dist"], predict(move_distcost), col=2)
lines(0:500, exp(-coef(move_distcost)*0:500), col=2, lwd=2)
lam_dist<-coef(move_distcost)

#Trying to correct for outlier movements
#Movement rate per day compared to total distance moved
plot(RSnap_Tag$Dist_Traveled/RSnap_Tag$days_at_large~RSnap_Tag$days_at_large, xlab="Days at Large", ylab="Daily Movement Rate", pch=16)

#Daily movement rate as a function of days at large, omitting fish that have been at large less than 200 days
plot(RSnap_Tag$Dist_Traveled[RSnap_Tag$days_at_large>200]/RSnap_Tag$days_at_large[RSnap_Tag$days_at_large>200]~RSnap_Tag$days_at_large[RSnap_Tag$days_at_large>200], xlab="Days at Large", ylab="Daily Movement Rate", pch=16)

#linear model 
lm_fit<-lm(log(RSnap_Tag$Dist_Traveled/RSnap_Tag$days_at_large)~RSnap_Tag$days_at_large)
lines(0:800, exp(-0.976303-0.004595*0:800)*exp(summary(lm_fit)$sigma^2/2), lwd=2, col=2)
lines(0:800, exp(-0.976303-0.004595*0:800+summary(lm_fit)$sigma^2/2), lwd=2, col=3)

#Movement per year
hist(RSnap_Tag$Dist_Traveled_peryr[RSnap_Tag$days_at_large>200], xlab="Distance Traveled (km)", breaks=seq(0,1010,10), main="Red Snapper Distance Traveled Per Year", las=1)

func_move<-function(theta){
  lam_dist<-exp(theta[1])

  NLL<--1*sum(dexp(RSnap_Tag$Dist_Traveled_peryr[RSnap_Tag$days_at_large>200],lam_dist,log=T))
  return(NLL)
}
fit_move<-optim(log(0.05), fn=func_move, method="BFGS")
exp(fit_move$par)

hist(rexp(length(RSnap_Tag$Dist_Traveled_peryr[RSnap_Tag$days_at_large>200]),exp(fit_move$par)), breaks=seq(0,2000,1), xlim=c(0,500), ylim=c(0,10), las=1)
hist(RSnap_Tag$Dist_Traveled_peryr[RSnap_Tag$days_at_large>200], xlab="Distance Traveled (km)", breaks=seq(0,1010,1), xlim=c(0,500), main="Red Snapper Distance Traveled Per Year", las=1, ylim=c(0,10), las=1)

#Proof of SSB0 from Red Snapper Assessment
M_vec<-c(2,1.2,0.19,0.15,0.129,0.115,0.106,0.099,0.095,0.091,0.088,0.086,0.085,0.083,0.082,0.081,0.081,0.08,0.08,0.079,0.078)
Adj_Mvec<-c(1,1.6,0.695,0.17,0.14,0.122,0.11,0.103,0.097,0.093,0.09,0.087,0.085,0.084,0.083,0.082,0.081,0.08,0.08,0.079,0.079)
Fec<-data.frame("Age"=0:20, "Fecundity"=c(0,0, 0.35E6, 2.62E6, 9.07E6, 20.3E6, 34.71E6, 49.95E6, 64.27E6, 76.76E6, 87.15E6, 95.53E6, 102.15E6, 107.3E6, 111.27E6, 114.3E6, 116.61E6, 118.36E6, 119.68E6, 120.67E6, 123.234591E6))

R0<-163212000
Num_age<-rep(NA, 21)
Num_age[1]<-R0
M<-Adj_Mvec
for(i in 2:21){
  Num_age[i]<-Num_age[i-1]*exp(-M[i-1]) #Unfished spawning biomass for Florida
}
Num_age[21]<-Num_age[21]/(1-exp(-M[21]))         #plus group
SSB0_FLA<-sum(Num_age*exp(-0.5*M)*Fec$Fecundity)  #How they say they calculate SSB
sum(Num_age*Fec$Fecundity)  #Actual value that lines up with SSB0

#########################################################
#Calculating the proportion of recruits going to Florida
#########################################################
setwd("C:/Users/nfisch/Documents/Mapping")
library(PBSmapping)
library(maptools)
#Downloading high resolution world coastlines
data(worldLLhigh)
clr <- .PBSclr()
gshhg = "./gshhg-bin-2.3.0/"
#Map parameters
limits <- list(x = c(265, 290), y = c(20, 35))
NWApolys_i <- importGSHHS("Gshhg/gshhs_i.b", xlim = limits$x, limits$y, maxLevel = 4)  # automatically converts to °W

plotMap(NWApolys_i, xlim=c(270,280), ylim=c(25, 31))
x<-ssdf[ssdf$Depth<70,]
pts<-as.EventData(data.frame("X"=x$Lon+360, "Y"=x$Lat, "EID"=seq(1,dim(x)[1]), "Z"=x$Depth), projection="LL")
addBubbles(pts, symbol.bg=2, min.size=0.0001, max.size=0.0001,legend.breaks=seq(0,ceiling(max(ssdf$Depth)),1), symbol.zero=1, legend.pos=NULL)

abline(v=-89+360)
abline(v=-87.5+360) #Florida Cutoff

Pred_Means[1]+2*sqrt(Pred_Vars[1])
Pref_sub[,1]

ssdf_FLAsub<-ssdf[ssdf$Lon > -87.5 & ssdf$Depth<70, ]
ssdf_EastMiss<-ssdf[ssdf$Lon > -89 & ssdf$Depth<70, ] 

dim(ssdf_FLAsub)[1]/dim(ssdf_EastMiss)[1]

#Plots for methods

#############################################################
#Making mega function faster by outsourcing some components
#############################################################
#Calculate euclidean distance in kilometers between two points
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.137   #Mean Earth Radius
  d <- R * c
  return(d)
}

dist_mat<-matrix(NA, nrow=num_cells, ncol=num_cells)
#Distance between cell i and cell j 
for (i in 1:num_cells){
  for (j in 1:num_cells){
    dist_mat[i,j]<-earth.dist(Cells[i,"Longneg"], Cells[i,"Latitude"], Cells[j,"Longneg"], Cells[j,"Latitude"])
  }
}

#write.table(dist_mat, file="C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/dist_mat.txt")

y<-as.matrix(read.table(file="C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/dist_mat.txt"))

####################################
#Getting County midpoints, etc.
####################################

#Crude county Gulf coast midpoints (county order south to north) There are 22 along the gulf coast
Gulf_County_Midpoints<-data.frame(County=c("Monroe","Collier","Lee","Charlotte","Sarasota","Manatee","Pinellas","Hillsborough","Pasco","Hernando","Citrus","Levy","Dixie","Taylor","Jefferson","Wakulla","Franklin","Gulf","Bay","Walton","Okaloosa","Santa Rosa","Escambia"),
                                  X=c(-80.94,-81.81,-82,-82.32,-82.53,-82.65,-82.85,-82.42,-82.73,-82.65,-82.64,-83.03,-83.31,-83.69,-84.02,-84.29,-84.68,-85.36,-85.72,-86.19,-86.59,-86.86,-87.19),
                                  Y=c(25.14,26.14,26.5,26.87,27.25,27.52,27.88,27.79,28.28,28.56,28.79,29.15,29.46,29.92,30.10,30.06,29.84,29.66,30.11,30.33,30.39,30.38,30.32),
                                  Pop=c(75027,378488,754610,184998,426718,394855,975280,1436888,539630,190865,147929,40770,16700,21623,14288,32461,11736,16164,185287,71375,207269,179349,315534))
#Adding column that has proportion of population in total (will use as proxy for total effort)
Gulf_County_Midpoints$Prop_pop<-Gulf_County_Midpoints$Pop/sum(Gulf_County_Midpoints$Pop)
#Number of ports
num_ports<-dim(Gulf_County_Midpoints)[1]                     #number of ports that effort goes out from

#Distance from each county midpoint to each cell 
County_Distance<-matrix(0,nrow=num_ports, ncol=num_cells)
for (i in 1:num_ports){
  for (j in 1:num_cells){
    County_Distance[i,j]<-earth.dist(Gulf_County_Midpoints$X[i], Gulf_County_Midpoints$Y[i], Cells$Longneg[j], Cells$Latitude[j])
  }
}

#write.table(County_Distance, file="C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/County_Distance.txt")

y<-as.matrix(read.table(file="C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/County_Distance.txt"))

sum(round(County_Distance)-round(y))

