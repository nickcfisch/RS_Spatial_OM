
##############################
#Function for running model
##############################
rm(list=ls(all=TRUE))

Spatial_Model<-function(save_wd,               #Working directory for saving output and data
                        seed,                  #random number seed
                        PE,                    #Process Error (turned on or off via T/F)
                        nyear,                 #Number of years in simulation, first 50 years will be unfished 
                        nyear_init,            #Number of unfished years to start simulation    
                        num_ages,              #Number of ages
                        lam_data,              #Whether or not fish movement cost is based on Data
                        lam_moveCost,          #The distance cost function for fish movement
                        DD_rate,               #DD rate for negative exponential
                        DD_Thresh_quant,       #Quantile of unfished cells used for DD threshold
                        sig_logR,              #Recruitment variation
                        q,                     #Fishery catchability
                        sel_grow,              #Fishery Selectivity (logistic growth rate)
                        sel_midpt,             #Fishery Selectivity (logistic midpoint)
                        lam_Costdist,          #Negative exponential parameter for fishing distance cost function 
                        Abunpref_grate,        #Logistic growth rate of fisher abundance preference
                        eff_scalar,            #Scalar for effort relationship
                        cv_totaleff,           #CV of total effort (puts spread around logistic) 
                        eff_midpt,             #midpoint of effort logistic
                        eff_grate,             #Logistic growth rate for effort logistic
                        cv_harv,               #CV of sampled harvest
                        cv_cpue,               #CV of CPUE
                        cv_effort,             #CV of sampled effort
                        Perc_eff_smpld,        #Percentage of effort sampled
                        Prop_sample,           #Are we proportionally sampling fish from the fishery catches or not
                        Perc_sample_pb,        #Percentage of samples taken per boat 
                        Sample_size_pb,        #Number of fish sampled from each trip
                        Perc_yrs_FIM,          #Percent of years with FIM data
                        FIM_q,                 #FIM q
                        num_FIM_samples,       #right now saying they take 20 trips a year
                        w_dist,                #power weight for fish movement distance 
                        w_depth,               #power weight for depth preference
                        w_substrate,           #power weight for substrate preference
                        w_DD,                  #power weight for DD preference
                        w_cost,                #power weight for fisher cost function
                        w_profit,              #power weight for fisher profit function
                        Super_Comp)            #Is this being run on the super computer?
{
  set.seed(seed)
  
  #Natural Mortality from table 2.1 in Red Snapper Assessment
  M_vec<-c(2,1.2,0.19,0.15,0.129,0.115,0.106,0.099,0.095,0.091,0.088,0.086,0.085,0.083,0.082,0.081,0.081,0.08,0.08,0.079,0.078)
  M<-data.frame("Age"=c(0:20),"M"=M_vec)
  
  #Fecundity Table
  Fec<-data.frame("Age"=0:20, "Fecundity"=c(0,0, 0.35E6, 2.62E6, 9.07E6, 20.3E6, 34.71E6, 49.95E6, 64.27E6, 76.76E6, 87.15E6, 95.53E6, 102.15E6, 107.3E6, 111.27E6, 114.3E6, 116.61E6, 118.36E6, 119.68E6, 120.67E6, 123.234591E6))
  
  #Von-Bert
  Linf<-85.64   #L-Infinity (cm)
  k<-0.19       #Brody Growth Coefficient
  tnot<--0.39   #T-not
  Lt<-Linf*(1-exp(-k*(0:20-tnot)))
  
  #W-L Relationship
  a<-1.7E-5
  b<-3
  Wt<-a*Lt^b
  
  #Recruitment
  h<-0.99        #Steepness
  Peast<-0.23    #Proportion of recruits allocated to east of the Mississippi (per RS assessment)
  R0<-1.63E8      #Unfished recruitment
  SSB0<-4.72E15  #Unfished Spawning Biomass size (Eggs) for the entire GOM, not used 
  
  #Actual Red Snapper Data from Zach and Fabio
  if (Super_Comp==FALSE){
    load("C:/Users/nfisch/Documents/RSnapper_Data/ssdf_drop_2.RDATA", verbose=T) 
  } else if (Super_Comp==TRUE){
    load("/blue/edvcamp/nfisch/Spatial_Model/ssdf_drop_2.Rdata", verbose=T) 
  }

  #Depth Data
  if (Super_Comp==FALSE){
    Cells<-read.delim("C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/GOMFLA_Depth_n_substrate.txt", sep=" ",header=TRUE)
  } else if(Super_Comp==TRUE){
    Cells<-read.delim("/blue/edvcamp/nfisch/Spatial_Model/GOMFLA_Depth_n_substrate.txt", sep=" ",header=TRUE)
  }
  num_cells<-dim(Cells)[1]                                     #number of spatial cells
  
  ssdf_FLAsubset<-ssdf[ssdf$Lon > -87.5 & ssdf$Depth<70, ]
  ssdf_EastMississipi<-ssdf[ssdf$Lon > -89 & ssdf$Depth<70, ] #Subsetting dataset for east of the mississippi river to apportion recruitment as the number of FLA cells / number of cells east of mississippi
  
  #Proportion of recruitment that goes into Florida cells (rough approximation)
  PropR_FLA<-dim(ssdf_FLAsubset)[1]/dim(ssdf_EastMississipi)[1]
  R0_FLA<-R0*Peast*PropR_FLA   #Unfished recruitment for Florida
  lR0_FLA<-signif(log(R0_FLA), digits=6)  #log R0 florida 
  lxo<-c(1,cumprod(exp(-M$M))[1:(num_ages-1)])   #survivorship
  lxo[num_ages]<-lxo[num_ages]/(1-exp(-M$M[num_ages]))   #plus group survivorship
  N0_FLA_age<-exp(lR0_FLA)*lxo
  SSB0_FLA<-sum(N0_FLA_age*Fec$Fecundity)  #Unfished SSB calc
  
  ################################################
  #Creating Preference functions
  ################################################
  
  #Aggregate, mean numbers caught at depth and sd
  Agg_list<-list()
  for (i in 1:11){
    Agg_list[[i]]<-do.call(data.frame, aggregate(ssdf[,50+i],by=list(ssdf$Depth), FUN=function(x) c(Av=mean(x), sig=sd(x))))
  }
  
  #Fitting mean depth at age and variance using VBF, then assuming those follow a normal
  library(Hmisc)
  Means<-sapply(Agg_list, FUN=function(x){weighted.mean(x[,1], w=x[,2])})
  Vars<-sapply(Agg_list, FUN=function(x){wtd.var(x[,1], x[,2])})
  
  mod_means<-nls(Means ~ C*(1-exp(-r*(0:10-a))), start=list(C=80, r=1, a=1))
  Pred_Means<-predict(mod_means)[1:11]
  
  mod_vars<-nls(Vars ~ C*(1-exp(-r*(0:10-a))), start=list(C=1000, r=1, a=1))
  Pred_Vars<-predict(mod_vars)[1:11]
  
  #Preference Densities
  Pref_mat_depth<-matrix(NA, ncol=num_ages, nrow=ceiling(max(Cells$Depth)))  #11 is the number of ages 
  Pref_mat_depth_cumd<-matrix(NA, ncol=num_ages, nrow=ceiling(max(Cells$Depth)))
  for (i in 1:length(Pred_Means)){
    Pref_mat_depth[,i]<-dnorm(seq(1,ceiling(max(Cells$Depth))),mean=Pred_Means[i], sd=sqrt(Pred_Vars[i]))
    for (j in 1:ceiling(max(Cells$Depth))){#Using CDF
      Pref_mat_depth_cumd[j,i]<-pnorm(j+0.5, mean=Pred_Means[i], sd=sqrt(Pred_Vars[i]))-pnorm(j-0.5, mean=Pred_Means[i], sd=sqrt(Pred_Vars[i]))
    }
  }
  Pref_mat_depth[,12:num_ages]<-Pref_mat_depth[,11] #preference for older ages is the same
  
  #Making it such that every column has a max at 1
  Pref_mat_depth_standardized<-t(t(Pref_mat_depth)/apply(Pref_mat_depth,2,max))
  row.names(Pref_mat_depth_standardized)<-seq(1,ceiling(max(Cells$Depth)))
  
  #Preference for substrate type
  if (Super_Comp==FALSE){
    Pref_sub<-read.delim("C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/substrate_pref.txt", sep=" ",header=TRUE)
  } else if(Super_Comp==TRUE){
    Pref_sub<-read.delim("/blue/edvcamp/nfisch/Spatial_Model/substrate_pref.txt", sep=" ",header=TRUE)
  }
  Pref_sub[,12:num_ages]<-Pref_sub[,11] #preference for older ages is the same
  Pref_sub<-as.matrix(Pref_sub)
  
  ################################################
  #Subsetting dataframe for only grids in Florida
  ################################################

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
  
  #Reading in distance matrix
  if (Super_Comp==FALSE){
    dist_mat<-as.matrix(read.table(file="C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/dist_mat.txt"))
  } else if (Super_Comp==TRUE){
    dist_mat<-as.matrix(read.table("/blue/edvcamp/nfisch/Spatial_Model/dist_mat.txt"))
  }
  
  #Crude county Gulf coast midpoints (county order south to north) There are 22 along the gulf coast
  Gulf_County_Midpoints<-data.frame(County=c("Monroe","Collier","Lee","Charlotte","Sarasota","Manatee","Pinellas","Hillsborough","Pasco","Hernando","Citrus","Levy","Dixie","Taylor","Jefferson","Wakulla","Franklin","Gulf","Bay","Walton","Okaloosa","Santa Rosa","Escambia"),
                                    X=c(-80.94,-81.81,-82,-82.32,-82.53,-82.65,-82.85,-82.42,-82.73,-82.65,-82.64,-83.03,-83.31,-83.69,-84.02,-84.29,-84.68,-85.36,-85.72,-86.19,-86.59,-86.86,-87.19),
                                    Y=c(25.14,26.14,26.5,26.87,27.25,27.52,27.88,27.79,28.28,28.56,28.79,29.15,29.46,29.92,30.10,30.06,29.84,29.66,30.11,30.33,30.39,30.38,30.32),
                                    Pop=c(75027,378488,754610,184998,426718,394855,975280,1436888,539630,190865,147929,40770,16700,21623,14288,32461,11736,16164,185287,71375,207269,179349,315534))
  #Adding column that has proportion of population in total (will use as proxy for total effort)
  Gulf_County_Midpoints$Prop_pop<-Gulf_County_Midpoints$Pop/sum(Gulf_County_Midpoints$Pop)
  #Number of ports
  num_ports<-dim(Gulf_County_Midpoints)[1]                     #number of ports that effort goes out from
  
  #Reading in Distance from each county midpoint to each cell
  if (Super_Comp==FALSE){
    County_Distance<-as.matrix(read.table(file="C:/Users/nfisch/Documents/Snapper_Simulation_Assessments/County_Distance.txt"))
  } else if (Super_Comp==TRUE){
    County_Distance<-as.matrix(read.table("/blue/edvcamp/nfisch/Spatial_Model/County_Distance.txt"))
  }
  
  #Exponential Distance function based on data
  if (Super_Comp==FALSE){
    RSnap_Tag<-read.csv("C:/Users/nfisch/Documents/RSnapper_Data/Tag_snapper.csv")
  } else if (Super_Comp==TRUE){
    RSnap_Tag<-read.csv("/blue/edvcamp/nfisch/Spatial_Model/Tag_snapper.csv")
  }
  RSnap_Tag<-RSnap_Tag[!is.na(RSnap_Tag[,1]),] #Taking out NA data
  RSnap_Tag$Dist_Traveled<-earth.dist(RSnap_Tag[,"Lon1"], RSnap_Tag[,"Lat1"], RSnap_Tag[,"Lon2"], RSnap_Tag[,"Lat2"])  #Calculating the distance a red snapper traveled
  RSnap_Tag$Dist_Traveled_peryr<-RSnap_Tag$Dist_Traveled*(365/RSnap_Tag$days_at_large) #Calculating the distance they would have traveled in a full year based on their days at large
  func_move<-function(theta){
    lam_dist<-exp(theta[1])
    NLL<--1*sum(dexp(RSnap_Tag$Dist_Traveled_peryr[RSnap_Tag$days_at_large>200],lam_dist,log=T))
    return(NLL)
  }
  move_distcost<-optim(log(0.05), fn=func_move, method="BFGS")
  if (lam_data==TRUE){
    lam<-exp(move_distcost$par)
  } else {
    lam<-lam_moveCost
  }
  
  #Movement Transition Matrix, without threshold preference (used for initial year parameterization)
  p_move<-array(0,dim=c(num_cells,num_cells,num_ages))
  #Probability of moving from cell i to cell j at age k
  for (j in 1:num_ages){
    p_move[,,j]<- t(t(exp(-lam*dist_mat)^w_dist) * (Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]^w_depth * Pref_sub[Cells[,"substrate_code"],j]^w_substrate))
    p_move[,,j]<-p_move[,,j]/rowSums(p_move[,,j])
  }
  
  ####################
  #Starting population
  ####################
  N_wSpace_preM<-array(0, dim=c(nyear+1,num_ages,num_cells))    #Actual population-at-age in each cell for each year  
  N_wSpace_postM<-array(0, dim=c(nyear+1,num_ages,num_cells))    #Actual population-at-age in each cell for each year
  SSB_space<-array(0,dim=nyear+1)                          #Spawning Biomass
  ExploitBiomass_space<-matrix(0, nrow=nyear+1, ncol=num_cells)  #Exploitable Biomass (for fishing preference)
  ExploitAbun_space<-matrix(0, nrow=nyear+1, ncol=num_cells)    #Exploitable Abundance (for fishing preference)
  Catch_num_space<-matrix(0,nrow=nyear, ncol=num_cells)              #Total catch in each cell
  Catch_numage_space<-array(0,dim=c(nyear,num_ages,num_cells))       #Catch at age in each cell 
  Catch_bio_space<-matrix(0,nrow=nyear, ncol=num_cells)              #Total Catch in each cell (biomass) 
  F_space<-array(0, dim=c(nyear,num_ages,num_cells))           #Fishing mortality over Space
  p_move_wDD<-array(0,dim=c(num_cells,num_cells,num_ages))     #Movement Transition Matrix
  
  P_Fish<-array(0,dim=c(num_ports,num_cells,nyear))      #Matrix describing the probability of fishing a cell 
#  Effort_space<-matrix(0,nrow=nyear, ncol=num_cells)     #Effort expended in each cell 
  Effort_space<-matrix(c(rep(0,num_cells*50),rep(1,num_cells*100)),nrow=nyear, ncol=num_cells, byrow=T)     #Effort expended in each cell 
  Selvec<-1/(1+exp(-sel_grow*(seq(0,num_ages-1)-sel_midpt))) #Selectivity, could name parameters
  Sel<-matrix(Selvec,nrow=num_ages, ncol=num_cells) #Sel Matrix
  
  ##################################
  #Statewide Angler hours each year
  ##################################
  #Effort as logistic increase and linear decrease
  Effort_logis_mean<-eff_scalar/(1+exp(-eff_grate*(seq(1,nyear,1)-eff_midpt)))
  if(PE==FALSE){
    Effort<-c(rep(0,nyear_init),Effort_logis_mean[1:50],eff_scalar+(eff_scalar/2)-seq(eff_scalar/2+2,eff_scalar+1,length.out=50))
  } else if(PE==TRUE | PE=="Hybrid"){
#    Effort<-c(rep(0,nyear_init),rnorm(nyear-nyear_init, mean=c(Effort_logis_mean[1:50],eff_scalar+(eff_scalar/2)-seq(eff_scalar/2+2,eff_scalar+1,length.out=50)), sd=cv_totaleff*c(Effort_logis_mean[1:50],eff_scalar+(eff_scalar/2)-seq(eff_scalar/2+2,eff_scalar+1,length.out=50))))
#For constant effort after asymptote
    Effort<-c(rep(0,nyear_init),rnorm(nyear-nyear_init, mean=c(Effort_logis_mean[1:50],rep(eff_scalar*0.75,50)), sd=cv_totaleff*c(Effort_logis_mean[1:50],rep(eff_scalar*0.75,50))))
  }
  Effort<-ifelse(Effort<0, 1, Effort)  #Error catcher to make negative efforts zero
  
  Effort_midpoints<-matrix(0,nrow=nyear,ncol=num_ports)
  #Effort going from each port each year (i), could add process error to this
  if(PE==FALSE){
    for (i in 1:nyear){
      Effort_midpoints[i,]<-Gulf_County_Midpoints$Prop_pop*Effort[i]
    }
  } else if(PE==TRUE | PE=="Hybrid"){
    for (i in 1:nyear){
      Effort_midpoints[i,]<-rmultinom(n=1, size=Effort[i], prob=Gulf_County_Midpoints$Prop_pop)
    }
  }
  
  #1st Dimension is years (nyear)
  #2nd Dimension is Age (0-20)
  #3rd Dimension of array is spatial cell
  
  ###########################################################################################################
  #Spatial Model with no stochasticity in fish movement or fisher effort distribution, used for DGM scenario
  ###########################################################################################################
    if (PE=="Hybrid"){
    #Initial Year Recruitment    
    init_recs<-rlnorm(num_ages-1,meanlog=lR0_FLA,sdlog=sig_logR) #Recruitments to initialize model 
    init_recs[num_ages]<-exp(lR0_FLA) #plus group variation should smooth out over time 
    N_wSpace_preM[1,1,]<-N_wSpace_postM[1,1,]<-init_recs[1]*((Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1]))
    
    #Setting initial pop over space (just spread according to their depth preference)
    for (j in 2:num_ages){  #Each j age
      N_wSpace_preM[1,j,]<-init_recs[j]*lxo[j]*((Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]*Pref_sub[Cells[,"substrate_code"],j])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]*Pref_sub[Cells[,"substrate_code"],j]))
    }
    
    #Initial Movement Loop
    for(j in 2:num_ages){
      N_wSpace_postM[1,j,]<-N_wSpace_preM[1,j,] %*% p_move[,,j]  #vectorized
    }
    
    #Exploitable Biomass, etc. in first year (in each cell), Spawning Biomass (Fecundity)
    ExploitBiomass_space[1,]<-colSums(Sel*N_wSpace_postM[1,,]*Wt)     #Exploitable Biomass, vectorized
    ExploitAbun_space[1,]<-colSums(Sel*N_wSpace_postM[1,,])     #Exploitable Abundance, vectorized
    SSB_space[1]<-sum(rowSums(N_wSpace_preM[1,,])*Fec$Fecundity)
    
    #Projecting population forward from first year
    for (i in 2:(nyear+1)){ #Pop Loop (i year)
      
      #Effort Loop 
      for (p in 1:num_ports){ #Each port
        #There is a probability of fishing matrix each year because abundance changes
        P_Fish[p,,i-1]<-exp(-lam_Costdist*County_Distance[p,])^w_cost * (1/(1+exp(-Abunpref_grate*(ExploitBiomass_space[i-1,]-median(ExploitBiomass_space[i-1,])))))^w_profit
        P_Fish[p,,i-1]<-P_Fish[p,,i-1]/sum(P_Fish[p,,i-1]) #Standardizing so each row sums to 1
        Effort_space[i-1,]<-Effort_space[i-1,]+P_Fish[p,,i-1]*Effort_midpoints[i-1,p]   #Effort in each cell in each year, sum of effort coming from each port
      }
      
      #Fishing Mortality
      F_space[i-1,,]<-q*t(t(Sel)*Effort_space[i-1,])   #Fishing mortality in each cell, year by age by cell, vectorized
      
      #Mortality loop 
      N_wSpace_preM[i,2:num_ages,]<-N_wSpace_postM[i-1,1:(num_ages-1),]*exp(-(M$M[1:(num_ages-1)]+F_space[i-1,1:(num_ages-1),]))
      #Plus Group
      N_wSpace_preM[i,num_ages,]<-N_wSpace_preM[i,num_ages,]+N_wSpace_postM[i-1,num_ages,]*exp(-(M$M[num_ages]+F_space[i-1,num_ages,]))
      
      #Spawning Biomass for SR function
      SSB_space[i]<-sum(rowSums(N_wSpace_preM[i,,])*Fec$Fecundity)
      #Recruitment in each year, since its age 0, same year index for SSB. Age zeros don't move
      N_wSpace_preM[i,1,]<-N_wSpace_postM[i,1,]<-rlnorm(n=1,meanlog=log((4*h*exp(lR0_FLA)*SSB_space[i])/(SSB0_FLA*(1-h)+SSB_space[i]*(5*h-1))), sdlog=sig_logR)*((Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1]))

      #Movement matrix, probability of moving from cell a to cell b at age j
      #It will be different for each year, however currently not saving it because memory allocation is too large
      ind_vec<-ifelse(colSums(t(t(N_wSpace_preM[i,,])*Lt^2)) > quantile(colSums(t(t(N_wSpace_preM[1,,])*Lt^2)),DD_Thresh_quant), 1, 0) #indicator vector for threshold
      for (j in 1:num_ages){
        #Threshold Preference function (which functions subset "b" cells for threshold condition)
        p_move_wDD[,,j]<-t(t(exp(-lam*dist_mat)^w_dist) * (Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]^w_depth * Pref_sub[Cells[,"substrate_code"],j]^w_substrate * ((1/(colSums(t(t(N_wSpace_preM[i,,])*Lt^2))/quantile(colSums(t(t(N_wSpace_preM[1,,])*Lt^2)),DD_Thresh_quant))^DD_rate)^w_DD)^ind_vec))
        p_move_wDD[,,j]<-p_move_wDD[,,j]/rowSums(p_move_wDD[,,j])
      }      
      
      #Movement loop, Age zeros excluded because they do not move... this is the one that can be moved to matrix algebra I think
      for(j in 2:num_ages){
        N_wSpace_postM[i,j,]<-N_wSpace_preM[i,j,] %*% p_move_wDD[,,j]
      }
      
      ExploitBiomass_space[i,]<-colSums(Sel*N_wSpace_postM[i,,]*Wt)     #Exploitable Biomass
      ExploitAbun_space[i,]<-colSums(Sel*N_wSpace_postM[i,,])           #Exploitable abundance
      Catch_bio_space[i-1,]<-colSums((F_space[i-1,,]/(F_space[i-1,,]+M$M))*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,])))*Wt)
      Catch_num_space[i-1,]<-colSums((F_space[i-1,,]/(F_space[i-1,,]+M$M))*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,]))))
      Catch_numage_space[i-1,,]<-F_space[i-1,,]/(F_space[i-1,,]+M$M)*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,])))
      
    }#End of Pop Loop 
    ##################################################################################################
  }
  #######################################
  #Spatial Model with no process error
  ####################################### 
  else if (PE==FALSE){
    #Initial Year Recruitment
    N_wSpace_preM[1,1,]<-N_wSpace_postM[1,1,]<-exp(lR0_FLA)*((Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1]))
    #Setting initial pop over space (according to depth and substrate preference)
    N_wSpace_preM[1,,]<-N0_FLA_age*(t(Pref_mat_depth_standardized[round(Cells[,"Depth"]),]*Pref_sub[Cells[,"substrate_code"],])/colSums(Pref_mat_depth_standardized[round(Cells[,"Depth"]),]*Pref_sub[Cells[,"substrate_code"],]))

    #Initial Movement Loop
    for(j in 2:num_ages){
      N_wSpace_postM[1,j,]<-N_wSpace_preM[1,j,] %*% p_move[,,j]  #vectorized
    }
    
    #Exploitable Biomass, etc. in first year (in each cell), Spawning Biomass (Fecundity)
    ExploitBiomass_space[1,]<-colSums(Sel*N_wSpace_postM[1,,]*Wt)     #Exploitable Biomass, vectorized
    ExploitAbun_space[1,]<-colSums(Sel*N_wSpace_postM[1,,])     #Exploitable Abundance, vectorized
    SSB_space[1]<-sum(rowSums(N_wSpace_preM[1,,])*Fec$Fecundity)
    
    #Projecting population forward from first year
    for (i in 2:(nyear+1)){ #Pop Loop (i year)
      
      #Effort Loop 
      for (p in 1:num_ports){ #Each port
        #There is a probability of fishing matrix each year because abundance changes
        P_Fish[p,,i-1]<-exp(-lam_Costdist*County_Distance[p,])^w_cost * (1/(1+exp(-Abunpref_grate*(ExploitBiomass_space[i-1,]-median(ExploitBiomass_space[i-1,])))))^w_profit
        P_Fish[p,,i-1]<-P_Fish[p,,i-1]/sum(P_Fish[p,,i-1]) #Standardizing so each row sums to 1
        Effort_space[i-1,]<-Effort_space[i-1,]+P_Fish[p,,i-1]*Effort_midpoints[i-1,p]   #Effort in each cell in each year, sum of effort coming from each port
      }
      
      #Fishing Mortality
      F_space[i-1,,]<-q*t(t(Sel)*Effort_space[i-1,])   #Fishing mortality in each cell, year by age by cell, vectorized
      
      #Mortality loop 
      N_wSpace_preM[i,2:num_ages,]<-N_wSpace_postM[i-1,1:(num_ages-1),]*exp(-(M$M[1:(num_ages-1)]+F_space[i-1,1:(num_ages-1),]))
      #Plus Group
      N_wSpace_preM[i,num_ages,]<-N_wSpace_preM[i,num_ages,]+N_wSpace_postM[i-1,num_ages,]*exp(-(M$M[num_ages]+F_space[i-1,num_ages,]))
      
      #Spawning Biomass for SR function
      SSB_space[i]<-sum(rowSums(N_wSpace_preM[i,,])*Fec$Fecundity)
      #Recruitment in each year, since its age 0, same year index for SSB. Age zeros don't move
      N_wSpace_preM[i,1,]<-N_wSpace_postM[i,1,]<-((4*h*exp(lR0_FLA)*SSB_space[i])/(SSB0_FLA*(1-h)+SSB_space[i]*(5*h-1)))*((Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1]))
      
      #Movement matrix, probability of moving from cell a to cell b at age j
      #It will be different for each year, however currently not saving it because memory allocation is too large
      ind_vec<-ifelse(colSums(t(t(N_wSpace_preM[i,,])*Lt^2)) > quantile(colSums(t(t(N_wSpace_preM[1,,])*Lt^2)),DD_Thresh_quant), 1, 0) #indicator vector for threshold
      for (j in 1:num_ages){
        #Threshold Preference function (which functions subset "b" cells for threshold condition)
        p_move_wDD[,,j]<-t(t(exp(-lam*dist_mat)^w_dist) * (Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]^w_depth * Pref_sub[Cells[,"substrate_code"],j]^w_substrate * ((1/(colSums(t(t(N_wSpace_preM[i,,])*Lt^2))/quantile(colSums(t(t(N_wSpace_preM[1,,])*Lt^2)),DD_Thresh_quant))^DD_rate)^w_DD)^ind_vec))
        p_move_wDD[,,j]<-p_move_wDD[,,j]/rowSums(p_move_wDD[,,j])
      }      
      
      #Movement loop, Age zeros excluded because they do not move... this is the one that can be moved to matrix algebra I think
      for(j in 2:num_ages){
        N_wSpace_postM[i,j,]<-N_wSpace_preM[i,j,] %*% p_move_wDD[,,j]
      }
      
      ExploitBiomass_space[i,]<-colSums(Sel*N_wSpace_postM[i,,]*Wt)     #Exploitable Biomass
      ExploitAbun_space[i,]<-colSums(Sel*N_wSpace_postM[i,,])           #Exploitable abundance
      Catch_bio_space[i-1,]<-colSums((F_space[i-1,,]/(F_space[i-1,,]+M$M))*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,])))*Wt)
      Catch_num_space[i-1,]<-colSums((F_space[i-1,,]/(F_space[i-1,,]+M$M))*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,]))))
      Catch_numage_space[i-1,,]<-F_space[i-1,,]/(F_space[i-1,,]+M$M)*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,])))
      
    }#End of Pop Loop 
    ##################################################################################################
  }
  #######################################
  #Spatial Model with Process Error
  ####################################### 
  else if(PE==TRUE){
    #Initial Year Recruitment
    init_recs<-rlnorm(num_ages-1,meanlog=lR0_FLA,sdlog=sig_logR) #Recruitments to initialize model 
    init_recs[num_ages]<-exp(lR0_FLA) #plus group variation should smooth out over time 
    N_wSpace_preM[1,1,]<-N_wSpace_postM[1,1,]<-rmultinom(n=1,size=round(init_recs[1]),prob=(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1]))
    #Setting initial pop over space (just spread according to their depth preference)
    for (j in 2:num_ages){  #Each j age
      N_wSpace_preM[1,j,]<-rmultinom(n=1,size=round(init_recs[j]*lxo[j]),prob=((Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]*Pref_sub[Cells[,"substrate_code"],j])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]*Pref_sub[Cells[,"substrate_code"],j])))
    }
    
    #Initial Movement Loop
    for (j in 2:num_ages){  
      N_wSpace_postM[1,j,]<-colSums( t( sapply(1:num_cells, function(b) rmultinom(n=1, size=round(N_wSpace_preM[1,j,b]), prob=p_move[b,,j])) ) )
    }
    
    #Exploitable Biomass, etc. in first year (in each cell), Spawning Biomass (Fecundity)
    ExploitBiomass_space[1,]<-colSums(Sel*N_wSpace_postM[1,,]*Wt)     #Exploitable Biomass, vectorized
    ExploitAbun_space[1,]<-colSums(Sel*N_wSpace_postM[1,,])     #Exploitable Abundance, vectorized
    SSB_space[1]<-sum(rowSums(N_wSpace_preM[1,,])*Fec$Fecundity)
    
    #Projecting population forward from first year
    for (i in 2:(nyear+1)){ #Pop Loop (i year)
      
      for (p in 1:num_ports){ #Each port
        #There is a probability of fishing matrix each year because abundance changes
        P_Fish[p,,i-1]<-exp(-lam_Costdist*County_Distance[p,])^w_cost * (1/(1+exp(-Abunpref_grate*(ExploitBiomass_space[i-1,]-median(ExploitBiomass_space[i-1,])))))^w_profit
        P_Fish[p,,i-1]<-P_Fish[p,,i-1]/sum(P_Fish[p,,i-1]) #Standardizing so each row sums to 1
        Effort_space[i-1,]<-Effort_space[i-1,]+rmultinom(n=1,size=Effort_midpoints[i-1,p], prob=P_Fish[p,,i-1])   #Effort in each cell in each year, sum of effort coming from each port
      }
      
      #Fishing Mortality
      F_space[i-1,,]<-q*t(t(Sel)*Effort_space[i-1,])   #Fishing mortality in each cell, year by age by cell, vectorized
      
      #Mortality Loop
      N_wSpace_preM[i,2:num_ages,]<-N_wSpace_postM[i-1,1:(num_ages-1),]*exp(-(M$M[1:(num_ages-1)]+F_space[i-1,1:(num_ages-1),]))
      #Plus group
      N_wSpace_preM[i,num_ages,]<-N_wSpace_preM[i,num_ages,]+N_wSpace_postM[i-1,num_ages,]*exp(-(M$M[num_ages]+F_space[i-1,num_ages,])) 
      
      #Spawning Biomass for SR function
      SSB_space[i]<-sum(rowSums(N_wSpace_preM[i,,])*Fec$Fecundity)
      #Recruitment in each year, since its age 0, same year index for SSB. Age zeros don't move
      N_wSpace_preM[i,1,]<-N_wSpace_postM[i,1,]<-rmultinom(n=1, size=rlnorm(n=1,meanlog=log((4*h*exp(lR0_FLA)*SSB_space[i])/(SSB0_FLA*(1-h)+SSB_space[i]*(5*h-1))), sdlog=sig_logR),prob=(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1])/sum(Pref_mat_depth_standardized[round(Cells[,"Depth"]),1]*Pref_sub[Cells[,"substrate_code"],1]))
      
      #Movement matrix, probability of moving from cell a to cell b at age j
      #It will be different for each year, however currently not saving it because memory allocation is too large
      ind_vec<-ifelse(colSums(t(t(N_wSpace_preM[i,,])*Lt^2)) > quantile(colSums(t(t(N_wSpace_preM[1,,])*Lt^2)),DD_Thresh_quant), 1, 0) #indicator vector for threshold
      for (j in 1:num_ages){
        #Threshold Preference function (which functions subset "b" cells for threshold condition)
        p_move_wDD[,,j]<-t(t(exp(-lam*dist_mat)^w_dist) * (Pref_mat_depth_standardized[round(Cells[,"Depth"]),j]^w_depth * Pref_sub[Cells[,"substrate_code"],j]^w_substrate * ((1/(colSums(t(t(N_wSpace_preM[i,,])*Lt^2))/quantile(colSums(t(t(N_wSpace_preM[1,,])*Lt^2)),DD_Thresh_quant))^DD_rate)^w_DD)^ind_vec))
        p_move_wDD[,,j]<-p_move_wDD[,,j]/rowSums(p_move_wDD[,,j]) #standardizing
      }
      
      #Movement loop, Age zeros excluded because they do not move
      for (j in 2:num_ages){  
        N_wSpace_postM[i,j,]<-colSums( t( sapply(1:num_cells, function(b) rmultinom(n=1, size=round(N_wSpace_preM[i,j,b]), prob=p_move_wDD[b,,j])) ) )
      }
      
      ExploitBiomass_space[i,]<-colSums(Sel*N_wSpace_postM[i,,]*Wt)     #Exploitable Biomass
      ExploitAbun_space[i,]<-colSums(Sel*N_wSpace_postM[i,,])           #Exploitable abundance
      Catch_bio_space[i-1,]<-colSums((F_space[i-1,,]/(F_space[i-1,,]+M$M))*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,])))*Wt)
      Catch_num_space[i-1,]<-colSums((F_space[i-1,,]/(F_space[i-1,,]+M$M))*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,]))))
      Catch_numage_space[i-1,,]<-F_space[i-1,,]/(F_space[i-1,,]+M$M)*N_wSpace_postM[i-1,,]*(1-exp(-(M$M+F_space[i-1,,])))
      
    }#End of Pop Loop 
  }
  
  #############################
  #Getting Data for dat file
  #############################
  #Sampling Catch, 5% CV
  Catch_obs<-rnorm(nyear, mean=rowSums(Catch_bio_space), sd=cv_harv*rowSums(Catch_bio_space)) #CV is 5%, making sd of catch 
  
  #CPUE data
  CPUE_true<-apply(Catch_bio_space[51:150,]/Effort_space[51:150,],1, sum, na.rm=T)
  CPUE_obs<-rnorm(nyear-nyear_init, CPUE_true, sd=cv_cpue*CPUE_true)

  #Observation error on Effort with a CV on 5%
  Effort_obs<-rnorm(nyear,mean=rowSums(Effort_space), sd=cv_effort*rowSums(Effort_space))
  
  #Composition sampling
  
  #Percentage of Effort sampled (currently 5%), divied up by each port. This assumes sampling is representative of port distribution
  #Effort_sample<-Effort_midpoints*Perc_eff_smpld
  Effort_sample<-matrix(0, nrow=nyear,ncol=num_ports)
  for (i in (nyear_init+1):nyear){
    if(sum(Effort_space[i,])>0){
     Effort_sample[i,]<-rmultinom(n=1, size=sum(Effort_space[i,])*Perc_eff_smpld, prob=Effort_midpoints[i,])  #adding in error to which ports are sampled
    }
  }
  
  #Choosing which cells to sample (matrix is filled in with cell indicator to be sampled with 1 unit of effort)
  Cells_sampled<-array(0, dim=c(nyear, num_ports, ceiling(max(Effort_sample))))
  #Bio_Samples contains year, ports, sample, and ages
  Bio_samples<-array(NA,dim=c(nyear,num_ports,ceiling(max(Effort_sample)),num_ages))
  #Sampling the Fishery
  for (i in 1:nyear){
    for (p in 1:num_ports){
      if (Effort_sample[i,p]!=0){  #If there are actual units of effort sampled from that port, then which cells do they go to
        Cells_sampled[i,p,1:ceiling(Effort_sample[i,p])]<-sample(1:num_cells,size=ceiling(Effort_sample[i,p]),prob=P_Fish[p,,i], replace=TRUE)
        for (k in 1:ceiling(Effort_sample[i,p])){
          if (Prop_sample==TRUE){
            if (Effort_space[i,Cells_sampled[i,p,k]]>0){ #If the effort going into that cell was higher than zero, can sample it
              #The number of samples taken number of fish caught per unit of effort that went to that cell, times the proportion sampled per unit of effort  
              if (sum(Catch_numage_space[i,,Cells_sampled[i,p,k]])>0){ #if there is no catch then can't sample
                Bio_samples[i,p,k,]<-rmultinom(n=1, size=round(sum(Catch_numage_space[i,,Cells_sampled[i,p,k]])/Effort_space[i,Cells_sampled[i,p,k]]*Perc_sample_pb), prob=Catch_numage_space[i,,Cells_sampled[i,p,k]]/sum(Catch_numage_space[i,,Cells_sampled[i,p,k]]))
              }
            }
          }
          else if (Prop_sample==FALSE){ #If we are not proportionally sampling
            if (sum(Catch_numage_space[i,,Cells_sampled[i,p,k]])>Sample_size_pb){ 
              #Lets say they take a set number of fish from each trip, Sample_size_pb
              Bio_samples[i,p,k,]<-rmultinom(n=1, size=Sample_size_pb, prob=Catch_numage_space[i,,Cells_sampled[i,p,k]]/sum(Catch_numage_space[i,,Cells_sampled[i,p,k]]))
            }
            # If the number of fish caught in a cell is less than the sample size take per boat, sample all the fish caught in that cell
            else if (sum(Catch_numage_space[i,,Cells_sampled[i,p,k]]) < Sample_size_pb & sum(Catch_numage_space[i,,Cells_sampled[i,p,k]]) > 0){
              Bio_samples[i,p,k,]<-rmultinom(n=1, size=round(sum(Catch_numage_space[i,,Cells_sampled[i,p,k]])), prob=Catch_numage_space[i,,Cells_sampled[i,p,k]]/sum(Catch_numage_space[i,,Cells_sampled[i,p,k]]))
            }
          }
        }
      }
    }
  }
  
  #Total sample for each year
  Bio_sample_yr<-apply(X=Bio_samples, MARGIN=c(1,4), FUN=sum, na.rm=T)
  
  #Now getting the data that would feed into assessment, the above but in composition form
  Age_Comp_Data<-Bio_sample_yr/rowSums(Bio_sample_yr)    #Data
  Age_Comp_Data<-ifelse(Age_Comp_Data==0,Age_Comp_Data+1E-5,Age_Comp_Data)  #Suppressing zeroes by adding a small constant
  Age_Comp_Data<-Age_Comp_Data/rowSums(Age_Comp_Data) #renormalizing
    
  #############################################
  #Adding FIM index CPUE and Catch Composition
  #############################################
  
  #FIM CPUE (only for latter 3/5ths of time series)
  FIM_years<-(nyear-nyear_init)*Perc_yrs_FIM
  
  #Start with Selectivity the same as the fishery 
  Sel_FIM<-Selvec
  
  #Randomly sample cells for bio data and FI-index
  FIM_sample<-matrix(0,nrow=FIM_years, ncol=num_FIM_samples)
  FIM_index_data<-vector(length=FIM_years)
  FIM_biosample<-array(0,dim=c(FIM_years,num_FIM_samples,num_ages)) 
  for (i in 1:FIM_years){
    FIM_sample[i,]<-sample(1:num_cells, size=num_FIM_samples, replace=TRUE) #Getting sampled cells each year
    for (b in 1:num_FIM_samples){
      if (sum(N_wSpace_postM[(i+nyear-FIM_years),,FIM_sample[i,b]])>0){ #If there are fish in a cell, then sample it
        FIM_biosample[i,b,]<-rmultinom(n=1, size=round(sum(Sel_FIM*N_wSpace_postM[(i+nyear-FIM_years),,FIM_sample[i,b]]*FIM_q)), prob=(Sel_FIM*N_wSpace_postM[(i+nyear-FIM_years),,FIM_sample[i,b]])/sum(Sel_FIM*N_wSpace_postM[(i+nyear-FIM_years),,FIM_sample[i,b]]))
      } 
    }
    FIM_index_data[i]<-sum(FIM_biosample[i,,])/num_FIM_samples  #Survey index is just the total catch over the number of samples
  }
  
  FIM_biodata<-t(apply(FIM_biosample,1,FUN=function(x) apply(x,2,sum)/sum(x)))
  FIM_biodata<-ifelse(FIM_biodata==0,FIM_biodata+1E-5,FIM_biodata)  #Suppressing zeroes by adding a small constant
  FIM_biodata<-FIM_biodata/rowSums(FIM_biodata) #renormalizing
  
  dir.create(save_wd)
  
  write(c(seed,PE,nyear,nyear_init,num_ages,lam_data,lam_moveCost,DD_rate,DD_Thresh_quant,sig_logR,q, sel_grow, sel_midpt, lam_Costdist,
          Abunpref_grate, eff_scalar, eff_midpt, eff_grate, cv_harv, cv_cpue, cv_effort, Perc_eff_smpld, Prop_sample, Perc_sample_pb,
          Sample_size_pb, Perc_yrs_FIM,FIM_q, num_FIM_samples, w_dist, w_depth, w_DD, w_cost, w_profit),
        file = paste0(save_wd, "/Parameters.txt"), ncolumns=32)

  #Writing Dat File
  dat_file<-list(1, nyear-nyear_init, 0, num_ages-1, M_vec, 
                 Catch_obs[(nyear_init+1):nyear],                 #Observed Catch
#                 Effort_obs[(nyear_init+1):nyear],                #Observed Effort
                 CPUE_obs,                                        #Observed CPUE
                 t(Age_Comp_Data[(nyear_init+1):nyear,]),             #Fishery Age Composition data
                 rowSums(Bio_sample_yr[(nyear_init+1):nyear,]),    #SS for Fishery composition
                 c(rep(0,nyear-nyear_init-(nyear-nyear_init)*Perc_yrs_FIM),FIM_index_data),        #Survey Index
                 t(rbind(matrix(0,nrow=(nyear-nyear_init)*(1-Perc_yrs_FIM), ncol=num_ages),FIM_biodata)),    #Survey Composition
                 c(rep(0,nyear-nyear_init-(nyear-nyear_init)*Perc_yrs_FIM),apply(FIM_biosample,1,FUN = sum)),  #Sample size of FIM Composition
                 Fec$Fecundity,                                    #Fecundity
                 Wt,                                               #Wt-at-age
                 c(11,22,33))                                      #test Vector
  
  lapply(dat_file, write, paste0(save_wd,"/dat_file.dat"), append=TRUE, ncolumns=21)
  
  #Saving Ouput from model
  saveRDS(N_wSpace_postM, file=paste0(save_wd, "/N_wSpace_postM.rds" ))
#  saveRDS(N_wSpace_preM, file=paste0(save_wd, "/N_wSpace_preM.rds" ))
  saveRDS(SSB_space, file=paste0(save_wd, "/SSB_space.rds" ))
#  saveRDS(ExploitBiomass_space, file=paste0(save_wd, "/ExploitBiomass_space.rds" ))
  saveRDS(Catch_bio_space, file=paste0(save_wd, "/Catch_bio_space.rds" ))
  saveRDS(Catch_num_space, file=paste0(save_wd, "/Catch_num_space.rds" ))
  saveRDS(Catch_numage_space, file=paste0(save_wd, "/Catch_numage_space.rds" ))
  saveRDS(F_space, file=paste0(save_wd, "/F_space.rds" ))
  saveRDS(Effort_space, file=paste0(save_wd, "/Effort_space.rds" ))
  saveRDS(Sel, file=paste0(save_wd, "/Selectivity.rds" ))
  saveRDS(Sel_FIM, file=paste0(save_wd, "/Sel_FIM.rds" ))
#  saveRDS(P_Fish, file=paste0(save_wd, "/P_Fish.rds" ))
  #Data section below
  saveRDS(Effort_midpoints, file=paste0(save_wd, "/Effort_midpoints.rds" ))
  saveRDS(Catch_obs, file=paste0(save_wd, "/Catch_obs.rds" ))
  saveRDS(Effort_obs, file=paste0(save_wd, "/Effort_obs.rds" ))
  saveRDS(CPUE_obs, file=paste0(save_wd, "/CPUE_obs.rds" ))
  saveRDS(Cells_sampled, file=paste0(save_wd, "/Cells_sampled.rds" ))
  saveRDS(Bio_samples, file=paste0(save_wd, "/Bio_samples.rds" ))
  saveRDS(Bio_sample_yr, file=paste0(save_wd, "/Bio_sample_yr.rds" ))
  saveRDS(Age_Comp_Data, file=paste0(save_wd, "/Age_Comp_Data.rds" ))
  saveRDS(FIM_sample, file=paste0(save_wd, "/FIM_sample.rds" ))
  saveRDS(FIM_biosample, file=paste0(save_wd, "/FIM_biosample.rds" ))
  saveRDS(FIM_index_data, file=paste0(save_wd, "/FIM_index_data.rds" ))
  saveRDS(FIM_biodata, file=paste0(save_wd, "/FIM_biodata.rds" ))
}

#Status Quo Model
for (i in 101:1000){
Spatial_Model(save_wd=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/SQ_PE_GE_00540_FQ_",i),
              seed=i,PE=TRUE,nyear=150,nyear_init=50,num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0.03, Abunpref_grate=0.00025, eff_scalar=1e5, cv_totaleff=0.25, eff_midpt=25, eff_grate=0.15,
              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.005, Prop_sample=TRUE, Perc_sample_pb=0.40, Sample_size_pb=20, Perc_yrs_FIM=0.6,
              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1, Super_Comp=TRUE)          
}

#Random Fishing Model
#for (i in 1:1000){
#Spatial_Model(save_wd=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/RF_GE_1010_",i), 
#              seed=i,PE=FALSE,nyear=150,nyear_init=50,num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3,
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0, Abunpref_grate=0, eff_scalar=1e5, cv_totaleff=0.25, eff_midpt=25, eff_grate=0.15,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.1, Prop_sample=TRUE, Perc_sample_pb=0.10, Sample_size_pb=20, Perc_yrs_FIM=0.6, 
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1, Super_Comp=TRUE) 
#}

#Hybrid Fishing Model
#for (i in 906:1000){
#Spatial_Model(save_wd=paste0("/blue/edvcamp/nfisch/Spatial_Model/OM_Runs/RF_Hybrid_GE_0210_FQ_",i), 
#              seed=i,PE="Hybrid",nyear=150,nyear_init=50,num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3,
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0, Abunpref_grate=0, eff_scalar=1e5, cv_totaleff=0.25, eff_midpt=25, eff_grate=0.15,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.02, Prop_sample=TRUE, Perc_sample_pb=0.10, Sample_size_pb=20, Perc_yrs_FIM=0.6, 
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1, Super_Comp=TRUE) 
#}

warnings()