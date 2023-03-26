
##############################################
#Sampling Model for Spatially Explicit OM
##############################################

rm(list=ls(all=TRUE))

Spatial_Model_data<-function(save_wd,               #Working directory for saving output and data
                        extract_wd, 
                        seed,                  #random number seed
                        PE,                    #Process Error (turned on or off via T/F)
                        nyear,                 #Number of years in simulation, first 50 years will be unfished 
                        nyear_init,            #Number of unfished years to start simulation    
                        num_ages,              #Number of ages
                        num_cells,
                        num_ports,
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
                        cv_harv,               #CV of sampled harvest
                        cv_cpue,               #CV of CPUE
                        cv_effort,             #CV of sampled effort
                        Perc_eff_smpld,        #Percentage of effort sampled
                        Prop_sample,           #Are we proportionally sampling fish from the fishery catches or not
                        Perc_sample_pb,        #Percentage of samples taken per boat 
                        FIM_q,                 #FIM q
                        num_FIM_samples,       #Number of surveys
                        w_dist,                #power weight for fish movement distance 
                        w_depth,               #power weight for depth preference
                        w_substrate,           #power weight for substrate preference
                        w_DD,                  #power weight for DD preference
                        w_cost,                #power weight for fisher cost function
                        w_profit)              #power weight for fisher profit function
{
  set.seed(seed)
  #Von-Bert
  Linf<-85.64   #L-Infinity (cm)
  k<-0.19       #Brody Growth Coefficient
  tnot<--0.39   #T-not
  Lt<-Linf*(1-exp(-k*(0:20-tnot)))
  M_vec<-c(2,1.2,(0.099*Lt[8])/Lt[3:21])

  #W-L Relationship
  a<-1.7E-5
  b<-3
  Wt<-a*Lt^b
 
  setwd(paste0("where you saved spatial OM output/",extract_wd))
  N_wSpace_postM<-readRDS("N_wSpace_postM.rds")
  Effort_midpoints<-readRDS("Effort_midpoints.rds")
  Catch_bio_space<-readRDS("Catch_bio_space.rds")
  Catch_num_space<-readRDS("Catch_num_space.rds")
  Effort_space<-readRDS("Effort_space.rds")
  Catch_numage_space<-readRDS("Catch_numage_space.rds")
  Sel<-readRDS("Selectivity.rds")
  F_space<-readRDS("F_space.rds")
  SSB_space<-readRDS("SSB_space.rds")
 
  County_Distance<-as.matrix(read.table("C:/Users/nfisch/Documents/GitHub/Spatial_RedSnapper_OM_and_SM/County_Distance.txt"))
  Fec<-data.frame("Age"=0:20, "Fecundity"=c(0,0, 0.35E6, 2.62E6, 9.07E6, 20.3E6, 34.71E6, 49.95E6, 64.27E6, 76.76E6, 87.15E6, 95.53E6, 102.15E6, 107.3E6, 111.27E6, 114.3E6, 116.61E6, 118.36E6, 119.68E6, 120.67E6, 123.234591E6))

  P_Fish<-array(0,dim=c(num_ports,num_cells,nyear))      #Matrix describing the probability of fishing a cell 
  ExploitBiomass_space<-matrix(0, nrow=nyear+1, ncol=num_cells)  #Exploitable Biomass (for fishing preference)
  ExploitBiomass_space[1,]<-colSums(Sel*N_wSpace_postM[1,,]*Wt)     #Exploitable Biomass, vectorized
  for (i in 2:(nyear+1)){ #Pop Loop (i year)     
    for (p in 1:num_ports){ #Each port
      #There is a probability of fishing matrix each year because abundance changes
      P_Fish[p,,i-1]<-exp(-lam_Costdist*County_Distance[p,])^w_cost * (1/(1+exp(-Abunpref_grate*(ExploitBiomass_space[i-1,]-median(ExploitBiomass_space[i-1,])))))^w_profit
      P_Fish[p,,i-1]<-P_Fish[p,,i-1]/sum(P_Fish[p,,i-1]) #Standardizing so each row sums to 1
    }
    ExploitBiomass_space[i,]<-colSums(Sel*N_wSpace_postM[i,,]*Wt)     #Exploitable Biomass
  }
  #############################
  #Getting Data for dat file
  #############################
  #Sampling Catch, 5% CV
  Catch_obs<-rnorm(nyear, mean=rowSums(Catch_bio_space), sd=cv_harv*rowSums(Catch_bio_space)) #CV is 5%, making sd of catch 
  
  #Nominal CPUE data, need to standardize if you want CPUE proportional to effort
  Index_nom_true<-rowSums(Catch_bio_space[(nyear_init+1):nyear,])/rowSums(Effort_space[(nyear_init+1):nyear,])
  Index_nom_obs<-rnorm(nyear-nyear_init, Index_nom_true, sd=cv_cpue*Index_nom_true)

  #Composition sampling
  
  #Percentage of Effort sampled, divied up by each port. This assumes sampling is representative of port distribution
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
  
  #################################################
  #Adding Survey index CPUE and Composition
  #################################################
  
  #Survey Selectivity the same as the fishery (unless the fishery contact selectivity is dome shaped)
  Sel_FIM<-1/(1+exp(-sel_grow*(seq(0,num_ages-1)-sel_midpt)))
  
  #Randomly sample cells for bio data and FI-index
  FIM_sample<-matrix(0,nrow=nyear, ncol=num_FIM_samples)
  FIM_index_data<-vector(length=nyear)
  FIM_biosample<-array(0,dim=c(nyear,num_FIM_samples,num_ages)) 
  for (i in (nyear_init+1):nyear){
    FIM_sample[i,]<-sample(1:num_cells, size=num_FIM_samples, replace=TRUE) #Getting sampled cells each year
    for (b in 1:num_FIM_samples){
      if (sum(N_wSpace_postM[i,,FIM_sample[i,b]])>0){ #If there are fish in a cell, then sample it
        FIM_biosample[i,b,]<-rmultinom(n=1, size=round(sum(Sel_FIM*N_wSpace_postM[i,,FIM_sample[i,b]]*FIM_q)), prob=(Sel_FIM*N_wSpace_postM[i,,FIM_sample[i,b]])/sum(Sel_FIM*N_wSpace_postM[i,,FIM_sample[i,b]]))
      } 
    }
    FIM_index_data[i]<-sum(FIM_biosample[i,,])/num_FIM_samples  #Survey index is just the total catch over the number of samples
  }
  
  FIM_biodata<-t(apply(FIM_biosample,1,FUN=function(x) apply(x,2,sum)/sum(x)))
  FIM_biodata<-ifelse(FIM_biodata==0,FIM_biodata+1E-5,FIM_biodata)  #Suppressing zeroes by adding a small constant
  FIM_biodata<-FIM_biodata/rowSums(FIM_biodata) #renormalizing
  
  #RANDOMLY sampling catch at age in each year
  catch_age<-rowSums(Catch_numage_space, dims=2)
  catch_age_smpld<-matrix(NA, nrow=nyear-nyear_init, ncol=num_ages)
  for (j in (nyear_init+1):nyear){
    catch_age_smpld[j-nyear_init,]<-rmultinom(n=1, size=rowSums(Bio_sample_yr)[j], prob=catch_age[j,]) 
  }
  rnd_age_comp<-catch_age_smpld/ifelse(rowSums(catch_age_smpld)==0,1,rowSums(catch_age_smpld))
  rnd_age_comp<-ifelse(rnd_age_comp==0,rnd_age_comp+1E-5,rnd_age_comp)  #Suppressing zeroes by adding a small constant
  rnd_age_comp<-rnd_age_comp/rowSums(rnd_age_comp) #renormalizing
        
  #RANDOMLY sampling survey in each year
  surv_cage<-t(t(rowSums(N_wSpace_postM, dims=2))*Sel_FIM)
  surv_cage_smpl<-matrix(NA, nrow=nyear-nyear_init, ncol=num_ages)
  for (k in (nyear_init+1):nyear){
    surv_cage_smpl[k-nyear_init,]<-rmultinom(n=1, size=apply(FIM_biosample,1,FUN = sum)[k], prob=surv_cage[k,])
  }
  rnd_surv_comp<-surv_cage_smpl/rowSums(surv_cage_smpl)
  rnd_surv_comp<-ifelse(rnd_surv_comp==0,rnd_surv_comp+1E-5,rnd_surv_comp)  #Suppressing zeroes by adding a small constant
  rnd_surv_comp<-rnd_surv_comp/rowSums(rnd_surv_comp) #renormalizing
  
  dir.create(save_wd)
  
  write(c(seed,PE,nyear,nyear_init,num_ages,lam_data,lam_moveCost,DD_rate,DD_Thresh_quant,sig_logR,q, sel_grow, sel_midpt, lam_Costdist,
          Abunpref_grate, eff_scalar, cv_harv, cv_cpue, cv_effort, Perc_eff_smpld, Prop_sample, Perc_sample_pb,
          FIM_q, num_FIM_samples, w_dist, w_depth, w_DD, w_cost, w_profit),
        file = paste0(save_wd, "/Parameters.txt"), ncolumns=32)
  
  #Writing Dat File
  #This one has clustered composition
  dat_file_Indexnom_compcorr<-list(1, nyear-nyear_init, 0, num_ages-1, M_vec, 
                 Catch_obs[(nyear_init+1):nyear],                             #Observed Catch
                 Index_nom_obs,                                               #Observed CPUE
                 t(Age_Comp_Data[(nyear_init+1):nyear,]),                     #Fishery Age Composition data
                 rowSums(Bio_sample_yr[(nyear_init+1):nyear,]),               #SS for Fishery composition
                 FIM_index_data[(nyear_init+1):nyear],                        #Survey Index
                 t(FIM_biodata[(nyear_init+1):nyear,]),                       #Survey Composition
                 apply(FIM_biosample[(nyear_init+1):nyear,,],1,FUN = sum),    #Sample size of FIM Composition
                 Fec$Fecundity,                                               #Fecundity
                 Wt,                                                          #Wt-at-age
                 c(11,22,33))                                                 #test Vector

  unlink(paste0(save_wd,"/dat_file_Indexnom_compcorr.dat"), recursive = TRUE)
  lapply(dat_file_Indexnom_compcorr, write, paste0(save_wd,"/dat_file_Indexnom_compcorr.dat"), append=TRUE, ncolumns=21)
  
  #This data file has composition data that were randomly sampled
  dat_file_Indexnom_comprnd<-list(1, nyear-nyear_init, 0, num_ages-1, M_vec, 
                Catch_obs[(nyear_init+1):nyear],                             #Observed Catch
                Index_nom_obs,                                               #Observed CPUE
                t(rnd_age_comp[1:(nyear-nyear_init),]),                      #Fishery Age Composition data
                rowSums(Bio_sample_yr[(nyear_init+1):nyear,]),               #SS for Fishery composition
                FIM_index_data[(nyear_init+1):nyear],                        #Survey Index
                t(rnd_surv_comp),                                            #Survey Composition
                apply(FIM_biosample[(nyear_init+1):nyear,,],1,FUN = sum),    #Sample size of FIM Composition
                Fec$Fecundity,                                               #Fecundity
                Wt,                                                          #Wt-at-age
                c(11,22,33))                                                 #test Vector
                 
  unlink(paste0(save_wd,"/dat_file_Indexnom_comprnd.dat"), recursive = TRUE)
  lapply(dat_file_Indexnom_comprnd, write, paste0(save_wd,"/dat_file_Indexnom_comprnd.dat"), append=TRUE, ncolumns=21)
    
  #Saving Ouput from model
  saveRDS(N_wSpace_postM, file=paste0(save_wd, "/N_wSpace_postM.rds" ))
  saveRDS(SSB_space, file=paste0(save_wd, "/SSB_space.rds" ))
  saveRDS(Catch_bio_space, file=paste0(save_wd, "/Catch_bio_space.rds" ))
  saveRDS(Catch_num_space, file=paste0(save_wd, "/Catch_num_space.rds" ))
  saveRDS(Catch_numage_space, file=paste0(save_wd, "/Catch_numage_space.rds" ))
  saveRDS(F_space, file=paste0(save_wd, "/F_space.rds" ))
  saveRDS(Effort_space, file=paste0(save_wd, "/Effort_space.rds" ))
  saveRDS(Sel, file=paste0(save_wd, "/Selectivity.rds" ))
  saveRDS(Sel_FIM, file=paste0(save_wd, "/Sel_FIM.rds" ))
  saveRDS(Effort_midpoints, file=paste0(save_wd, "/Effort_midpoints.rds" ))
  #Data section below
  saveRDS(Catch_obs, file=paste0(save_wd, "/Catch_obs.rds" ))
  saveRDS(Index_nom_true, file=paste0(save_wd, "/Index_nom_true.rds" ))
  saveRDS(Index_nom_obs, file=paste0(save_wd, "/Index_nom_obs.rds" ))
  saveRDS(Cells_sampled, file=paste0(save_wd, "/Cells_sampled.rds" ))
  saveRDS(Bio_samples, file=paste0(save_wd, "/Bio_samples.rds" ))
  saveRDS(Bio_sample_yr, file=paste0(save_wd, "/Bio_sample_yr.rds" ))
  saveRDS(Age_Comp_Data, file=paste0(save_wd, "/Age_Comp_Data.rds" ))
  saveRDS(FIM_sample, file=paste0(save_wd, "/FIM_sample.rds" ))
  saveRDS(FIM_biosample, file=paste0(save_wd, "/FIM_biosample.rds" ))
  saveRDS(FIM_index_data, file=paste0(save_wd, "/FIM_index_data.rds" ))
  saveRDS(FIM_biodata, file=paste0(save_wd, "/FIM_biodata.rds" ))
  }

#Examples for running the functions

#GM Model
#for (i in 1:100){
#Spatial_Model_data(save_wd=paste0("where you want to save SM output",i), extract_wd=paste0("where your OM is located",i),
#              seed=i,PE=TRUE,nyear=90,nyear_init=50,num_cells=1559, num_ports=23, num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0.03, Abunpref_grate=0.00025, eff_scalar=1e5, cv_totaleff=0.25,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.0025, Prop_sample=TRUE, Perc_sample_pb=0.80,
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1)          
#}

#RF Model
#for (i in 1:100){
#Spatial_Model_data(save_wd=paste0("where you want to save SM output_",i), extract_wd=paste0("where your OM is located",i),
#              seed=i,PE=TRUE,nyear=90,nyear_init=50,num_cells=1559, num_ports=23, num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0, Abunpref_grate=0, eff_scalar=1e5, cv_totaleff=0.25,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.0025, Prop_sample=TRUE, Perc_sample_pb=0.80,
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1)          
#}

#Hybrid RF Model
#for (i in 1:100){
#Spatial_Model_data(save_wd=paste0("where you want to save SM output",i), extract_wd=paste0("where your OM is located",i),
#              seed=i,PE="Hybrid_CE",nyear=130,nyear_init=50,num_cells=1559, num_ports=23, num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0, Abunpref_grate=0, eff_scalar=1e5, cv_totaleff=0.25,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.000625, Prop_sample=TRUE, Perc_sample_pb=0.20,
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1)          
#}

#Hybrid GM Model
#for (i in 1:100){
#Spatial_Model_data(save_wd=paste0("where you want to save SM output",i), extract_wd=paste0("where your OM is located",i),
#              seed=i,PE="Hybrid",nyear=130,nyear_init=50,num_cells=1559, num_ports=23, num_ages=21,lam_data=TRUE,lam_moveCost=0.02,DD_rate=0.5,DD_Thresh_quant=0.75,sig_logR=0.3, 
#              q=0.005, sel_grow=2, sel_midpt=2, lam_Costdist=0.03, Abunpref_grate=0.00025, eff_scalar=1e5, cv_totaleff=0.25,
#              cv_harv=0.05, cv_cpue=0.25, cv_effort=0.05, Perc_eff_smpld=0.000625, Prop_sample=TRUE, Perc_sample_pb=0.20,
#              FIM_q=0.001, num_FIM_samples=50, w_dist=1, w_depth=1, w_substrate=1, w_DD=1, w_cost=1, w_profit=1)          
#}
