######################################################################
# PROJECT: SPATIAL SIMULATION OF VIRAL SPREAD WITHIN A 2D LAYER OF CELLS
#
# R-file to create the overview-files used in the simulations
# @Authors: Peter Kumberger & Frederik Graw
######################################################################

# FUNCTION TO CREATE THE OVERVIEW FILE CONTAINING ALL VARIABLES USED FOR SIMULATION
# (originally designed to run several scenarios on a cluster)

create.overviewHCVspread <- function(directory,scen,length.side,sim,timesteps,prob_infectCC,prob_infectDIFF,frac_HCV,prob_dieUI,prob_dieJ,prob_dieI,prob_infinf,max_init_nbr,init_exp_rate,trunc_time,media_change,prolif_time,perc_empty,coupling,r_dieV,runV,ovID,R0,prodRNA,Rcap,rate_exportRNA,rate_degradeRNA,E2_0,degr_e2){

	out <- data.frame(length.side,sim,timesteps,prob_infectCC,prob_infectDIFF,frac_HCV,prob_dieUI,prob_dieJ,prob_dieI,prob_infinf,max_init_nbr,init_exp_rate,trunc_time,media_change,prolif_time,perc_empty,coupling,r_dieV,runV,R0,prodRNA,Rcap,rate_exportRNA,rate_degradeRNA,E2_0,degr_e2)

	names(out) <- c("length.side","sim","time","p_infectCC","p_infectDIFF","frac_HCV","p_dieUI","p_dieJ","p_dieI","p_infinf","max_initnbr","init_exprate","trunc_time","media_change","prolif_time","perc_empty","coupling","r_dieV","run_V","R0","prodRNA","Rcap","rate_exportRNA","rate_degradeRNA","E2_0","degr_e2")
	write.table(out,paste0(directory,"/hcvCluster_overview_",scen,".csv"),sep=",")
	return(out)
}


########################################################################################################
# FIRST SIMULATION RUN - REVISED ( 1TIMESTEP=1min)

create.parms <- function(directory,file.path,pICC,pIDIFF,coupling,pII,degrE2,scen,nbr.runs){
  
  # Parameterization for single run (  # assume 1 timestep = 1 min

  length.side.vec <- 90                         # Dimensions of the viral grid
  time.vec <- 3960								# corresponds to 5 days if 1ts=1min
  sim.vec <- nbr.runs                           # Number of simulation runs
  # assume virion production rate of 100 day^-1 and probability to budd of p=0.25 (can be changed)
  pICC.vec <- pICC
  pIDIFF.vec <- pIDIFF
  frac_HCV.vec <- 1                             # dummy variable from previous simulations
  pdieUI.vec <- 0								# no death rate for uninfected hepatocytes assumed in this simulation
  pdieI.vec <- 0								# no death rate for infected hepatocytes assumed in this simulation
  pdieJ.vec <- 0								# no death rate for infectious hepatocytes assumed in this simulation
  pii.vec <- pII                                # probability of an infected cell of becoming infectious
  max.init.vec <- 25									# maximal nbr of infectious hepatocytes after initiation time
  init.exp.rate.vec <- 0.076/60          # Rate for initiation of infectious hepatocytes (loss of infectivity as estimated from Exp 49C, Uprichard)
  trunc.time.vec <- 17*60            # Time at which initiation of infection is over, either trough administration of anti-E2 or through change of media
  media.change.vec <- 72*60             # Time at which media change removing all extracellular virus in culture (Experimental design)
  prolif.time.vec <-  24*60              # Average time an uninfected hepatocyte needs to proliferate (Timpe Hepatology 2008)
  perc.empty.vec <-   0.6                 # Percentage of empty space (1/(1-perc.empty)=possible fold-increase)
  coupling.vec <- coupling            	# coupling between lattice sites determining how virus can diffuse to neighbouring sites (see manuscript)

  rdieV.vec <- 0.076/60         # Loss of infectivity per extracellular RNA as estimated from Experiment 49C (Uprichard)
  runV.vec <- 1									# how many diffusion runs
  # INTRACELLULAR VIRAL REPLICATION
  R0.vec <- 80.07                   # R.0 as estimated from Experiment 49C (see Durso-Cain et al. 2021)
  prodRNA.vec <- 0.24/60            # Viral reproduction rate per cell as estimated from Experiment 49C (see Durso-Cain et al. 2021)
  rate_exportRNA.vec <- 2.09e-3/60   # Viral export rate per cell as estimated form Experiment 49C (see Durso-Cain et al. 2021)
  #rate_degradeRNA.vec <- 0/60      # No intracellular degradation of viral RNA (assumption from fit to Exp. 49C)
  Rcap.vec <- 7.736e4						# Capacity of positive strand RNA, i.e., R_C, (see Durso-Cain et al. 2021)
  E2_0.vec <- 1e4                  # Initial anti-E2 concentration after media change
  if(scen=="BOTH"){E2_0.vec <- 0}   # If both transmission modes occur, there is no addition of anti-E2
  degr_e2.vec <- degrE2               # Degradation rate of anti-E2 due to viral neutralization
  
  simIFNa <- create.overviewHCVspread(directory,scen,length.side.vec,sim.vec,time.vec,pICC.vec,pIDIFF.vec,frac_HCV.vec,pdieUI.vec,pdieI.vec,pdieJ.vec,pii.vec,max.init.vec,init.exp.rate.vec,trunc.time.vec,media.change.vec,prolif.time.vec,perc.empty.vec,coupling.vec,rdieV.vec,runV.vec,ovID,R0.vec,prodRNA.vec,Rcap.vec,rate_exportRNA.vec,rate_degradeRNA.vec,E2_0.vec,degr_e2.vec)
    return(simIFNa)
}

