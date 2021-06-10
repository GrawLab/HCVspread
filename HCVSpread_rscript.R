######################################################################
# PROJECT: SPATIAL SIMULATION OF VIRAL SPREAD WITHIN A 2D LAYER OF CELLS
#
# R-file used to define the global variables and to start the simulations
# @Authors: Peter Kumberger & Frederik Graw
######################################################################

# Function to create a shared library with the actual parameter combination which should be tested so that it
# can be called by R
##################################################
# (only needed for HCV in vitro experiment)
# extra function to have different starting points
loss.inf <- function(l,a,t){
  inf <- function(x){l*exp(-a*x)}
  integrate(inf,lower=0,upper=t)$value
}

# (a) Truncated cumulative probability to get infected
trunc.infLI <- function(l,a,interval,t){
  cum.inf <- function(x){loss.inf(l=l,a=a,t=x)}
  (cum.inf(t)-cum.inf(interval[1]))/(cum.inf(interval[2])-cum.inf(interval[1]))
}

# (b) Random sample for time to infection (n=number of samples)
rSample.cdfLI <- function(l,a,interval,n){
  out <- NULL
  for (r in 1:n){
    u <- runif(1,min=0,max=1)
    cif <- function(t){trunc.infLI(l=l,a=a,interval,t=t)-u}
    out <- c(out,uniroot(cif,interval)$root)
  }
  out
}


##################################################
# directory : path to the folder containing the files used for simulation
# scen : Name of the simulation scenario used for simulation
Sim.prohep <- function(directory,scen){
  ptm <- proc.time()
  
  # read the simulation overview file containing all the different parameter combinations
  inter <- read.table(paste0(directory,"/HCVSpread_overview_",scen,".csv"),sep=",",header=TRUE)
  # choose the corresponding simulation
  simu <- inter
  print(simu)
  
  # define the matrix for the global parametrs
  globvar <- matrix(data=NA,nrow=8,ncol=1)
  
  # write the filename(s) into the parameter file containing the global parameters (HEPstat,HEPinf,HEPinfectious,HEPdeath)
  globvar[1,1] <- paste0("ofstream outHEPstat(\"",paste0(directory,"/HEPstat_",scen,".csv"),"\", ios::out);")
  globvar[2,1] <- paste0("ofstream outHEPinf(\"",paste0(directory,"/HEPinf_",scen,".csv"),"\", ios::out);")
  globvar[3,1] <- paste0("ofstream outHEPinfectious(\"",paste0(directory,"/HEPinfectious_",scen,".csv"),"\", ios::out);")
  globvar[4,1] <- paste0("ofstream outInfectClusterType(\"",paste0(directory,"/HEPInfectClusterType_",scen,".csv"),"\", ios::out);")
  globvar[5,1] <- paste0("ofstream outEmpty(\"",paste0(directory,"/Empty_",scen,".csv"),"\"",", ios::out);")
   globvar[6,1]<- paste0("ofstream outTEST(\"",paste0(directory,"/TEST_",scen,".csv"),"\", ios::out);")
  # determine the size of the simulated tissue layer
  globvar[7,1] <- paste0("#define LENGTH_SIDE ",simu[,"length.side"])
  
  # write the parameter file
  write.table(globvar,paste0(directory,"/HCVSpread_parameter.h"),row.names=FALSE,col.names=FALSE,quote=FALSE)
  
  #####################
  # # Load gnu-compiler
  # create the shared library, the old one has to be deleted first
  system(paste0("rm ",directory,"/HCVSpread_mainfile.dll"))
  system(paste0("rm ",directory,"/HCVSpread_mainfile.o"))

  system(paste0("R CMD SHLIB ",directory,"/HCVSpread_mainfile.cc -o ",directory,"/HCVSpread_mainfile.dll"))

  
  # LOAD THE COMPILED FUNCTION
  dyn.load(paste0(directory,"/HCVSpread_mainfile.dll"))
  

  #####################
  # BASIC SIMULATION FUNCTION
  SimBase <- function(sim,timesteps,dummy,dummyDIST,prob_infectCC,prob_infectDIFF,frac_HCV,prob_dieUI,prob_dieJ,prob_dieI,prob_infinf,max_init_nbr,init_exp_rate,trunc_time,media_change,prolif_time,perc_empty,coupling,rate_dieV,runsV,R0,prodRNA,Rcap,rate_exportRNA,rate_degradeRNA,E2_0,degr_e2,initINF)
    
    .C("Simulate_infection",as.integer(sim),as.integer(timesteps),as.double(dummy),as.double(dummyDIST),as.double(prob_infectCC),
       as.double(prob_infectDIFF),as.double(frac_HCV),as.double(prob_dieUI),as.double(prob_dieJ),as.double(prob_dieI),
       as.double(prob_infinf),as.integer(max_init_nbr),as.double(init_exp_rate),as.integer(trunc_time),as.integer(media_change),as.integer(prolif_time),as.double(perc_empty),as.double(coupling),as.double(rate_dieV),as.double(runsV),as.double(R0),as.double(prodRNA),as.double(Rcap),as.double(rate_exportRNA),as.double(rate_degradeRNA),as.double(E2_0),as.double(degr_e2),as.integer(initINF))
  
  SIMULATION <- function(sim,timesteps,dummy,dummyDIST,prob_infectCC,prob_infectDIFF,frac_HCV,prob_dieUI,prob_dieJ,prob_dieI,prob_infinf,max_init_nbr,init_exp_rate,trunc_time,media_change,prolif_time,perc_empty,coupling,rate_dieV,runsV,R0,prodRNA,Rcap,rate_exportRNA,rate_degradeRNA,E2_0,degr_e2,initINF){
    
    set.seed(unclass(Sys.time()))
    SimBase(sim=sim,timesteps=timesteps,dummy=dummy,dummyDIST=dummyDIST,prob_infectCC=prob_infectCC,prob_infectDIFF=prob_infectDIFF,frac_HCV=frac_HCV,prob_dieUI=prob_dieUI,prob_dieJ=prob_dieJ,prob_dieI=prob_dieI,prob_infinf=prob_infinf,max_init_nbr=max_init_nbr,init_exp_rate=init_exp_rate,trunc_time=trunc_time,media_change=media_change,prolif_time=prolif_time,perc_empty=perc_empty,coupling=coupling,rate_dieV=rate_dieV,runsV=runsV,R0=R0,prodRNA=prodRNA,Rcap=Rcap,rate_exportRNA=rate_exportRNA,rate_degradeRNA=rate_degradeRNA,E2_0=E2_0,degr_e2=degr_e2,initINF=initINF)
  }
  
  # definition of the parameters for the actual simulation
  timesteps <- as.integer(simu["time"])							# number of timesteps considered in each simulation
  sim <- as.integer(simu["sim"])                      # Number of simulations
  prob_infectCC <- as.double(simu["p_infectCC"])					# probability of a virion to infect by cell-to-cell transmission
  prob_infectDIFF <- as.double(simu["p_infectDIFF"])				# probability of a virion to infect as a free diffusing viral particle
  frac_HCV <- as.double(simu["frac_HCV"])							# fraction of infectious material per viral RNA in hepatocyte
  prob_dieUI <- as.double(simu["p_dieUI"])						# probability of an uninfected cell to die
  prob_dieJ <- as.double(simu["p_dieJ"])							# probability of an infected cell to die
  prob_dieI <- as.double(simu["p_dieI"])							# probability of an infectious cell to die
  prob_infinf <- as.double(simu["p_infinf"])						# probability of an infected cell to become infectious
  max_init_nbr <- as.integer(simu["max_initnbr"])			# maximal nbr of infectious hepatocytes after initiation time
  init_exp_rate <- as.double(simu["init_exprate"])    # Rate for initiation of infectious hepatocytes
  trunc_time <- as.integer(simu["trunc_time"])       # Time at which initiation of infection is over
  media_change <- as.integer(simu["media_change"])    # Time of secound media change (after trunc_time)
  prolif_time <- as.integer(simu["prolif_time"])      # Average time for uninfected hepatocyte to proliferate
  perc_empty <- as.double(simu["perc_empty"])         # Percentage of empty space at start of experiment
  coupling <- as.double(simu["coupling"])							# rate of coupling between adjacent sites in the viral grid
  rate_dieV <- as.double(simu["r_dieV"])							# death rate of free diffusing virus
  runsV <- as.double(simu["run_V"])								# how often does the virus diffuse?
  R0 <- as.double(simu["R0"])                         # Initial positive-strand RNA in a cell after becoming infectious
  prodRNA <- as.double(simu["prodRNA"])						    # Production rate of intracellular positive strand HCV RNA
  Rcap <- as.double(simu["Rcap"])									    # Capacity of intracellular positive strand HCV RNA
  rate_exportRNA <- as.double(simu["rate_exportRNA"])	# Budding rate of intracellular positive strand HCV RNA
  rate_degradeRNA <- as.double(simu["rate_degradeRNA"])	# Degradation rate of intracellular positive strand HCV RNA
  E2_0 <- as.double(simu["E2_0"])                       # Initial anti-E2 concentration after media change
  degr_e2 <- as.double(simu["degr_e2"])                 # Degradation of anti-E2 due to extracellular virus neutralization
  
  t0.vec <- 0
  
  # actual call of the simulation function
  for (s in c(1:sim)){
    time_start <- as.integer(t0.vec)
    SIMULATION(s,timesteps,0,0,prob_infectCC,prob_infectDIFF,frac_HCV,prob_dieUI,prob_dieJ,prob_dieI,prob_infinf,max_init_nbr,init_exp_rate,trunc_time,media_change,prolif_time,perc_empty,coupling,rate_dieV,runsV,R0,prodRNA,Rcap,rate_exportRNA,rate_degradeRNA,E2_0,degr_e2,time_start)
  }
  proc.time() - ptm
}

