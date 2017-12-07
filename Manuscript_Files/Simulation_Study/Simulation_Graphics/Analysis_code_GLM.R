############# Parameters ###################
simulated.QTL.info.dir = "C://Users//shensto2//Desktop//Simulation_Study//old//Simulations_Model2//Setting_6";
summary.dir = "C://Users//shensto2//Desktop//Simulation_Study//old//Simulations_Model2//Setting_6";
simulation.results.dir = "C://Users//shensto2//Desktop//Simulation_Study//old//Simulations_Model2//Setting_6";
GLM_results = "setting6e_GLM";
file.name.prefix = "Top.100.SNPs.for.All.Traits";
number.of.simulated.additive.QTN = 1;
number.of.simulated.epistatic.QTN.pairs = 1;
setting.number = 6
include.additive.QTL = TRUE;
include.epistatic.QTL = FALSE;

############################################# 
  setwd("C:/Users/shensto2/Desktop/Simulation_Study/Simulation_Graphic/")
  home.dir = getwd() 
  ##Set up some global variables
  flanking.bp.dist.from.QTN = 250000
  number.of.traits.per.setting = 100
  chm.to.analyze = c(1,2,3,4,5,6,7,8,9,10)
  number.of.simulations = 100
  size.of.maize.chromosomes = c(301331005, 237007151, 232096209, 241402927, 217748528,
                                 169132432, 176664055, 175775017, 156591811, 150168908) #This will make the Manhattan plots look good

################Create FDR corrected results ##################### 
#source("/Users/chen398/Documents/Epistasis Simulation/New_simulated_settings/Setting_4_1_Add_QTN4_Epi_QTN_h.2_0.5_add.eff_0_epis.eff_0.6_reps_100/GLM/Create_FDR_corrected_GLM.R")

################Create statistics and plots ##################### 
  source("R_for_Epistasis_GLM_Analysis_no_parameters_AEL.R")
  time.start <- Sys.time()
  assess.the.results()
  time.stop <- Sys.time()
  
  time.stop-time.start
