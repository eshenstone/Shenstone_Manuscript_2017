#source("https://bioconductor.org/biocLite.R")
#biocLite("multtest")
library(multtest)

### Post-process the GLM results 
setwd("C:/Users/shensto2/Desktop/Simulation_Study/old/Simulations_Model2/Setting_6//")
home.dir <- getwd()
number.of.simulated.traits.per.setting <- 100
file.name.prefix <- "GAPIT.binomial.rvs."
file.name.suffix<-".GWAS.Results"
#this.FDR.rate <- 0.5
top.number.of.SNPs.to.extract <- 100
simulation.results.dir <- "//Simulations_Model2//Setting_6//"
#Set the working directory to the setting where all of the simulation results are kept
setwd(paste(home.dir, simulation.results.dir, sep = ""))

count.traits <- 0
#Read in all of the results, and put them all in one big file


for(i in 1:number.of.simulated.traits.per.setting){
  print(paste("Currently working on trait ", i, sep = ""))
  these.trait.results <- read.csv(paste(file.name.prefix,i,file.name.suffix,".csv", sep = ""), head = TRUE)
  #these.trait.results <- these.trait.results[which(these.trait.results$PWIP.FDR_Adjusted_P.values <= this.FDR.rate),]
  these.trait.results <- these.trait.results[order(these.trait.results$FDR_Adjusted_P.values),]
  these.trait.results <- these.trait.results[1:top.number.of.SNPs.to.extract,]
  #remove all SNPs that have identical chromosomal/bp positions
  these.trait.results <- these.trait.results[order(these.trait.results[,3]),]
  these.trait.results <- these.trait.results[order(these.trait.results[,2]),]
  vector.of.identical.phys.pos <- rep(FALSE, nrow(these.trait.results))
  for(j in 2:nrow(these.trait.results)){
    if((these.trait.results[j,2] == these.trait.results[(j-1),2])&
       (these.trait.results[j,3] == these.trait.results[(j-1),3])){
      vector.of.identical.phys.pos[j] = TRUE
    }   
  }#end for(j in 1:nrow(these.trait.results))
  these.trait.results <- these.trait.results[which(vector.of.identical.phys.pos==FALSE),]
  
  these.trait.results <- these.trait.results[order(these.trait.results[,4]),]
  
  #If nothing is significant at the given false discovery rate, then move on to the iteration of the loop
  if(nrow(these.trait.results)==0){
    print(paste("--------------No signficiant P-values for trait ",i, "!!!!!---------------------", sep = ""))
    next;
  }
  these.trait.results <- data.frame(c(rep(paste("Trait",i, sep = "")), these.trait.results))
  colnames(these.trait.results)[1] <- "Trait.Number"
  these.trait.results[,1] <- as.character(these.trait.results[,1])
  these.trait.results[,2] <- as.character(these.trait.results[,2])
  if(count.traits == 0){
    all.trait.results <- these.trait.results
  }else{
    all.trait.results <- rbind(all.trait.results, these.trait.results)
  }#end if(count.traits == 0) 

  count.traits <- count.traits + 1
}  #end for(i in 1:number.of.simulated.traits.per.setting)


#write.csv(all.trait.results,paste("All.Sig.SNPs.At.",this.FDR.rate,".FDR.csv",sep = ""),
#          eol = "\n", na = "NA", row.names = FALSE)
write.csv(all.trait.results,paste("Top.",top.number.of.SNPs.to.extract, ".SNPs.for.All.Traits.csv",sep = ""),
          eol = "\n", na = "NA", row.names = FALSE)


#####################################################
######### Don't run this, this is some code code that I "ripped off" to write this script
PWI<- read.table(paste(GLM_results,".txt",sep=""), header = TRUE) 

#Parse results by traits (100 traits totally; 100*1106=110600 phenotype values)
PWI.list <- lapply(1:100, function(i){
  return(PWI[((i-1)*1106+1):((i-1)*1106+1106),])
})

###Create a modified version of the BH_FDR function 
#####################################################
"GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure"<-
  function(PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH"){
    #Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
    #Output: PWIP, number.of.significant.SNPs
    #Authors: Alex Lipka and Zhiwu Zhang 
    # Last update: May 5, 2011 
    ##############################################################################################
    #Make sure that your compouter has the latest version of Bioconductor (the "Biobase" package) and multtest
    
    # if(is.null(PWI))
    # {
    #   PWIP=NULL
    #   number.of.significant.SNPs = 0
    # }
    
    # if(!is.null(PWI))
    # {  
    
    #library(multtest)
    
    #  if(dim(PWI)[1] == 1){
    #    PWIP <- cbind(PWI, PWI[4])
    #    colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
    #  }
    
    #if(dim(PWI)[1] > 1){  # number of row 
    #mt.rawp2adjp Performs the Simes procedure.  The output should be two columns, Left column: originial p-value
    #Right column: Simes corrected p-value
    
    res <- mt.rawp2adjp(PWI[,6], proc=FDR.Procedure, alpha=FDR.Rate, na.rm = FALSE)
    
    #This command should order the p-values in the order of the SNPs in the data set
    adjp <- res$adjp[order(res$index), ]
    
    #round(adjp[1:7,],4)
    #Logical statment: 0, if Ho is not rejected; 1, if  Ho is rejected, by the Simes corrected p-value
    #  temp <- mt.reject(adjp[,2], FDR.Rate)
    
    #Lists all number of SNPs that were rejected by the BY procedure
    #temp$r
    
    #Attach the FDR adjusted p-values to AS_Results
    
    PWIP <- cbind(PWI, adjp[,2])
    
    #Sort these data by lowest to highest FDR adjusted p-value
    PWIP <- PWIP[order(PWIP[,14]),]
    
    colnames(PWIP)[14] <- "FDR_Adjusted_P-values"
    #  number.of.significant.SNPs = temp$r
    #}
    #print("GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure accomplished successfully!")
    # }  
    #return(list(PWIP=PWIP, number.of.significant.SNPs = number.of.significant.SNPs))
    return(PWIP)
  }#GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure ends here

###Execute the BH-FDR function for GLM 
#####################################################

##Markers are adjusted within each trait (100 traits totally)
##PWIP[[1]] stores FDR.adjusted p values that are smaller than 0.05 from the first trait
PWIP<- lapply(1:100, function(i){
  tmp_PWIP<-GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.list[[i]], FDR.Rate = 0.05, FDR.Procedure = "BH")
  return(tmp_PWIP[which(tmp_PWIP[,14]<0.05),])
})

# PWIP<- lapply(1:100, function(i){
#   tmp_PWIP<-GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.list[[i]], FDR.Rate = 0.05, FDR.Procedure = "BH")
#   return(tmp_PWIP)
# })
###Post-processing GLM results 
#####################################################

## This concatenate each data frame. Using Reduce() to "rbind" is faster than using a for loop!
GLM <- Reduce(rbind, PWIP) 

## To check whether all data frames are merged together correctly
#nrow=c()
#for (i in 1:100){
#  nrow[i] <- nrow(PWIP[[i]])
#}
#sum(nrow)

## Write a table for the GLM_results  
#write.table(GLM, file=paste("FDR_adjusted_", GLM_results, ".txt", sep = ""), sep="", eol="\n", na="NA", row.names = TRUE, col.names = TRUE) 

write.csv(GLM, file=paste("FDR_adjusted_", GLM_results, ".csv", sep = ""),
          eol = "\n", na = "NA", row.names = TRUE)

## Note: the generated csv file contains 1 extra column [,1]

# ### This function is from Code_for_Assessing_Epi_Sim_Study_20160512.R 
# #######