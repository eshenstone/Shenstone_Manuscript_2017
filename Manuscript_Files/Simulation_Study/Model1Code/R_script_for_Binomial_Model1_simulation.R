

#########Obtain the list of paths to the simulated traits
#set the working directory
setwd("C:/Users/shensto2/Desktop/Simulation_Study/Simulations_Model1")
home.dir <- getwd()

#Write down the path to the directory where you want to put your results
results.dir <- "C:/Users/shensto2/Desktop/Simulation_Study/Simulations_Model1/Model1_Simulation_Results_09022017"


#Set count == 1
the.full.path.to.folders <- list.dirs(home.dir)



#Set the working directory

setwd("C:/Users/shensto2/Desktop/Binomial_Data")
home.dir <- getwd()

#Read in the GAPIT code and prerequisite files
library(rrBLUP)
library('MASS')
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/previous/gapit_functions20160408.txt")
#source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.R")
source("http://zzlab.net/GAPIT/emma.txt")
source("Shenstone_GAPIT_FDR_for_LR_Code.R")

#Input parameters for the genotype files
file.G="SNP55K_maize282_AGPv2_20100513_"
file.Ext.G="hmp.txt"
file.from=1
file.to = 1
SNP.fraction=1
group.from = 252 
group.to = 252
SNP.fraction = 0.1
file.fragment = 100000
PCA.total=3
#Read in the phenotypic data


#Reset the working directory to home directory
setwd(home.dir)


#Read in the genotypic data, use GAPIT to convert it into numeric format

myFRG=GAPIT.Fragment(file.path=NULL,file.from=file.G, file.to=file.to,file.total=NULL,file.G=file.G,
                     file.Ext.G=file.Ext.G,seed=123,SNP.effect="Add",SNP.impute="Middle",
                     genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file.fragment=file.fragment,
                     LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)


#the above code is still only reading file 11 (my fake one)
#I would like to read all the hapmap files in and then save it as an image so that during
#the K folds we dont need to reload them for each fold- it took 4 days to do 2 folds......
#At that rate I won't finish my masters before 2020
#GD=myFRG
#Old code
#myFRG=GAPIT.Fragment(file.path=NULL,file.from=file.G, file.to=file.to,file.total=NULL,file.G=file.G,
#                     file.Ext.G=file.Ext.G,seed=123,SNP.fraction=SNP.fraction,SNP.effect="Add",SNP.impute="Middle",
#                     genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file.fragment=file.fragment,
#                     file=1,frag=1,LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)


#End old code

hm=GAPIT.HapMap(G = myFRG$G,SNP.effect="Add",SNP.impute="Major")

#####################################
#Obtain the mafs of all SNPs

#Total number of lines
ns <- nrow(hm$GD)

#Sum of the allele scores for each SNP
ss <- apply(hm$GD, 2, sum)

#Combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
maf.matrix <- rbind((.5*ss/ns), (1-(0.5*ss/ns)))

#Copy the minor allele frequencies for all SNPs
maf <- apply(maf.matrix, 2, min)

#Find out which SNPs have MAF < 0.05
snps.below.0.05.maf <- which(maf < 0.05)

# Remove these SNPs from hm$GD

hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf]
the.SNP.info <- hm$GI[-snps.below.0.05.maf,]


###############################

#something I stole from rrblup source code
#not sure what CV is but i hope this solve my error from #
#line : qc=GAPIT.QC(Y = myY, GT = hm$GT, CV = CV, GK = GK)
#which gives error: Error in GAPIT.RemoveDuplicate(CV) : object 'CV' not found



GK <- cbind(hm$GT, hm.GD.without.snps.below.0.05.maf)


#for loop for each of the simulated traits
for(j in the.full.path.to.folders){
  
  #Go into a given folder, extract the simulation number
  setwd(j)
  
  this.simulation.number <-  substr(unlist(strsplit(j, "_"))[15], 4,1000)
  
  #try to read in 
  all.Y <- try(read.table("Simulated.Data.100.Reps.no.trials.100.txt", head = TRUE))
  
  #If the folder does not contain any simulated output data, go to the next directory
  if(inherits(all.Y, "try-error")) next;

 
  #count <- 0
  
  #for loop through all of the traits
  for(k in 3:ncol(all.Y)){
    this.SNP.name <- NULL
    this.chromosome <- NULL
    this.position <- NULL
    this.LR.P.value <- NULL
    this.slope.estimate <- NULL
    Y <- data.frame(all.Y[,c(1,2,k)])
    for(i in 1:nrow(the.SNP.info)){#For loop through all markers
    
      if(i %% 100 == 0)print(paste("Simulation: ",  this.simulation.number, "Trait: ", (k-2),
                                   "Number of Markers: ", i, sep = ))
      #Extract the i'th SNP
      this.SNP <- data.frame(GK[,1], as.numeric(GK[,(i+1)]))
    
    
      #Merge the i'th SNP to the phenotypes
      pheno.and.this.SNP <- merge(Y, this.SNP, by.x = "the.data...1.", 
                                  by.y = "GK...1.", all.x = TRUE)
    
      #Give the SNP column a name that will consistent throughout all tested SNPs
      colnames(pheno.and.this.SNP)[length(colnames(pheno.and.this.SNP))] = "this.SNP"
      colnames(pheno.and.this.SNP)[which(colnames(pheno.and.this.SNP) == paste("binomial.rvs.",(k-2),sep =""))] = "number.of.successes"
     
      #Remove all rows that have at least one missing value
      pheno.and.this.SNP <- pheno.and.this.SNP[complete.cases(pheno.and.this.SNP),]
      
      #Calculate the number of failures
      number.of.failures <- pheno.and.this.SNP$total.no.trials - pheno.and.this.SNP$number.of.successes
      
      pheno.and.this.SNP.a <- data.frame(pheno.and.this.SNP[,c(1,2)], number.of.failures, 
                                         pheno.and.this.SNP[,3:ncol(pheno.and.this.SNP)])
      
      pheno.and.this.SNP <- pheno.and.this.SNP.a
      
      #Fit a null model
      nullmod=glm(cbind(number.of.successes, number.of.failures)~1,family=binomial(logit),data=pheno.and.this.SNP)	
      
      
      #Fit the logistic regression model
      glm.out = glm(cbind(number.of.successes, number.of.failures)~  this.SNP, family=binomial(logit), data=pheno.and.this.SNP)
      
      this.output <- summary(glm.out)
    
      the.slope.estimate <- this.output$coefficients[nrow(this.output$coefficients)]
     
      #Obtain the likelihood ratio test statistic (and corresponding P-value, of course)
      LR.test <- anova(nullmod,glm.out,test="Chisq")
      
      the.LR.P.value <- LR.test$`Pr(>Chi)`[2]
    
      
      #Collate Chromsome, position, (effect estimate), LRTS and P-value into an output file
      the.SNP.name <- as.character(the.SNP.info[i,1]) 
      the.chromosome <- as.character(the.SNP.info[i,2])
      the.position <- as.character(the.SNP.info[i,3])
      
      #Append all of the results for this SNP into an output file
     # if(count == 0){
    #    the.results <- data.frame(c(the.SNP.name, the.chromosome, the.position,  the.LR.P.value, the.slope.estimate))
    #  }else{
    #    the.results <- rbind(the.results, c(the.SNP.name, the.chromosome, the.position, the.LR.P.value, the.slope.estimate))
    #  }#end if(count == 0)
      
      this.SNP.name <- c(this.SNP.name, the.SNP.name)
      this.chromosome <- c(this.chromosome, the.chromosome)
      this.position <- c(this.position, the.position)
      this.LR.P.value <- c(this.LR.P.value, the.LR.P.value)
      this.slope.estimate <- c(this.slope.estimate, the.slope.estimate)
      
      
      #count <- count + 1
      
    }#End for(i in 1:nrow(the.SNP.info))
    
    the.results <- data.frame(this.SNP.name, as.numeric(this.chromosome), 
                              as.numeric(this.position), this.LR.P.value, this.slope.estimate)
    
    #Sort the results from smallest to largest P-value
    the.results <- the.results[order(the.results[,4]),]
    
    #Label the columns
    colnames(the.results) <- c("SNP.ID", "Chr", "Pos",  "P.value", "Slope.Estimate")
    
    
    #To-do list
    # shoot holes into what we just did (e.g., look over the.results, the input files, the models 
    #  that were fitted), see if we overlooked anything
    #
    # 2.) Incoporate the Benjamini-Hochberg procedure to obtain FDR-adjusted P-values
    
    the.results.with.FDR <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI=the.results, 
                                                       FDR.Rate = 0.05, FDR.Procedure = "BH")
    
    #Change the working directory to the 55K logistic regression results
    setwd(results.dir)
  
    
    write.csv(the.results.with.FDR, paste("Model_1_Results_Setting",this.simulation.number, "_Trait_", (k-2),".csv", sep = ""), 
              row.names = FALSE)
    
    #
    # 3.) Create a Manhattan plot of the results
    the.results.for.manhattan.plot <- the.results.with.FDR$PWIP[,-1]
    
    #Get rid of any P-values that are "NA"
    the.results.for.manhattan.plot <- the.results.for.manhattan.plot[complete.cases(the.results.for.manhattan.plot),]
    
    
    GAPIT.Manhattan(GI.MP = the.results.for.manhattan.plot, name.of.trait =  
                    paste("Model_1_Results_Setting_",this.simulation.number,"_Trait_", (k-2), sep = ""), 
                    plot.type = "Genomewise",DPP=50000,cutOff=0.01)
    
  }# end for(k in 3:ncol(all.Y))  
   
}#end for(j in the.full.path.to.folders)              




