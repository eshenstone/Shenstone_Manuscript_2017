rm(list = ls())
#########################################

#This function will create simulated data for a trait with genes that following a geometric series.
# The user will input the number of QTL, the heritabilities to be explored, the additive effect,
# and the output directory

create.simluated.data <- function(){
  #Experimental code
  setwd(home.dir)
  myFRG=GAPIT.Fragment(file.path=NULL,file.from=file.from, file.to=file.to,file.total=NULL,file.G=file.G,
                       file.Ext.G=file.Ext.G,seed=123,SNP.effect="Add",SNP.impute="Middle",
                       genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file.fragment=file.fragment,
                       LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE)
  
  
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
  
  #hm.GD.without.snps.below.0.05.maf <- hm$GD[,-snps.below.0.05.maf]
  
  ###############################
  #temp code#
  hm.GD.without.snps.below.0.05.maf <- hm$GD
  
  genotypes <- data.frame(hm$GI[,1], rep(NA,nrow(hm$GI)),hm$GI[,2:3],rep(NA,nrow(hm$GI)),
                              t(hm.GD.without.snps.below.0.05.maf))
  
  colnames(genotypes) <- c("Snp", "allele", "chr", "pos", "cm", t(as.character(hm$GT)))
  
  #End experimental code
  
  #Create a working directory for the output results:
  dir.create(paste(home.dir,"/", output.dir, sep = ""))
  
  #Set the working directory
  setwd(paste(home.dir,"/", output.dir, sep = ""))
  
  #Randomly select (without replacement) k additive QTN, and assign an effect size
  #seed.number <- sample(-1000000:1000000, 1)
  seed.number <- -875677
  
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.add.QTN <- sample(1:nrow(genotypes), Additive.QTN.number, replace = FALSE)
  Add.QTN.genotypic.information <- genotypes[vector.of.add.QTN,]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Additive.QTN.number,"Add.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Add.QTN.genotypic.information, paste("Genotypic.information.for.", Additive.QTN.number,".Additive.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
  #seed.number <- sample(-1000000:1000000, 1)
  seed.number <- 562991
  
  # Use the seed number to generate uniform random variables
  set.seed(seed.number)
  
  vector.of.epi.QTN <- sample(1:nrow(genotypes), (2*Epistatic.QTN.number), replace = FALSE)
  Epi.QTN.genotypic.information <- genotypes[vector.of.epi.QTN,]
  
  
  #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
  # number used
  write.table(paste("Here_is_the_seed_number:_",seed.number, sep = ""), paste("Seed.number.for.", Epistatic.QTN.number,"Epi.QTN",
                                                                              ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(Epi.QTN.genotypic.information, paste("Genotypic.information.for.", Epistatic.QTN.number,".Epistatic.QTN",
                                                   ".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
  #Create a "base line" trait, which is basically just the additive effects; this is what we would see if 
  # the heritability of the simulated trait were 1
  additive.effect.trait.object <- t(Add.QTN.genotypic.information[,-c(1:5)]) #this was originally the base.line.trait.object
  
  epistatic.effect.trait.object <-t(Epi.QTN.genotypic.information[,-c(1:5)])
  #AEL Changed: - epistatic.effect.trait.object<- epistatic.effect.trait.object[,number.of.epistasis]
  
  colnames(additive.effect.trait.object) <- paste("Chr_",Add.QTN.genotypic.information[,3], "_",Add.QTN.genotypic.information[,4],sep = "")
  colnames(epistatic.effect.trait.object) <- paste("Chr_",Epi.QTN.genotypic.information[,3], "_",Epi.QTN.genotypic.information[,4],sep = "")
  
  #make base.line.trait additive.component and epistatic.component
  additive.component<- as.data.frame(matrix(0, nrow = nrow(additive.effect.trait.object), ncol = 1))
  epistatic.component<- as.data.frame(matrix(0, nrow = nrow(epistatic.effect.trait.object), ncol = 1))
  #base.line.trait <- as.data.frame(matrix(0, nrow = nrow(base.line.trait.object), ncol = 1)) 
  for(i in 1:Additive.QTN.number) additive.component <- additive.component + (additive.effect.trait.object[,i]*(additive.effect^i))
  rownames(additive.component) <- rownames(additive.effect.trait.object)
  colnames(additive.component) <- "Additive.effect"
  additive.genetic.variance <- var(additive.component)
  
  last.number.of.this.loop <- Epistatic.QTN.number - 1
  for(i in 0:last.number.of.this.loop) epistatic.component <- epistatic.component + ((epistatic.effect.trait.object[,((2*i)+1)]*epistatic.effect.trait.object[,((2*i)+2)])*(epistatic.effect^(i+1)))
  rownames(epistatic.component) <- rownames(epistatic.effect.trait.object)
  colnames(epistatic.component) <- "Epistatic.effect"
  epistatic.genetic.variance<- var(epistatic.component)
  
  #Set the row names of the base.line.trait object to the new names
  base.line.trait <- additive.component+epistatic.component
  #base.line.trait.with.new.taxa <- merge(base.line.trait, taxa.name.converter, by.x = "row.names", 
  #                                       by.y = "Old_Taxa_ID")
  
  #the.new.taxa.ids <- as.character(base.line.trait.with.new.taxa[,2])
  #base.line.trait <- as.matrix(base.line.trait.with.new.taxa[,2], nrow = nrow(base.line.trait.with.new.taxa))
  #rownames(base.line.trait) <- as.character(base.line.trait[,3])
  
  #calculate the probability of each genotype lodging
  pi.success <- exp(the.intercept + base.line.trait)/(1+exp(the.intercept + base.line.trait))
  
  #Using this probability, conduct replicates of binomial random variables
  for(i in 1:replicates){
    print(paste("Now simulating trait ", i, sep = ""))
    binomial.rvs <- rep(NA, nrow(pi.success))
    the.seed.number <- NULL
     for(j in 1:nrow(pi.success)){
      seed.number <- sample(-1000000:1000000, 1)
      set.seed(seed.number)
      binomial.rvs[j] <- rbinom(1, num.plants.per.plot, pi.success[j,1])
      the.seed.number <- c(the.seed.number, seed.number)
    }#end for(j in 1:nrow(pi.success)) 
    binomial.rvs <- data.frame(binomial.rvs)
    rownames(binomial.rvs) <- rownames(pi.success)
    
    if(i == 1){
      simulated.data <- binomial.rvs
      seed.numbers <- the.seed.number
    }else{
      simulated.data <- cbind(simulated.data, binomial.rvs)
      seed.numbers <- cbind(seed.numbers, the.seed.number)
      colnames(simulated.data)[i] <- paste(colnames(simulated.data)[i],".",i,sep = "")
      colnames(seed.numbers)[i] <- paste(colnames(seed.numbers)[i],".",i,sep = "")
    }#end if(i == 1)
  }#End for(i in 1:replicates)
  
  #Format the output file for the simulated phenotypes
  simulated.data <- cbind(rownames(base.line.trait),simulated.data)
  colnames(simulated.data)[1] <- "<Trait>"
  colnames(simulated.data)[2] <- paste(colnames(simulated.data)[2],".",1,sep = "")
  
  #Format the output file for the seed numbers
  seed.numbers <- cbind(rownames(base.line.trait),seed.numbers)
  colnames(seed.numbers)[1] <- "Taxa"
  colnames(seed.numbers)[2] <- "the.seed.number.1"
  
  #Add the number of plants per plot to the output file
  simulated.data.a <- data.frame(simulated.data[,1], 
                                 rep(num.plants.per.plot, nrow(simulated.data)),
                                 simulated.data[,2:ncol(simulated.data)])
  simulated.data <- simulated.data.a
  
  
  #Output the m replicates and the seed numbers, formatted for TASSEL
  write.table(seed.numbers, paste("Seed.numbers.for.", replicates,".Reps",
                                            ".Herit.",i,".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  write.table(simulated.data, paste("Simulated.Data.", replicates,".Reps",
                                    ".Herit.",i,".txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
  
  
}#end "create.simluated.data()"





###########################################################################################
###########################################################################################
###########################################################################################
#setwd("/Users/adminuser/Box Sync/Lipka_Mainzer_Chen_Epistasis_Shared_Folder/Simulation_Study")
setwd("C:/Users/shensto2/Desktop/Simulation_Study")
home.dir <- getwd()
#dir.of.GBS.SNPs <- "/Users/adminuser/Desktop/Work/Tocos_NAM_2009_2010/Joint_Linkage_Analysis/GBS_SNPs/"


#Read in the 1,106 markers that are genotyped on the NAM familes scored for kernel color in Chandler et al. 2013

#setwd(dir.of.GBS.SNPs)
#genotypes <- read.table("4K_SNPsmdp_genotype_test1.hmp.txt", head = TRUE)
#setwd(home.dir)
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
#library(rrBLUP)


###############
#User input

#Number of additive QTN (k)
Additive.QTN.number <- 1


#Number of epistatic QTN (m)
Epistatic.QTN.number <- 1

#Vector of heritabilities to investigate
#heritabilities.vector <- 0

#Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))
additive.effect <- .9


#Size of the epistatic effect of the largest QTL (must be (-1,1) but preferably (0,1)|
epistatic.effect <- 0

#Number of replicates of simulated phenotypes for each heritability (m)
replicates <- 100


file.G="4K_SNPsmdp_genotype_test" 
file.Ext.G = "hmp.txt"


file.from = 1
file.to = 1
file.fragment = 100000 #Note: if your data set has more that 100,000 SNPs, please change
                       # this number.


#Output directory
output.dir <- paste(Additive.QTN.number,"_Add_QTN",Epistatic.QTN.number,"_Epi_QTN_h.2_",
                    "_add.eff_", additive.effect,"_epis.eff_", epistatic.effect,"_reps_", replicates,"_sim_6_15", sep = "")

the.intercept <- 0 #Es, we should talk about the intercept
num.plants.per.plot <- 15
################
#Create the simulated data
create.simluated.data()









