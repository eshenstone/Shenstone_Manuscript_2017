# GAPIT - Genome Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# Written by Zhiwu Zhang, Alex Lipka and Feng Tian 
# Last update: August 10, 2011 

#Step 0: Set directory to load GAPIT source files (this step is omited for using package)
#######################################################################################

setwd("C:/Users/shensto2/Desktop/Binomial_Data/55K_GAPIT")

library('MASS') # required for ginv
#source("https://bioconductor.org/biocLite.R")
#biocLite("multtest")
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
#source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")


#Step 1: Set data directory and import files
#######################################################################################

mydataPath="/Users/shensto2/Desktop/Binomial_Data/55K_GAPIT/"

#Phenotypic Data
myY  <- read.table(paste(mydataPath,"BLUPs_2016_Inoculated.txt",sep=""), head = TRUE)
head(myY)
dim(myY)


myG <- read.delim(paste(mydataPath,"SNP55K_maize282_AGPv2_20100513_1.hmp.txt",sep=""), head = FALSE)
head(myG[,1:7])

myKI <- read.csv(paste(mydataPath,"Comprehensive_K.csv",sep=""), head = FALSE)
head(myKI[,1:7])

myCV <- read.csv(paste(mydataPath,"GAPIT.PCA.csv",sep=""), head = TRUE)
#myCV <- myCV[,c(1,ncol(myCV))]
head(myCV)


#your genotype file needs to be in hapmap format
###################################################################################
#Step 2:  Run GAPIT on the first trait (In this sitation I had one trait with 97 observations,
#and then six traits with 150 obersvations which is why I do this in two steps)
#######################################################################################
#Create a new phenotype file that only has the taxa names and trait 1 


myGAPIT <- GAPIT(
Y=myY,				#This is phenotype data
G=myG,				#This is genotype data,set it to NULL with multiple genotype files
#GD=myGD,
#GM=myGM,

file.path=mydataPath,		#The location of genotype files


group.from=262,			#my file had 97 obersvations, since I don't want to do any sort of compression i set both groups to 97
group.to=262,			#Upper bound for number of group
group.by=1,				#range between 1 and number of individuals, smaller the finner optimization 
#kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
#kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)

kinship.algorithm = "Loiselle",

SNP.impute = "Major",
SNP.MAF = 0.05,
cutOff = 0.00,
KI=myKI,  			#This is kinship data, set it to NULL in case that geneotype files are used for estimation
CV=myCV,				#This is the covariate variables of fixed effects, such as population structure

Model.selection = TRUE


)
#





