/*Read in the data from the field*/
data lodging;
length plot$20;  
infile 'C:\Users\shensto2\Desktop\Binomial_Data\Model 3\Input_file_for_Model_3_Pheno_2795_SNPs.txt' delimiter='09'x MISSOVER DSD lrecl=35000 firstobs=2;
input plot$ block score n Y PC1 PC2 PC3 SNP1-SNP2795; /*? about column labels */
run;

/*Read in the kinship matrix*/
data Kin;
        filename in "C:\Users\shensto2\Desktop\Binomial_Data\Model 3\Comprehensive_K _col_parm_new.txt";
/*		length Parm$15;*/
        infile in  dlm='09'x linesize=9000 firstobs=2;
        input Parm Row Col1-Col271;
run;

data Results;
run;

/*Begin macro - Esperanza, please highlight and run this macro. This will "read this macro"
into SAS*/
%macro logistic(num);
/*Do loop*/
%do i = 1 %to &num;
	%put 'Working on iteration' &i;
	/*Fit the model*/
	proc glimmix data=lodging;
	class plot block ;
	model Y/ n= block PC1 PC2 PC3 SNP&i / link=logit dist=binomial solution;
	random plot /type=lin(1) lDATA=Kin solution;
	ods output Tests3= Results_SNP&i;
	run;

	data Results_SNP&i; set Results_SNP&i;
	if Effect = "block" then delete;
	if Effect = "PC1" then delete; 
	if Effect = "PC2" then delete; 
	if Effect = "PC3" then delete;  
	run; 

	/*Obtain the output you need*/
	data Results; set Results Results_SNP&i;
	run;
%end;
/*End do loop*/
%mend;
/*End macro*/

/*Esperanza, once you have ran this macro, all you need to do is change the number
in parantheses after "%logistic(THIS IS THE NUMBER I AM TALKING ABOUT)" so that
it equals the number of SNPs you want to test*/
%logistic(2795)

/*The next line of code is a housekeeping step.
If you take a look at the first row of "Results", 
you will see that that it is a blank row. So 
the next line of code will get rid of that blank row*/
data Results; set Results;
if NumDF = "." then delete;
run;

/*Final step is to export the results*/
proc export data=Results
	
   outfile="C:\Users\shensto2\Desktop\Binomial_Data\Model 3\Test_Results.csv"
 	dbms=csv replace; 
run;
