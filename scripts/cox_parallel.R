#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(tidyr)
library(survival)
library(parallel)
library(compiler)
library(data.table)

ds_folder= args[1] #'/mnt/work/hunt/dosage/'  # dosage folder
phenofile=  args[2] #'/mnt/work/hunt/pheno/HUNT_PROM_surv_moms' # path to phenotype file
IDfile= args[3] #'./raw_data/samples_ID' # file with sample IDs with same order as in dosage files
ID= args[4] #'MOR_PID' # ID name
outfile= args[5] #'/mnt/work/pol/gwas/res/survival/momsHUNT_PROM_chr' # path to output file
time_t= 'SVLEN_DG' # time variable name
outcome= 'PROM' # outcome variable name
covars= c('FAAR', 'PC1','PC2','PC3', 'PC4', 'PC5', 'PC6','PARITY0') # covariate variable name (multiple are accepted)

options(stringsAsFactors=FALSE)

flist= list.files(path= ds_folder, pattern='.gz')

cox_funk= function(snp){
	X= cbind(snp, covars_m)
        cox_coef= try(unlist(summary(coxph(Surv( time_vec, outcome_vec)~ X, na.action = na.omit))[c(4,5,7)])[c(1,3,4,3+ (length(covars)+1)*2 + 1, 3 + (length(covars)+1)*4 +1)])
	if(class(cox_coef) == "try-error") {
    cox_coef= c(NA,NA,NA,NA)
    }
        return(cox_coef)
}

cox_zph_funk= function(snp){
	X= cbind(snp, covars_m)
        cox_zph= cox.zph(coxph( Surv( time_vec, outcome_vec)~ X, na.action = na.omit ) )$table[1,c(1,3)]
	return(cox_zph)
}

cox_funk_cmp= cmpfun(cox_funk)

pheno= read.table(phenofile, h= T, sep='\t')

pheno[[ID]]= as.character(pheno[[ID]])
phenoID= pheno[[ID]]

d= read.table(IDfile, h=F)

pheno= pheno %>% filter(MOR_PID %in% d$V1)

d= as.character(d$V1)

colnames= c('variant','x1', 'x2', 'x3', 'x4')
colnames= append(colnames, d)

time_vec= pheno[, time_t]
outcome_vec= pheno[, outcome]
covars_m= as.matrix(pheno[,covars])

loglik_cox= summary(coxph( Surv( time_vec, outcome_vec)~ covars_m, na.action = na.omit ) )$loglik[2]

chunkSize <- 1000
sampleData <- read.table(gzfile(paste0(ds_folder, flist[1]), 'r'), h=F, nrows = 5, col.names= colnames, sep= '\t') #### ADD EXAMPLE DATA SET #######################################################
classes <- sapply(sampleData, class)
classes= gsub('integer','numeric',classes)
classes[2:5]= c('NULL','NULL','NULL','NULL')

# Function to run in parallel



funk= function(transactFile){

options(stringsAsFactors=FALSE)

library(dplyr)
library(tidyr)
library(survival)
library(compiler)

#print(args)
#ds_folder= args[1]  # dosage folder
#phenofile= args[2] # path to phenotype file
#IDfile= args[3] # file with sample IDs with same order as in dosage files
#ID= args[4] # ID name
#outfile= args[5] # path to output file
#time_t= args[6] # time variable name
#outcome= args[7] # outcome variable name
#covars= args[8:length(args)] # covariate variable name (multiple are accepted)

# read pheno file; each row is one participant, and column represents one variable
print(paste0(ds_folder,transactFile))
#falsecon = gzfile(paste0(ds_folder, transactFile))

#if (grepl('X',transactFile)){
#classes[1]= 'character'
#}

d= read.table(gzfile(paste0(ds_folder,transactFile)), nrows= 1, skip=0, header=F, fill = TRUE, sep="\t", col.names= colnames, colClasses= classes)
chr= strsplit(d[1,1], ':', fixed= T)[1]
out= paste0(outfile, chr)

if (basename(out) %in% list.files(dirname(outfile))){
return('Chromosome already analysed.')
}

con = gzfile(paste0(ds_folder, transactFile), 'r')

lskiped= 0
chunksize= 1000

tryCatch(repeat {
dataChunk= read.table(con, nrows= chunkSize, skip=0, header=F, fill = TRUE, sep="\t", col.names= colnames, colClasses= classes)
        genvars= dataChunk$variant
        if (length(genvars) == 0) break
	
	dataChunk= subset( dataChunk, select = -c(variant) )


        dataChunk= as.data.frame(t(dataChunk))
        dataChunk$id= gsub('X','',rownames(dataChunk))
	names(dataChunk)[1:length(genvars)]= genvars
        geno= inner_join(data.frame(MOR_PID= phenoID), dataChunk, by= c('MOR_PID' = 'id'))
	cox_coef= mclapply(geno[,-1], 2, cox_funk)
	cox_coef= do.call("rbind", cox_coef)
        cox_coef= data.frame(cox_coef)
        cox_coef[,2]= (-2*loglik_cox) - (-2*cox_coef[,2])
	cox_coef$variant= rownames(cox_coef) 
	
        write.table(cox_coef, out, append=T, row.names=F, col.names=F, quote=F, sep= '\t')
}, error= function(e) {})

close(con)
print(paste('Chromosome', chr,'finished!', sep= ' '))

}

mclapply(flist, funk, mc.cores= 5)
