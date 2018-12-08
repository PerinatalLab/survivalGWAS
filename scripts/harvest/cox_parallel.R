ptm <- proc.time()

library(dplyr)
library(tidyr)
library(survival)
library(parallel)
library(compiler)
library(data.table)

ds_folder= '/mnt/work/pol/harvest/dosage/maf01/fets/' #'/mnt/work/hunt/dosage/'  # dosage folder
phenofile=  '/mnt/work/pol/harvest/pheno/harvest_PROM_surv_fets_sens' #'/mnt/work/hunt/pheno/HUNT_PROM_surv_moms' # path to phenotype file
IDfile= '/mnt/work/pol/harvest/pheno/fets_IDs_harvest' #'./raw_data/samples_ID' # file with sample IDs with same order as in dosage files
ID= 'SentrixID_1' #'MOR_PID' # ID name
outfile= '/mnt/work/pol/gwas/results/surv/fets/PROM/HARVEST_PROM_fets_sent' #'/mnt/work/pol/gwas/res/survival/momsHUNT_PROM_chr' # path to output file
time_t= 'SVLEN_UL_DG' # time variable name
outcome= 'PROM' # outcome variable name
covars= c('BATCH', 'PC1','PC2','PC3', 'PC4', 'PC5', 'PC6','PARITY0') # covariate variable name (multiple are accepted)

options(stringsAsFactors=FALSE)

flist= list.files(path= ds_folder, pattern='.gz')

pheno= read.table(phenofile, h= T, sep='\t')

pheno[[ID]]= as.character(pheno[[ID]])
phenoID= pheno[[ID]]

d= read.table(IDfile, h=F)

pheno= pheno %>% filter(SentrixID_1 %in% d$V1)

d= as.character(d$V1)

colnames= c('variant','x1', 'x2', 'x3', 'x4')
colnames= append(colnames, d)

time_vec= pheno[, time_t]
outcome_vec= pheno[, outcome]
covars_m= as.matrix(pheno[,covars])

sampleData <- read.table(gzfile(paste0(ds_folder, flist[1]), 'r'), h=F, nrows = 5, col.names= colnames, sep= '\t') #### ADD EXAMPLE DATA SET #######################################################
classes= lapply(sampleData, class)
classes= gsub('integer','numeric',classes)

# Function to run in parallel

funk= function(transactFile){

if (grepl('X', transactFile)){
classes[2]= 'character'
}

options(stringsAsFactors=FALSE)

# read pheno file; each row is one participant, and column represents one variable
print(paste0(ds_folder,transactFile))

d= read.table(gzfile(paste0(ds_folder,transactFile)), nrows= 1, skip=0, header=F, fill = TRUE, sep="\t", col.names= colnames, colClasses= classes)
chr= unlist(strsplit(as.character(d[1,1]), ':', fixed= T))[1]

out= paste0(outfile, chr)

if (basename(out) %in% list.files(dirname(outfile))){
return('Chromosome already analysed.')
}

con = gzfile(paste0( ds_folder, transactFile))

lskiped= 0
chunkSize= 150
open(con)

repeat {
block.text= readLines(con, n= chunkSize)

if (length(block.text) == 0){ # if there's nothing left, leave the loop
        break
}

dataChunk= fread(text=block.text, sep="\t", col.names= colnames, colClasses= classes)
        genvars= dataChunk$variant
        if (length(genvars) == 0) break
	
	dataChunk= subset( dataChunk, select = -c(variant,x1,x2,x3,x4))


        dataChunk= as.data.frame(t(dataChunk))
        dataChunk$id= gsub('X','',rownames(dataChunk))
	names(dataChunk)[1:length(genvars)]= genvars
        geno= inner_join(pheno, dataChunk, by= c('SentrixID_1' = 'id'))
	cox_coef= mclapply(names(geno[,-c(1:dim(pheno)[2])]), mc.cores= 3, function(snp){cox_coef= coxph(Surv( geno$SVLEN_UL_DG, geno$PROM)~ geno[,snp] + geno$BATCH + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6, na.action = na.omit)
	coef = summary( cox_coef)$coefficients[1,1]
	sd= summary(cox_coef)$coefficient[1,3]
	n= summary(cox_coef)$n
	pvalue= summary(cox_coef)$coefficient[1,5]
	txt = sprintf( "%s\t%e\t%e\t%e\t%e\n", snp, n, coef, sd, pvalue)
	cat(txt, file= out, append= T)
    
}
)

}

close(con)
print(paste('Chromosome', chr,'finished!', sep= ' '))

}

mclapply(flist, mc.cores= 15, funk)
