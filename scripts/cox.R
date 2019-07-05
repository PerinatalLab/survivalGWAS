


library(dplyr)
library(survival)
library(data.table)

options(stringsAsFactors=FALSE)

funk= function(pheno, geno, outcome, outfile){
        cox_coef= lapply(names(geno[,-c(1:ncol(pheno))]), function(snp) {
	cox_coef= coxph(Surv(geno$SVLEN_UL_DG, geno[, outcome])~ geno[,snp] + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6, na.action = na.omit)
        coef = summary( cox_coef)$coefficients[1,1]
        sd= summary(cox_coef)$coefficient[1,3]
        n= summary(cox_coef)$n
        pvalue= summary(cox_coef)$coefficient[1,5]
        zph= cox.zph(cox_coef)
        correlation= zph$table[1, 1]
        corr_pvalue= zph$table[1, 3]
        txt = sprintf( "%s\t%e\t%e\t%e\t%e\t%e\t%e\n", snp, n, coef, sd, pvalue, correlation, corr_pvalue)
cat(txt, file= outfile, append= T)
}
)
}


moms_phenofile= snakemake@input[[1]]
fets_phenofile= snakemake@input[[2]]
ids=  snakemake@input[[3]]
moms_ids=  snakemake@input[[4]]
fets_ids=  snakemake@input[[5]]



outfile_spont_moms= snakemake@output[[1]]
outfile_PROM_moms= snakemake@output[[2]]
outfile_spont_fets= snakemake@output[[3]]
outfile_PROM_fets= snakemake@output[[4]]

if (grepl('harvestm12', moms_phenofile)){
infile= snakemake@input[[6]]
} else if (grepl('harvestm24', moms_phenofile)) {
infile= snakemake@input[[7]]
} else if (grepl('rotterdam1', moms_phenofile)) {
infile= snakemake@input[[8]]
} else if (grepl('rotterdam2', moms_phenofile)) {
infile= snakemake@input[[9]]
} else if (grepl('normentfeb', moms_phenofile)) {
infile= snakemake@input[[10]]
} else if (grepl('normentmay', moms_phenofile)) {
infile= snakemake@input[[11]]
}


moms_pheno= fread(moms_phenofile)
fets_pheno= fread(fets_phenofile)
ids= readLines(ids)
moms_ids= readLines(moms_ids)
fets_ids= readLines(fets_ids)


format_pheno= function(pheno) {
if (grepl('harvest', moms_phenofile)){
pheno= pheno %>% mutate(spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 &
                INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 &
                INDUKSJON_AMNIOTOMI==0),
                PARITY0= as.numeric(PARITET_5==0),
                PROM= as.numeric(!is.na(VANNAVGANG)))
} else if (!grepl('harvest', moms_phenofile)){
pheno= pheno %>% mutate(spont= as.numeric(FSTART=='Spontan' | FSTART== '' & (KSNITT=='' | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' &
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' &
                INDUKSJON_AMNIOTOMI=='Nei'),
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'),
                PROM= as.numeric(!is.na(VANNAVGANG)))
names(pheno)[names(pheno) == 'SentrixID'] <- 'SentrixID_1'
}
return(pheno)
}

moms_pheno= format_pheno(moms_pheno)
fets_pheno= format_pheno(fets_pheno)

con = file(infile, "r")

repeat {
z= scan(con, n=100, what= 'character', sep= '\n', comment.char= '#')

if ( length(z) == 0 ) {
      break
}

df= fread(text= z, sep= ':')
x= Filter(is.numeric, df)

varnames= apply(df[,1], 1, function(x) paste(unlist(strsplit(as.character(x), '\t'))[c(1,2,4,5)], collapse= ':'))

#varnames= lapply(z, function(x) substr(x, 1, 100))
#varnames= lapply(varnames, function(x) paste(unlist(strsplit(x, '\t'))[c(1,2,4,5)], collapse= ':'))

#z= lapply(z, function(x) as.numeric(sapply(strsplit(unlist(strsplit(x, '\t'))[-2:-9], ':'), '[', 3)))

#x= data.frame(z)

x= as.data.frame(t(x))
names(x)= unlist(varnames)
x= x[which(apply(!is.na(x), 1, all)),]

x= cbind(id= ids, x)




geno_moms= inner_join(moms_pheno, x, by= c('SentrixID_1' = 'id'))
geno_fets= inner_join(fets_pheno, x, by= c('SentrixID_1' = 'id'))

funk(moms_pheno, geno_moms, 'spont', outfile_spont_moms)
funk(moms_pheno, geno_moms, 'PROM', outfile_PROM_moms)

funk(fets_pheno, geno_fets, 'spont', outfile_spont_fets)
funk(fets_pheno, geno_fets, 'PROM', outfile_PROM_fets)

}

