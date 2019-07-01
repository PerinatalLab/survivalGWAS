ptm <- proc.time()

library(dplyr)
library(tidyr)
library(survival)
library(parallel)
library(compiler)
library(data.table)

options(stringsAsFactors=FALSE)

infile= snakemake@input[[1]]
phenofile= snakemake@input[[2]]
outfile= snakemake@ouput[[1]]

pheno= read.table(phenofile, h= T, sep='\t')

if (grepl('harvest', phenofile)){
pheno= mutate(pheno, spont= as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                INDUKSJON_PROSTAGLANDIN==0 &
                INDUKSJON_ANNET==0 &
                INDUKSJON_OXYTOCIN==0 &
                INDUKSJON_AMNIOTOMI==0),
                PARITY0= as.numeric(PARITET_5==0),
		PROM= as.numeric(!is.na(VANNAVGANG)))
} else if (!grepl('harvest', phenofile)){
pheno= mutate(pheno, spont= as.numeric(FSTART=='Spontan' | FSTART== '' & (KSNITT=='' | KSNITT== 'Uspesifisert' | KSNITT== 'Akutt keisersnitt') &
                INDUKSJON_PROSTAGLANDIN=='Nei' &
                INDUKSJON_ANNET=='Nei' &
                INDUKSJON_OXYTOCIN=='Nei' &
                INDUKSJON_AMNIOTOMI=='Nei'),
                PARITY0= as.numeric(PARITET_5=='0 (førstegangsfødende)'),
		PROM= as.numeric(!is.na(VANNAVGANG)))
names(pheno)[names(pheno) == 'SentrixID'] <- 'SentrixID_1'
}


if (grepl('spont', outfile)){
funk= function(block.text){
dataChunk= fread(text=block.text, sep="\t", col.names= colnames)
        genvars= paste(dataChunk$CHR, dataChunk$BP, sep=':')

        if (length(genvars) == 0) break

        dataChunk= subset( dataChunk, select = -c(CHR,BP))


        dataChunk= as.data.frame(t(dataChunk))
        dataChunk$id= gsub('X','',rownames(dataChunk))
        names(dataChunk)[1:length(genvars)]= genvars
        geno= inner_join(pheno, dataChunk, by= c('SentrixID_1' = 'id'))
        cox_coef= lapply(names(geno[,-c(1:dim(pheno)[2])]), function(snp){cox_coef= coxph(Surv( geno$SVLEN_UL_DG, geno$spont)~ geno[,snp] + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6, na.action = na.omit)
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

} else if (grepl('PROM', outfile)){

funk= function(block.text){
dataChunk= fread(text=block.text, sep="\t", col.names= colnames)
        genvars= paste(dataChunk$CHR, dataChunk$BP, dataChunk$REF, dataChunk$ALT sep=':')

        if (length(genvars) == 0) break

        dataChunk= subset( dataChunk, select = -c(CHR,BP,REF,ALT))


        dataChunk= as.data.frame(t(dataChunk))
        dataChunk$id= gsub('X','',rownames(dataChunk))
        names(dataChunk)[1:length(genvars)]= genvars
        geno= inner_join(pheno, dataChunk, by= c('SentrixID_1' = 'id'))
        cox_coef= lapply(names(geno[,-c(1:dim(pheno)[2])]), function(snp){cox_coef= coxph(Surv( geno$SVLEN_UL_DG, geno$PROM)~ geno[,snp] + geno$PARITY0 + geno$PC1 + geno$PC2 + geno$PC3 + geno$PC4 + geno$PC5 + geno$PC6, na.action = na.omit)
        coef = summary(cox_coef)$coefficients[1,1]
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

}

colnames= readLines(gzfile(infile), n= 1)
colnames= unlist(strsplit(colnames, '\t'))

con= gzfile(infile)
open(con)

repeat {
block.text= readLines(con, n= 200)

if (length(block.text) == 0){ # if there's nothing left, leave the loop
        break
}

funk(block.text)

}

close(con)

