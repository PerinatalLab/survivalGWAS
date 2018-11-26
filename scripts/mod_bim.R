#!/usr/bin/R

options(stringsAsFactors= F)

print('Modifying .bim file')

flist= list.files('/mnt/hdd/common/survivalGWAS/results/VEGAS2/1KG/plink_files/', '*.bim')


for (i in flist){
d= read.table(i, h=F)
d[(d$V5>d$V6), c('V5','V6')]= d[(d$V5>d$V6), c('V6','V5')]
d$V2= paste(d$V1, d$V4, d$V5, d$V6, sep=':')
write.table(d, i, row.names=F, col.names=F, quote=F, sep= '\t')
}

print('All .bim file modified')
