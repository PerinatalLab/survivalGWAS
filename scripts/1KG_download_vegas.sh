!/usr/bin/bash

# Script to download 1KG vcf files, convert to plink .bed .bim .fam files, and calculate LD score. 

# First step, download cM map from 1KG

cd /mnt/hdd/common/survivalGWAS/results/VEGAS2/1KG/

while read i; do
        wget "$i"
        ~/software/plink2 --vcf ./*.vcf.gz --keep ./EUR_tokeep2 --maf 0.01 --make-bed --out ./plink_files/$(basename "$i" .vcf.gz)
        rm ./*.vcf.gz        
done < ~/survivalGWAS/raw_data/1KG_links

Rscript ~/survivalGWAS/scripts/mod_bim.R
