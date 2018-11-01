#!/bin/bash

# Script to calculate kinship in maternal and fetal samples. Requires PLINK v1.90b6.6 and PLINK v2.00a1LM

mkdir /mnt/work/pol/harvest/
mkdir /mnt/work/pol/harvest/genotyped/
mkdir /mnt/work/pol/harvest/relatedness/
mkdir ./raw_data/

# 1. Obtain a list of maternal and fetal samples

zgrep -v '##' /mnt/archive/jjuod/merging/moms_1.vcf.gz | head -1 | cut -f10- > ./raw_data/samples_ID_moms

sed -i -e 's/\t/\n/g' ./raw_data/samples_ID_moms

zgrep -v '##' /mnt/archive/jjuod/merging/fets_1.vcf.gz | head -1 | cut -f10- > ./raw_data/samples_ID_fets

sed -i -e 's/\t/\n/g' ./raw_data/samples_ID_fets

# 2. Merge batches (m12 and m24)

../soft/plink --bfile /mnt/archive/HARVEST/delivery-fhi/data/genotyped/m12/m12-genotyped --bmerge /mnt/archive/HARVEST/delivery-fhi/data/genotyped/m24/m24-genotyped --out /mnt/work/pol/harvest/genotyped/m24-m12-genotyped

# 3. Run KING 

../soft/plink2 --bfile /mnt/work/pol/harvest/genotyped/m24-m12-genotyped --keep-fam ./raw_data/samples_ID_moms --make-king-table --king-table-filter 0.03125 --out /mnt/work/pol/harvest/relatedness/moms

../soft/plink2 --bfile /mnt/work/pol/harvest/genotyped/m24-m12-genotyped --keep-fam ./raw_data/samples_ID_fets --make-king-table --king-table-filter 0.03125 --out /mnt/work/pol/harvest/relatedness/fets
