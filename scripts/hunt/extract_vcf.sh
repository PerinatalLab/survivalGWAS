#!/bin/bash

echo 'Script to extract maternal and fetal vcf'

mkdir /mnt/work/hunt/vcf/
mkdir /mnt/work/hunt/vcf/moms/
mkdir /mnt/work/hunt/vcf/fets/

parallel -j 7 '../soft/bcftools/bin/bcftools view -Oz -S /mnt/work/pol/HUNT_PCA/mother_samples {} >> /mnt/work/hunt/vcf/moms/moms_$(basename {})' ::: /mnt/archive/hunt/genotypes/vcf/*

echo 'Maternal vcf extracted'

parallel -j 7 '../soft/bcftools/bin/bcftools view -Oz -S /mnt/work/pol/HUNT_PCA/fetal_samples {} >> /mnt/work/hunt/vcf/fets/fets_$(basename {})' ::: /mnt/archive/hunt/genotypes/vcf/*

echo 'Fetal vcf exctracted'





