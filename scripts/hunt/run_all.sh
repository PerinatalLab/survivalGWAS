#!/bin/bash

# Rscript scripts/pheno_HUNT.R

mkdir /mnt/work/pol/gwas/res/survival/moms/
mkdir /mnt/work/pol/gwas/res/survival/fets/

# Rscript --vanilla scripts/cox_parallel.R /mnt/work/hunt/dosage/moms/ /mnt/work/hunt/pheno/HUNT_PROM_surv_moms ./raw_data/samples_ID_moms MOR_PID /mnt/work/pol/gwas/res/survival/moms/momsHUNT_PROM_chr

Rscript --vanilla scripts/cox_parallel.R /mnt/work/hunt/dosage/moms/ /mnt/work/hunt/pheno/HUNT_PROM_surv_moms ./raw_data/samples_ID_moms MOR_PID /mnt/work/pol/gwas/res/survival/moms/momsHUNT_PROM_chr

# Rscript --vanilla scripts/cox_parallel.R /mnt/work/hunt/dosage/fets/ /mnt/work/hunt/pheno/HUNT_PROM_surv_fets ./raw_data/samples_ID_fets BARN_PID /mnt/work/pol/gwas/res/survival/fets/fetsHUNT_PROM_chr

# Rscript --vanilla scripts/cox_parallel.R /mnt/work/hunt/dosage/fets/ /mnt/work/hunt/pheno/HUNT_PROM_surv_fets_sens ./raw_data/samples_ID_fets BARN_PID /mnt/work/pol/gwas/res/survival/fets/fetsHUNT_PROM_sens_chr


echo 'End!!!'
