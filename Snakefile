import pandas as pd
import numpy as np

cohort_nms= ['harvestm12', 'harvestm24', 'rotterdam1', 'rotterdam2', 'normentfeb', 'normentmay']
sample_nms= ['maternal', 'fetal']
CHR_nms= ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
pheno_nms= ['spont', 'PROM']

rule all:
	'Collect all inputs.'
	input:
		expand('/mnt/work/pol/survivalGWAS/results/{cohort}/results_{pheno}_{sample}.txt', cohort= cohort_nms, pheno= pheno_nms, sample= sample_nms)


rule extract_vcf_samples:
        'Extract samples id included in the VCF file, for each cohort.'
        input:
                '/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m12/1.vcf.gz',
                '/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m24/1.vcf.gz',
                '/mnt/archive/ROTTERDAM1/delivery-fhi/data/imputed/1.vcf.gz',
                '/mnt/archive/ROTTERDAM2/delivery-fhi/data/imputed/1.vcf.gz',
                '/mnt/archive/NORMENT1/delivery-fhi/data/imputed/feb18/1.vcf.gz',
                '/mnt/archive/NORMENT1/delivery-fhi/data/imputed/may16/1.vcf.gz'
        output:
                temp('/mnt/work/pol/survivalGWAS/raw_data/{cohort}/vcf_ids')
        run:
                if 'harvestm12' == wildcards.cohort: vcf= input[0]
                if 'harvestm24' == wildcards.cohort: vcf= input[1]
                if 'rotterdam1' == wildcards.cohort: vcf= input[2]
                if 'rotterdam2' == wildcards.cohort: vcf= input[3]
                if 'normentfeb' == wildcards.cohort: vcf= input[4]
                if 'normentmay' == wildcards.cohort: vcf= input[5]
                shell("set +o pipefail; zgrep -v '##' {vcf} | head -1 | cut -f10- | sed 's/\\t/\\n/g'  > {output[0]} ")

rule extract_samples:
        'Samples for filtering VCF files.'
        input:
                '/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv',
                '/mnt/work/pol/survivalGWAS/raw_data/{cohort}/vcf_ids'
        output:
                '/mnt/work/pol/survivalGWAS/raw_data/{cohort}/{sample}_toextract'
        run:
                if 'harvest' in wildcards.cohort:
                        d= pd.read_csv(input[0], delim_whitespace= True, header= 0)
                        Sentrix= 'SentrixID_1'
                if 'harvest' not in wildcards.cohort:
                        d= pd.read_csv(input[0], delim_whitespace= True, header= 0)
                        Sentrix= 'SentrixID'
                if 'maternal' in wildcards.sample:
                        d= d.loc[d.Role=='Mother', :]
                if 'fetal' in wildcards.sample:
                        d= d.loc[d.Role=='Child', :]
                x= [line.strip() for line in open(input[1], 'r')]
                d= d.loc[d[Sentrix].isin(x), :]
                d.drop_duplicates(subset= [Sentrix], inplace= True)
                d.to_csv(output[0], header= False, columns= [Sentrix], index= False, sep= '\t')

rule phenofile:
        'Merge all data necessary to create a phenotype file for spontaneous delivery and PROM.'
        input:
                '/mnt/work/pol/{cohort}/pheno/{cohort}_mfr.csv',
                '/mnt/work/pol/{cohort}/pheno/{cohort}_linkage.csv',
                '/mnt/work/pol/{cohort}/pca/{cohort}_pca.txt',
                '/mnt/work/pol/{cohort}/relatedness/all_{cohort}.kin0',
                '/mnt/work/pol/{cohort}/pheno/flag_list.txt',
		'/mnt/work/pol/survivalGWAS/raw_data/{cohort}/{sample}_toextract',
		'/mnt/work/pol/harvestm12/pca/all_pca_outliers_hapmap.txt'
        output:
                '/mnt/work/pol/survivalGWAS/pheno/{cohort}/pheno_{sample}.txt'
        script:
                'scripts/pheno_file.py'

rule survival_analysis:
	''
	input:
		'/mnt/work/pol/survivalGWAS/pheno/{cohort}/pheno_maternal.txt',
		'/mnt/work/pol/survivalGWAS/pheno/{cohort}/pheno_fetal.txt',
		'/mnt/work/pol/survivalGWAS/raw_data/{cohort}/vcf_ids',
		'/mnt/work/pol/survivalGWAS/raw_data/{cohort}/maternal_toextract',
		'/mnt/work/pol/survivalGWAS/raw_data/{cohort}/fetal_toextract',
		'/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m12/{CHR}.vcf.gz',
		'/mnt/archive/HARVEST/delivery-fhi/data/imputed/imputed_m24/{CHR}.vcf.gz',
		'/mnt/archive/ROTTERDAM1/delivery-fhi/data/imputed/{CHR}.vcf.gz',
		'/mnt/archive/ROTTERDAM2/delivery-fhi/data/imputed/{CHR}.vcf.gz',
		'/mnt/archive/NORMENT1/delivery-fhi/data/imputed/feb18/{CHR}.vcf.gz',
		'/mnt/archive/NORMENT1/delivery-fhi/data/imputed/may16/{CHR}.vcf.gz'
	output:
		temp('/mnt/work/pol/survivalGWAS/results/{cohort}/spont_maternal_{CHR}_results.txt'),
		temp('/mnt/work/pol/survivalGWAS/results/{cohort}/PROM_maternal_{CHR}_results.txt'),
		temp('/mnt/work/pol/survivalGWAS/results/{cohort}/spont_fetal_{CHR}_results.txt'),
		temp('/mnt/work/pol/survivalGWAS/results/{cohort}/PROM_fetal_{CHR}_results.txt')
	script:
		'scripts/cox.R'

rule concat_results:
	'Concat results from GWAS.'
	input:
		expand('/mnt/work/pol/survivalGWAS/results/{{cohort}}/{{pheno}}_{{sample}}_{CHR}_results.txt', CHR= CHR_nms)
	output:
		'/mnt/work/pol/survivalGWAS/results/{cohort}/results_{pheno}_{sample}.txt'
	shell:
		'cat {input} > {output[0]}'

