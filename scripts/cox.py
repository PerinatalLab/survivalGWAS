#!/usr/bin/python3

ds_folder= '/mnt/work/hunt/dosage/moms/'  # dosage folder
phenofile=  '/mnt/work/hunt/pheno/HUNT_PROM_surv_moms' # path to phenotype file
IDfile= '/mnt/work/hunt/pheno/maternal_ids' # file with sample IDs with same order as in dosage files
ID= 'MOR_PID' # ID name
outfile= '/mnt/work/pol/gwas/res/survival/momsHUNT_PROM_chr' # path to output file
time_t= 'SVLEN_DG' # time variable name
outcome= 'PROM' # outcome variable name
covars= ['FAAR', 'PC1','PC2','PC3', 'PC4', 'PC5', 'PC6','PARITY0'] # covariate variable name (multiple are accepted)

import gzip
import glob
import sys
import pandas as pd
import multiprocessing as mp
import argparse

parser = argparse.ArgumentParser(description='Survival analysis in GWAS. Inputs a dosage folder and phenotype file, outputs results (beta, SD, p-value, ...')
parser.add_argument('--ds_folder', help='Folder and path storing dosage.', default= None, action='store', required=True)
parser.add_argument('--phenofile', help='Folder and path of phenotype file.', default= None, action='store', required=True)
parser.add_argument('--IDfile', help='Folder and path sample IDs.', default= None, action='store', required=True)
parser.add_argument('--ID', help='Phenotype sample ID name.', default= 'MOR_PID', action='store', required=True)
parser.add_argument('--outfile', help='Path and output file name.', default= 'survival_chr', action='store', required=True)
#parser.add_argument('--time_t', help='Time variable name.', default= SVLEN_DG, action='store', required=True)
#parser.add_argument('--outcome', help='Outcome variable name.', default= PROM, action='store', required=True)
#parser.add_argument('--covars', help='Covariate/s name/s.', nargs='+', default= None, action='store', required=True)

args = parser.parse_args()


# To call R from python3

from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2 import rinterface
#from rpy2.robjects import r, pandas2ri
from rpy2.robjects import Formula
pandas2ri.activate()
survival = importr('survival')
coxph = survival.coxph
Surv = survival.Surv
stats = importr('stats')
base = importr('base')

flist= glob.glob(args.ds_folder+'*.gz')

pheno= pd.read_csv(args.phenofile, sep='\t')

d= pd.read_csv(args.IDfile, header= None)

pheno= pheno[pheno[args.ID].isin(d[0])]

pheno.loc[:, args.ID]= pheno.loc[:, args.ID].astype(str)

cnames= ['variant','x1', 'x2', 'x3', 'x4']
cnames.extend(d[0].astype(str))

loglik_cox= coxph(Formula("Surv(SVLEN_DG, PROM) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + FAAR + PARITY0"), data= pheno, method= 'efron')
loglik_cox= loglik_cox.rx2('loglik')[1]

sampleData= pd.read_csv(flist[22], compression= 'gzip', names=cnames, nrows = 5, sep= '\t') 
classes= sampleData.dtypes
classes= classes.replace('int64','float64')

def processChunk ( nl ):
	new_list = nl.replace('\n', '').split('\t')
	d= pd.DataFrame(new_list, columns= ['var'])
	genvars= d.iloc[0,0]
	d[args.ID]= cnames
	d.drop(d.index[0:5], inplace= True)
	d['var']= d['var'].astype('float64')
	d= d.merge(pheno, how='right', on= args.ID)
	d.drop([args.ID], axis=1, inplace= True)
	try:
		ch= coxph(Formula("Surv(SVLEN_DG, PROM) ~ var + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + FAAR + PARITY0"), data = d, method = 'efron')
	except rpy2.rinterface.RRuntimeError:
		coeffs= NA
		beta= NA
		se= NA
		pvalue= NA
		loglik= NA
	coeffs = base.summary(ch).rx2('coefficients')
	df= pd.DataFrame(pandas2ri.ri2py(coeffs), index=coeffs.names[0], columns=coeffs.names[1])
	beta= df.iloc[0,0]
	se= df.iloc[0,2]
	pvalue= df.iloc[0,4]
	loglik= (-2*loglik_cox) - (-2*ch.rx2('loglik')[1])
	#print('%s\t%s\t%s\t%s\t%s' % (genvars, beta, se, pvalue, loglik), sep='', end='\n', file= open(out, "a"))
	return str(genvars) + '\t' + str(beta) + '\t' + str(se) + '\t' + str(pvalue) + '\t' + str(loglik)

pool = mp.Pool(25)

for i in flist:
	d= pd.read_csv(i, compression= 'gzip', nrows= 1, skiprows=0, names= cnames, sep="\t")
	CHROM= d.iloc[0,0].split(':')[0]
	out= args.outfile+CHROM
	with gzip.open(i, 'rt') as f:
		for r in pool.imap_unordered(processChunk, (line for line in f), chunksize= 100):
			print(r, file= open(out, 'a'), end= '\n', sep= '')


