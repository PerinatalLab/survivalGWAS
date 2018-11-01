#!/bin/bash

# The aim of this script is to obtain a list of samples in the order of the .vcf file and a dosage file from a multi-genome .vcf using bcftools. Schema of the dosage output file (for each .vcf): 
# CHR \t POS \t REF \t ALT \t SAMPLE_ID{1...X}

echo 'First we obtain a list of samples in the same order as in the .vcf file (we assume all files have the same order).'

cd ~/survivalGWAS/
zgrep -v '##' /mnt/work/hunt/vcf/moms/moms_CHR01_PID106764.vcf.gz | head -1 | cut -f10- > ./raw_data/samples_ID_moms

sed -i -e 's/\t/\n/g' ./raw_data/samples_ID_moms

echo 'Number of samples: ' $(cat ./raw_data/samples_ID_moms | wc -l)

zgrep -v '##' /mnt/work/hunt/vcf/fets/fets_CHR01_PID106764.vcf.gz | head -1 | cut -f10- > ./raw_data/samples_ID_fets

sed -i -e 's/\t/\n/g' ./raw_data/samples_ID_fets

echo 'Number of samples: ' $(cat ./raw_data/samples_ID_fets | wc -l)


echo 'Obtaining dosage from each .vcf file...'

mkdir /mnt/work/hunt/dosage/
mkdir /mnt/work/hunt/dosage/moms/
mkdir /mnt/work/hunt/dosage/fets/

bcftools=../soft/bcftools/bin/bcftools
vcffolder_moms=/mnt/work/hunt/vcf/moms/
vcffolder_fets=/mnt/work/hunt/vcf/fets/
outfolder_moms=/mnt/work/hunt/dosage/moms/
outfolder_fets=/mnt/work/hunt/dosage/fets/

if [ "$(ls "$outfolder_moms" | wc -l > 0)" ]; then
	echo 'Make sure the script has not been run previously and that the folder is empty...'
	exit
fi

parallel -j 7 '../soft/bcftools/bin/bcftools view -f PASS -i "MAF[0]>=0.01" {} | ../soft/bcftools/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n" -o /mnt/work/hunt/dosage/moms/dosage_$(basename {} _PID106764.vcf.gz); gzip /mnt/work/hunt/dosage/moms/dosage_$(basename {} _PID106764.vcf.gz)' ::: "$vcffolder_moms"*

echo 'Maternal dosages successfully extracted!'

if [ "$(ls "$outfolder_fets" | wc -l > 0)" ]; then
        echo 'Make sure the script has not been run previously and that the folder is empty...'
        exit
fi

parallel -j 7 '../soft/bcftools/bin/bcftools view -f PASS -i "MAF[0]>=0.01" {} | ../soft/bcftools/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n" -o /mnt/work/hunt/dosage/fets/dosage_$(basename {} _PID106764.vcf.gz); gzip /mnt/work/hunt/dosage/fets/dosage_$(basename {} _PID106764.vcf.gz)' ::: "$vcffolder_fets"*

echo 'Fetal dosages successfully extracted!'


