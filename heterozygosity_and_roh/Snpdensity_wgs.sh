#!/bin/bash

VCF="28samples_ref1_snps_pass_bial.vcf.gz"
SAMPLES=$(bcftools query -l $VCF)

echo "Make individual VCFs"
for sample in $SAMPLES; do
    echo "Processing: $sample"
    vcftools --gzvcf $VCF --recode --recode-INFO-all --indv $sample --out $sample
done

echo "Heterozygous sites only"
for sample in $SAMPLES; do
    echo "Filtering heterozygous sites for: $sample"
    vcftools --vcf ${sample}.recode.vcf --recode --maf 0.1 --out ${sample}_het
done

echo "Calculate SNP density"
for sample in $SAMPLES; do
    echo "Calculating SNP density for: $sample"
    vcftools --vcf ${sample}_het.recode.vcf --SNPdensity 1000000 --out ${sample}_hetsites

    if [[ -f ${sample}_hetsites.snpden ]]; then
        awk -v sample=$sample 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' ${sample}_hetsites.snpden > ${sample}>    else
        echo "Warning: ${sample}_hetsites.snpden not found!"
    fi
done

# all SNP density files
if ls *_id.snpden 1> /dev/null 2>&1; then
    tail -q -n +2 *_id.snpden > all_hetsites2.snpden
    echo "Merged SNP density files into all_hetsites2.snpden"
else
    echo "No SNP density files found to merge!"
fi
