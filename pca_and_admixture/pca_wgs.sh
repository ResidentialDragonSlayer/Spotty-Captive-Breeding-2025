# code for PCA from WGS dataset

# drop all first-degree relatives
vcftools --gzvcf ../28samples_ref1_snps_pass_bial.vcf.gz --recode --recode-INFO-all \
--remove-indv NN114393 \
--remove-indv NNE_114394 \
--remove-indv SAMN32324417 \
--remove-indv SAMN32324419 \
--out 28samples_pca_admix

# drop Sunda clouded leopards
vcftools --gzvcf ../28samples_pca_admix.recode.vcf.gz --recode --recode-INFO-all --remove-indv SAMN32324430 --remove-indv SRR13774415 \
--out take3_admix_no_sunda

# filter for LD and HWE
plink --vcf take3_admix_no_sunda.recode.vcf.gz --double-id --allow-extra-chr --freq --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --hwe 0.001 --out pruned_no_sunda_take3


# calculate PCA
plink --vcf take3_admix_no_sunda.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract pruned_no_sunda_take3.prune.in \
--make-bed --pca --read-freq pruned_no_sunda_take3.frq --out no_sunda_take3
