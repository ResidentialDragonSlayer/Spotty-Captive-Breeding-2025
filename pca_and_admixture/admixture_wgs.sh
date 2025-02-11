# ADMIXTURE commands for WGS dataset

# filter for LD, HWE, genotype frequency, and minor allele frequency

plink --vcf take3_admix_no_sunda.recode.vcf.gz \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.5 \
--hwe 0.01 \
--geno 0.1 \
--mind 0.15 \
--maf 0.05 \
--out pruned_no_sunda_take3

# create bed/bim files
plink --vcf take3_admix_no_sunda.recode.vcf.gz \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--extract pruned_no_sunda_take3.prune.in \
--make-bed \
--out no_sunda_take3

awk '{$1=0;print $0}' no_sunda_take3.bim > no_sunda_take3.bim.tmp
mv no_sunda_take3.bim.tmp no_sunda_take3.bim

# ADMIXTURE command
# NOTE: this command was run as a job array with task IDs 2-6, thus k=2 through k=6 was evaluated in total
admixture --cv no_sunda_take3.bed $SGE_TASK_ID > no_sunda_all_take3_log${SGE_TASK_ID}.out

# format CV error and ID files
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > cloudy.cv.error
awk '{split($1,name,"."); print $1,name[2]}' cloudy.nosex > cloudy.list
