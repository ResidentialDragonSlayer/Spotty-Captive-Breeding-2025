# 🧬 Variant Calling Pipeline

This repository contains a full pipeline for variant calling, filtration, and annotation of genomic data. It is designed to process high-throughput sequencing reads using standard bioinformatics tools.

---

## 📋 Pipeline Overview

The pipeline consists of three main stages:

1. **Step 1: Preprocessing**
   - Quality control with `FastQC` and `MultiQC`
   - Read trimming with `fastp`
   - Alignment to two reference genomes using `BWA-MEM2`
   - Duplicate marking and metrics collection with `Picard` and `samtools`

2. **Step 2: Variant Calling**
   - Variant calling with `GATK HaplotypeCaller` for each reference
   - Generation of GVCF files per sample

3. **Step 3: Joint Genotyping and Filtering**
   - Combining GVCFs with `GATK CombineGVCFs`
   - Genotyping with `GATK GenotypeGVCFs`
   - SNP selection and hard filtering
   - Extraction of high-quality biallelic SNPs and summary statistics

4. **Step 4: Annotation**
   - Sorting and indexing GFF files
   - Variant annotation with 

---

## ⚙️ SLURM Batch Script

Each step is implemented as a separate SLURM batch script for parallel execution across samples or references.

---
### 🔹 Step 1: Preprocessing and Alignment

```bash
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=14-0
#SBATCH --array=1-28%12

# 📁 Base directories
samplesheet="path/to/28samples.txt"
reads_dir="path/to/reads"
ref1="path/to/reference1.fasta"
ref2="path/to/reference2.fasta"
picard_jar="path/to/picard.jar"

# 📁 Unified output base
base_output="path/to/output"
bam_dir="${base_output}/bam"
md_dir="${base_output}/marked_duplicates"
metrics_dir="${base_output}/metrics"

# 📁 Reference-specific subfolders
ref1_tag="ref1"
ref2_tag="ref2"

# 🔧 Thread count
threads=8

# 🧬 Sample names
r1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samplesheet" | awk '{print $1}')
r2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samplesheet" | awk '{print $2}')

# 📁 QC and trimming directories
qc_raw_dir="${base_output}/qc/raw"
qc_trimmed_dir="${base_output}/qc/trimmed"
fastp_json_dir="${base_output}/fastp/json"
fastp_html_dir="${base_output}/fastp/html"
multiqc_dir="${base_output}/qc/multiqc"

# 📊 Raw read QC
fastqc "${reads_dir}/${r1}_1.fastq.gz" "${reads_dir}/${r1}_2.fastq.gz" -t 2 -o "$qc_raw_dir"

# ✂️ Read trimming
fastp \
  --in1 "${reads_dir}/${r1}_1.fastq.gz" \
  --in2 "${reads_dir}/${r1}_2.fastq.gz" \
  --out1 "${reads_dir}/${r1}_1.trimmed.fastq.gz" \
  --out2 "${reads_dir}/${r1}_2.trimmed.fastq.gz" \
  --detect_adapter_for_pe \
  -l 20 -q 20 \
  --trim_poly_g \
  --overrepresentation_analysis \
  --json "${fastp_json_dir}/${r1}.json" \
  --html "${fastp_html_dir}/${r1}.html" \
  --thread "$threads" \
  --report_title "${r1}"

# 📊 Trimmed read QC
fastqc "${reads_dir}/${r1}_1.trimmed.fastq.gz" "${reads_dir}/${r1}_2.trimmed.fastq.gz" -o "$qc_trimmed_dir"

# 📊 MultiQC summary
multiqc "$qc_raw_dir" "$qc_trimmed_dir" "$fastp_json_dir" -o "$multiqc_dir"

# 🧷 Alignment to ref1
bwa mem -t "$threads" -R "@RG\tID:${r1}\tPL:ILLUMINA\tPU:${r1}\tLB:${r1}\tSM:${r2}" \
    "$ref1" "${reads_dir}/${r1}_1.trimmed.fastq.gz" "${reads_dir}/${r1}_2.trimmed.fastq.gz" \
    | samtools view -u - | samtools sort -@ "$threads" -o "${bam_dir}/${ref1_tag}/${r2}.bam"
samtools index -@ "$threads" "${bam_dir}/${ref1_tag}/${r2}.bam"

# 🧷 Alignment to ref2
bwa mem -t "$threads" -R "@RG\tID:${r1}\tPL:ILLUMINA\tPU:${r1}\tLB:${r1}\tSM:${r2}" \
    "$ref2" "${reads_dir}/${r1}_1.trimmed.fastq.gz" "${reads_dir}/${r1}_2.trimmed.fastq.gz" \
    | samtools view -u - | samtools sort -@ "$threads" -o "${bam_dir}/${ref2_tag}/${r2}.bam"
samtools index -@ "$threads" "${bam_dir}/${ref2_tag}/${r2}.bam"

# 🧼 Mark duplicates
java -jar "$picard_jar" MarkDuplicates \
    -I "${bam_dir}/${ref1_tag}/${r2}.bam" \
    -O "${md_dir}/${ref1_tag}/${r2}.bam" \
    -M "${metrics_dir}/${ref1_tag}/${r2}_MD_metrics.txt" \
    --CREATE_INDEX true

java -jar "$picard_jar" MarkDuplicates \
    -I "${bam_dir}/${ref2_tag}/${r2}.bam" \
    -O "${md_dir}/${ref2_tag}/${r2}.bam" \
    -M "${metrics_dir}/${ref2_tag}/${r2}_MD_metrics.txt" \
    --CREATE_INDEX true

# 📊 Metrics collection
java -jar "$picard_jar" CollectWgsMetrics \
    -I "${bam_dir}/${ref1_tag}/${r2}.bam" \
    -O "${metrics_dir}/${ref1_tag}/${r2}_wgs_metrics.txt" \
    -R "$ref1"

java -jar "$picard_jar" CollectAlignmentSummaryMetrics \
    -I "${bam_dir}/${ref1_tag}/${r2}.bam" \
    -O "${metrics_dir}/${ref1_tag}/${r2}_alignment_metrics.txt" \
    -R "$ref1"

java -jar "$picard_jar" CollectWgsMetrics \
    -I "${bam_dir}/${ref2_tag}/${r2}.bam" \
    -O "${metrics_dir}/${ref2_tag}/${r2}_wgs_metrics.txt" \
    -R "$ref2"

java -jar "$picard_jar" CollectAlignmentSummaryMetrics \
    -I "${bam_dir}/${ref2_tag}/${r2}.bam" \
    -O "${metrics_dir}/${ref2_tag}/${r2}_alignment_metrics.txt" \
    -R "$ref2"

```


### 🔹 Step 2: Variant Calling with GATK

```bash
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=14-0
#SBATCH --array=1-28%12

# 📁 Base directories
samplesheet="path/to/28samples.txt"
gatk_path="path/to/gatk"
ref1="path/to/reference1.fasta"
ref2="path/to/reference2.fasta"

# 📁 Unified output base
base_output="path/to/output"
md_dir="${base_output}/marked_duplicates"
vcf_dir="${base_output}/vcf"

# 📁 Reference-specific subfolders
ref1_tag="ref1"
ref2_tag="ref2"

# 📁 Interval files
intervals_ref1="${vcf_dir}/${ref1_tag}.intervals"
intervals_ref2="${vcf_dir}/${ref2_tag}.intervals"

# 🔧 Thread count
threads=8

# 🧬 Sample names
r1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samplesheet" | awk '{print $1}')
r2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samplesheet" | awk '{print $2}')

# 🧬 HaplotypeCaller for ref1
"${gatk_path}/gatk" HaplotypeCaller --java-options "-Xmx60g" \
  -I "${md_dir}/${ref1_tag}/${r2}.bam" \
  -O "${vcf_dir}/${ref1_tag}/${r2}.g.vcf.gz" \
  -ERC GVCF \
  -R "$ref1" \
  -L "$intervals_ref1" \
  --native-pair-hmm-threads "$threads"

# 🧬 HaplotypeCaller for ref2
"${gatk_path}/gatk" HaplotypeCaller --java-options "-Xmx60g" \
  -I "${md_dir}/${ref2_tag}/${r2}.bam" \
  -O "${vcf_dir}/${ref2_tag}/${r2}.g.vcf.gz" \
  -ERC GVCF \
  -R "$ref2" \
  -L "$intervals_ref2" \

```

### 🔹 Step 3: Joint Genotyping and SNP Filtering & Stats

```bash

#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=14-0
#SBATCH --array=1-2%2

# 📁 Base directories
genome_sheet="path/to/2genomes.txt
ref_dir="path/to/references"
gatk_path="path/to/gatk"
vcf_dir="path/to/output/vcf"

# 📁 Reference-specific subfolders
ref1_tag="ref1"
ref2_tag="ref2"

# 🔧 Thread count
threads=4

# 🧬 Reference names
ref_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$genome_sheet" | awk '{print $1}')
ref_tag=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$genome_sheet" | awk '{print $2}')

# 📁 Reference and interval files
ref="${ref_dir}/${ref_name}"
intervals="${vcf_dir}/${ref_tag}.intervals"

# 📁 Output files
combined_gvcf="${vcf_dir}/${ref_tag}/28samples_${ref_tag}.g.vcf.gz"
raw_vcf="${vcf_dir}/${ref_tag}/28samples_${ref_tag}_raw.vcf.gz"
snps_vcf="${vcf_dir}/${ref_tag}/28samples_${ref_tag}_snps.vcf.gz"
filtered_vcf="${vcf_dir}/${ref_tag}/28samples_${ref_tag}_snps_HF.vcf.gz"
biallelic_vcf="${vcf_dir}/${ref_tag}/28samples_${ref_tag}_snps_pass_bial.vcf.gz"
stats_file="${vcf_dir}/${ref_tag}/28samples_${ref_tag}_snps_pass_bial.stats"

echo "🧬 Processing reference ${SLURM_ARRAY_TASK_ID}: ${ref_name} (${ref_tag})"

# 🧬 CombineGVCFs
"${gatk_path}/gatk" CombineGVCFs --java-options "-Xmx60g" \
  -R "$ref" \
  $(for sample in NN114296 NN114297 NN114393 NN115950 NN190240 NNE_114394 \
    SAMN32324407 SAMN32324408 SAMN32324409 SAMN32324410 SAMN32324411 SAMN32324412 \
    SAMN32324413 SAMN32324414 SAMN32324415 SAMN32324416 SAMN32324417 SAMN32324418 \
    SAMN32324419 SAMN32324420 SAMN32324421 SAMN32324422 SAMN32324423 SAMN32324424 \
    SAMN32324425 SAMN32324426 SAMN32324430 SRR13774415; do
      echo "--variant ${vcf_dir}/${ref_tag}/${sample}.g.vcf.gz"
  done) \
  -O "$combined_gvcf" \
  -L "$intervals"

# 🧬 GenotypeGVCFs
"${gatk_path}/gatk" GenotypeGVCFs --java-options "-Xmx60g" \
  -R "$ref" \
  -V "$combined_gvcf" \
  -O "$raw_vcf" \
  -L "$intervals"

# 🧬 Select SNPs
"${gatk_path}/gatk" SelectVariants --java-options "-Xmx60g" \
  -V "$raw_vcf" \
  -select-type SNP \
  -O "$snps_vcf"

# 🧬 Hard Filtering
"${gatk_path}/gatk" VariantFiltration \
  -V "$snps_vcf" \
  -filter "DP < 10.0" --filter-name "DP10" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O "$filtered_vcf"

# 🧬 Extract PASS biallelic SNPs
bcftools view -f PASS -m2 -M2 "$filtered_vcf" -Oz -o "$biallelic_vcf" --threads "$threads"
tabix "$biallelic_vcf"

# 📊 Summary stats
bcftools stats -s - "$biallelic_vcf" > "$stats_file"

```

🔹 Step 4: Annotation with VEP

```
#!/bin/bash

# 📁 Input files
gff_file="path/to/reference.gff"
fasta_file="path/to/reference.fasta"
vcf_input="path/to/output/vcf/ref1/28samples_ref1_snps_pass_bial.vcf"

# 📁 Output files
sorted_gff="path/to/output/ref1/genomic_sort.gff.gz"
vcf_annotated="path/to/output/vcf/ref1/28samples_ref1_snps_pass_bial_vep.vcf"

# 📊 Sort and index GFF
grep -v "#" "$gff_file" | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > "$sorted_gff"
tabix -p gff "$sorted_gff"

# 🧬 Run VEP annotation
vep \
  -i "$vcf_input" \
  --gff "$sorted_gff" \
  --fasta "$fasta_file" \
  --format vcf \
  --everything \
  --vcf \
  --species Neofelis_nebulosa \
  -o "$vcf_annotated"

```

