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
   - Extraction of high-quality biallelic SNPs and summary 

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
samplesheet="path/to/samplesheet.txt"
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
samplesheet="path/to/samplesheet.txt"
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


