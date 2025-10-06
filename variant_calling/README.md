# 🧬 Variant Calling Pipeline

This repository contains a full pipeline for variant calling, filtration, and annotation of genomic data. It is designed to process high-throughput sequencing reads using standard bioinformatics tools.

---

## 📋 Pipeline Overview

The pipeline includes the following steps:

- **Quality Control**: `FastQC`, `MultiQC`
- **Read Trimming**: `fastp`
- **Alignment**: `BWA-MEM2`
- **Post-processing**: `Picard`, `samtools`
- **Variant Calling**: `GATK HaplotypeCaller`

---

## ⚙️ SLURM Batch Script

Below is the SLURM job script used to process samples in parallel using an array job. It performs alignment to two reference genomes, marks duplicates, and collects metrics.

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
qc_raw_dir="${base_output}/qc/raw"
qc_trimmed_dir="${base_output}/qc/trimmed"
fastp_json_dir="${base_output}/fastp/json"
fastp_html_dir="${base_output}/fastp/html"
multiqc_dir="${base_output}/qc/multiqc"
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

# 📊 FastQC для сырых ридов
fastqc "${reads_dir}/${r1}_1.fastq.gz" "${reads_dir}/${r1}_2.fastq.gz" -t 2 -o "$qc_raw_dir"

# ✂️ Обрезка ридов с fastp
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

# 📊 FastQC для обрезанных ридов
fastqc "${reads_dir}/${r1}_1.trimmed.fastq.gz" "${reads_dir}/${r1}_2.trimmed.fastq.gz" -o "$qc_trimmed_dir"

# 📊 MultiQC для объединения всех QC-отчётов
multiqc "$qc_raw_dir" "$qc_trimmed_dir" "$fastp_json_dir" -o "$multiqc_dir"

echo "Processing sample ${SLURM_ARRAY_TASK_ID}: ${r1}, ${r2}"

# 🧷 Alignment for ref1 (по trimmed ридам)
bwa mem -t "$threads" -R "@RG\tID:${r1}\tPL:ILLUMINA\tPU:${r1}\tLB:${r1}\tSM:${r2}" \
    "$ref1" "${reads_dir}/${r1}_1.trimmed.fastq.gz" "${reads_dir}/${r1}_2.trimmed.fastq.gz" \
    | samtools view -u - | samtools sort -@ "$threads" -o "${bam_dir}/${ref1_tag}/${r2}.bam"

samtools index -@ "$threads" "${bam_dir}/${ref1_tag}/${r2}.bam"

# 🧷 Alignment for ref2 (по trimmed ридам)
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

# 📊 Collect metrics
java -jar "$picard_jar" CollectWgsMetrics \
    -I "${bam_dir}/${ref1_tag}/${r2}.bam" \
    -O "${metrics_dir}/${ref1_tag}/${r2}_wgs_metrics.txt" \
    -R "$ref1"

java -jar "$picard_jar" CollectAlignmentSummaryMetrics 
    -I "${bam_dir}/${ref1_tag}/${r2}.bam" 
    -O "${metrics_dir}/${ref1_tag}/${r2}_alignment_metrics.txt" 
    -R "$ref1"

java -jar "$picard_jar" CollectWgsMetrics 
    -I "${bam_dir}/${ref2_tag}/${r2}.bam" 
    -O "${metrics_dir}/${ref2_tag}/${r2}_wgs_metrics.txt" 
    -R "$ref2"

java -jar "$picard_jar" CollectAlignmentSummaryMetrics 
    -I "${bam_dir}/${ref2_tag}/${r2}.bam" 
    -O "${metrics_dir}/${ref2_tag}/${r2}_alignment_metrics.txt" 
    -R "$ref2"
