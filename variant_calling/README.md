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

echo "Processing sample ${SLURM_ARRAY_TASK_ID}: ${r1}, ${r2}"

# 🧷 Alignment for ref1
bwa mem -t "$threads" -R "@RG\tID:${r1}\tPL:ILLUMINA\tPU:${r1}\tLB:${r1}\tSM:${r2}" \
    "$ref1" "${reads_dir}/${r1}_1.fastq.gz" "${reads_dir}/${r1}_2.fastq.gz" \
    | samtools view -u - | samtools sort -@ "$threads" -o "${bam_dir}/${ref1_tag}/${r2}.bam"

samtools index -@ "$threads" "${bam_dir}/${ref1_tag}/${r2}.bam"

# 🧷 Alignment for ref2
bwa mem -t "$threads" -R "@RG\tID:${r1}\tPL:ILLUMINA\tPU:${r1}\tLB:${r1}\tSM:${r2}" \
    "$ref2" "${reads_dir}/${r1}_1.fastq.gz" "${reads_dir}/${r1}_2.fastq.gz" \
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
