# 🔹 Mitogenome Assembly and Coverage Analysis

This workflow includes two main steps:

## 🧬 Step 1: Mitogenome Assembly Using GetOrganelle
Trimmed paired-end reads are assembled against a mitochondrial reference to reconstruct circular contigs.

## 📊 Step 2: Coverage Analysis
Assembled contigs are indexed and aligned with BWA, duplicates are removed with Picard, and coverage is calculated using mosdepth.

## 🖥️ SLURM Script

```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G or 250G
#SBATCH --time=14-0
#SBATCH --array=1-26%8

# 📁 Base directories
base_dir="path/to"
samplesheet="${base_dir}/26samples.txt"
reads_dir="${base_dir}/reads"
mito_dir="${base_dir}/mito"
ref_fasta="${mito_dir}/ref/Neofelis_nebulosa_4samples.fasta"
coverage_dir="${mito_dir}/coverage"
bam_dir="${coverage_dir}/bam"
picard_jar="path/to/picard.jar"

# 🔧 Thread count
threads=16

# 🧬 Sample name
r1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samplesheet" | awk '{print $1}')
echo "SLURM Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${r1}"

# 🧬 Step 1: Mitogenome assembly
get_organelle_from_reads.py \
  -1 "${reads_dir}/${r1}_1.trimmed.fastq.gz" \
  -2 "${reads_dir}/${r1}_2.trimmed.fastq.gz" \
  -t $threads \
  -o "${mito_dir}/results/${r1}" \
  -s "$ref_fasta" \
  -F animal_mt

# 🧬 Step 2: Coverage analysis
sample_fasta="${coverage_dir}/${r1}.fasta"
raw_bam="${bam_dir}/${r1}_raw.bam"
dedup_bam="${bam_dir}/${r1}.bam"
metrics_file="${bam_dir}/${r1}_MD.txt"

bwa index "$sample_fasta"
samtools faidx "$sample_fasta"

bwa mem -t $threads -R "@RG\tID:${r1}\tPL:ILLUMINA\tPU:${r1}\tLB:${r1}\tSM:${r1}" \
  "$sample_fasta" \
  "${reads_dir}/${r1}_1.trimmed.fastq.gz" \
  "${reads_dir}/${r1}_2.trimmed.fastq.gz" \
  | samtools view -u - \
  | samtools sort -@ $threads -o "$raw_bam"

samtools index -@ $threads "$raw_bam"

java -Xmx40g -jar "$picard_jar" MarkDuplicates \
  -I "$raw_bam" \
  -O "$dedup_bam" \
  -M "$metrics_file" \
  --CREATE_INDEX true \
  --REMOVE_DUPLICATES true

mosdepth "${r1}" -n -t $threads --fast-mode "$dedup_bam"

```
