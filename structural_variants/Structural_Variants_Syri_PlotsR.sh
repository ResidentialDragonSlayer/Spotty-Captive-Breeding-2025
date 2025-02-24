### Structural Variants - Synteny Mapping using Syri and PlotsR ###

# Syri web page and tutorial: https://schneebergerlab.github.io/syri/
# Summary: SyRI is a comprehensive tool for predicting genomic differences between related genomes using whole-genome assemblies (WGA). The assemblies are aligned using whole-genome alignment tools, and these alignments are then used as input to SyRI. SyRI identifies syntenic path (longest set of co-linear regions), structural rearrangements (inversions, translocations, and duplications), local variations (SNPs, indels, CNVs etc) within syntenic and structural rearrangements, and un-aligned regions.

# SyRI uses an unprecedented approach where it starts by identifying longest syntenic path (set of co-linear regions). Since, all non-syntenic regions corresponds to genomic regions which have rearranged between the two genomes, identification of syntenic simultaneously identifies all structural rearrangements as well. After this step, all aligned non-syntenic regions are then classified as either inversion, translocation, or duplication based on the conformation of the constituting alignments. This approach transforms the challenging problem of SR identification to a comparatively easier problem of SR classificaiton.

# Further, SyRI also identifies local variations within all syntenic and structurally rearranged regions. Local variations consists of short variations like SNPs, and small indels as well as structural variations like large indels, CNVs (copy-number variations), and HDRs. Short variations are parsed out from the constituting alignments, where as structural variations are predicting by comparing the overlaps and gaps between consecutive alignments of a syntenic or rearranged region.

### Steps:
## 1) Prepare whole genome files
		# 1a) Use DGenies (https://dgenies.toulouse.inra.fr) to visualize potential chromosomal scaffolds between files. Scaffolds that are reversed will need to be reverse complemented in Geneious.
		# 1b) Use Geneious Prime to reverse complement scaffolds if needed, and rename so the scaffolds match in both files.
## 2) Perform whole genome alignment with minimap2
## 3) Run Syri 
## 4) Plot genomic structures predicted by SyRI using PlotsR

## 1) Prepare Whole Genome Files: 
# Ideally, syri expects that the homologous chromosomes in the two genomes would have exactly same chromosome id. Therefore, it is recommended that the user pre-processes the fasta files to ensure that homologous chromosomes have exactly the same id in both fasta files corresponding to the two genomes. In case, that is not the case, syri would try to find homologous genomes using whole genome alignments, but that method is heuristical and can result in suboptimal results. Also, it is recommended that the two genomes (fasta files) should have same number of chromosomes.

# Genomes were run against each other using D-Genies. We took note of the aligned scaffolds using the interactive D-Genies dot plot, and renamed scaffolds to match in geneious.

## 2) Whole genome alignment using Minimap2:
# Code was taken from Syri "working example" page on their website (https://schneebergerlab.github.io/syri/pipeline.html).

##Installation
mamba create -n syri -c conda-forge -c bioconda cython numpy scipy pandas=0.23.4 biopython psutil matplotlib=3.0.0
mamba install -c conda-forge -c bioconda cython numpy
mamba install -c conda-forge -c bioconda python=3.5
mamba install -c conda-forge -c bioconda scipy pandas=0.23.4 
mamba install -c conda-forge -c bioconda biopython psutil matplotlib=3.0.0
mamba install -c conda-forge python-igraph
mamba install -c bioconda pysam
mamba install -c bioconda longestrunsubsequence

python3 setup.py install		            # Install syri
chmod +x syri/bin/syri syri/bin/chroder	syri/bin/plotsr	# Make files executable

### Slurm job for Syri:

#!/bin/bash
#SBATCH -J Syri_CL2 # job name
#SBATCH -o log_slurm.o%j # output and error file name (%j expands to jobID)
#SBATCH -n 28 # total number of tasks requested
#SBATCH -N 1 # number of nodes you want to run on
#SBATCH -p bsudfq
#SBATCH -t 12:00:00 # run time (hh:mm:ss) - 12.0 hours in this example.
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#
#Activate environment
. ~/.bashrc
conda activate syri
#
syri -c mapped_clouded_leopard.sam -r /bsuscratch/morgancalahan/SMSC_Clouded_Leopard/new/structural_variants/SNNU_for_svs.fasta -q /bsuscratch/morgancalahan/SMSC_Clouded_Leopard/new/structural_variants/mNeoNeb_for_svs.fasta -F S
#
#Syri Explanation
#-r: reference genome
#-q: query genome
#-F S: indicates that input file alignment is in SAM format

#resulting files:
-rw-rw-r-- 1 morgancalahan morgancalahan 6.4K Jun 30 13:39 invOut.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 1.1K Jun 30 13:39 TLOut.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 1.3K Jun 30 13:39 invTLOut.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 1.2K Jun 30 13:39 dupOut.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 2.0K Jun 30 13:39 invDupOut.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 1.9K Jun 30 13:39 ctxOut.txt
-rw-rw-r-- 1 morgancalahan morgancalahan  97K Jun 30 13:39 synOut.txt
-rw-rw-r-- 1 morgancalahan morgancalahan  19M Jun 30 13:39 sv.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 5.7K Jun 30 13:39 notAligned.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 807M Jun 30 13:40 snps.txt
-rw-rw-r-- 1 morgancalahan morgancalahan 180M Jun 30 13:41 syri.out
-rw-rw-r-- 1 morgancalahan morgancalahan 250M Jun 30 13:41 syri.vcf
-rw-rw-r-- 1 morgancalahan morgancalahan  584 Jun 30 13:41 syri.summary
-rw-rw-r-- 1 morgancalahan morgancalahan 8.5K Jun 30 13:41 syri.log


### Slurm job for PlotsR:
#
#!/bin/bash
#SBATCH -J Syri_Plot_Struct # job name
#SBATCH -o log_slurm.o%j # output and error file name (%j expands to jobID)
#SBATCH -n 28 # total number of tasks requested
#SBATCH -N 1 # number of nodes you want to run on
#SBATCH -p bsudfq
#SBATCH -t 12:00:00 # run time (hh:mm:ss) - 12.0 hours in this example.
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=morgancalahan@boisestate.edu
#
#Activate environment
. ~/.bashrc
conda activate syri
#
plotsr --sr /bsuscratch/morgancalahan/SMSC_Clouded_Leopard/new/structural_variants/syri.out --genomes /bsuscratch/morgancalahan/SMSC_Clouded_Leopard/new/structural_variants/genomes.txt -H 8 -W 5 -o /bsuscratch/morgancalahan/SMSC_Clouded_Leopard/new/structural_variants/plotsr.pdf
#
# Plotsr Explanation
#--sr: specifies input file (generated by syri)
#-H 8: specifies height of plot to be 8 inches
#-W 5: specifies width of plot to be 5 inches
#-o : output
