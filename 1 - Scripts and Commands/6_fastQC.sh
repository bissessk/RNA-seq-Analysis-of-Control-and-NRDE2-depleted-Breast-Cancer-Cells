#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=job2.slurm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kb2784@nyu.edu

# purge modules
module purge
echo "purged any already loaded models"
echo ""
# load modules
module load fastqc/0.11.8
echo "loaded fastqc"
echo ""


fastqc /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819995_out.fastq \
	--outdir=/scratch/kb2784/final.proj/2_multiqc
echo "fastqc SRR7819995 done"
echo ""
fastqc /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819994_out.fastq \
	--outdir=/scratch/kb2784/final.proj/2_multiqc
echo "fastqc SRR7819994 done"
echo ""
fastqc /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819993_out.fastq \
	--outdir=/scratch/kb2784/final.proj/2_multiqc
echo "fastqc SRR7819993 done"
echo ""
fastqc /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819992_out.fastq \
	--outdir=/scratch/kb2784/final.proj/2_multiqc
echo "fastqc SRR7819992 done"
echo ""
fastqc /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819991_out.fastq \
	--outdir=/scratch/kb2784/final.proj/2_multiqc
echo "fastqc SRR7819991 done"
echo ""
fastqc /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819990_out.fastq \
	--outdir=/scratch/kb2784/final.proj/2_multiqc
echo "fastqc SRR7819990 done"
echo ""

# purge modules
module purge
echo "purged fastqc"
echo ""