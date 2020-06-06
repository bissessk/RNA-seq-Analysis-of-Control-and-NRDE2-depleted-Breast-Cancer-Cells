#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=job4.slurm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kb2784@nyu.edu

module purge
echo ""
echo "moduled purged"

module load salmon/1.0.0
echo ""
echo "salmon model loaded"

salmon index \
	-t /scratch/kb2784/final.proj/3_fasta/Homo_sapiens.GRCh38.cdna.all.normalized.shuffled.fa \
	-i Homo_sapiens.GRCh38.cdna.all.normalized.shuffled_index \
	-k 31
echo ""
echo "indexed"

module purge
echo ""
echo "moduled purged"