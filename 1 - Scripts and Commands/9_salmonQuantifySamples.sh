#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=job5.slurm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kb2784@nyu.edu

module purge
echo ""
echo "moduled purged"

module load salmon/1.0.0
echo ""
echo "salmon model loaded"

salmon quant \
	-i /scratch/kb2784/final.proj/4_salmon_index/Homo_sapiens.GRCh38.cdna.all.normalized.shuffled_index \
	-l A \
	-r /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819990_out.fastq \
	--validateMappings \
	--gcBias \
	--threads ${SLURM_CPUS_PER_TASK} \
	-o SRR7819990.transcripts_quant

salmon quant \
	-i /scratch/kb2784/final.proj/4_salmon_index/Homo_sapiens.GRCh38.cdna.all.normalized.shuffled_index \
	-l A \
	-r /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819991_out.fastq \
	--validateMappings \
	--gcBias \
	--threads ${SLURM_CPUS_PER_TASK} \
	-o SRR7819991.transcripts_quant

salmon quant \
	-i /scratch/kb2784/final.proj/4_salmon_index/Homo_sapiens.GRCh38.cdna.all.normalized.shuffled_index \
	-l A \
	-r /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819992_out.fastq \
	--validateMappings \
	--gcBias \
	--threads ${SLURM_CPUS_PER_TASK} \
	-o SRR7819992.transcripts_quant

salmon quant \
	-i /scratch/kb2784/final.proj/4_salmon_index/Homo_sapiens.GRCh38.cdna.all.normalized.shuffled_index \
	-l A \
	-r /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819993_out.fastq \
	--validateMappings \
	--gcBias \
	--threads ${SLURM_CPUS_PER_TASK} \
	-o SRR7819993.transcripts_quant

salmon quant \
	-i /scratch/kb2784/final.proj/4_salmon_index/Homo_sapiens.GRCh38.cdna.all.normalized.shuffled_index \
	-l A \
	-r /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819994_out.fastq \
	--validateMappings \
	--gcBias \
	--threads ${SLURM_CPUS_PER_TASK} \
	-o SRR7819994.transcripts_quant

salmon quant \
	-i /scratch/kb2784/final.proj/4_salmon_index/Homo_sapiens.GRCh38.cdna.all.normalized.shuffled_index \
	-l A \
	-r /scratch/kb2784/final.proj/1_trim_fastq/1.2_fastp/SRR7819995_out.fastq \
	--validateMappings \
	--gcBias \
	--threads ${SLURM_CPUS_PER_TASK} \
	-o SRR7819995.transcripts_quant

module purge
echo ""
echo "moduled purged"