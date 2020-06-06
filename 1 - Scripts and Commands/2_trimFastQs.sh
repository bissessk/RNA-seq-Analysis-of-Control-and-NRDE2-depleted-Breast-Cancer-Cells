#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=job1.slurm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kb2784@nyu.edu

# purge modules
module purge
echo "purged any already loaded models"
echo ""
# load modules
module load fastp/20190829
echo "fastp is loaded"
echo ""

# SRR7819995
echo ""
echo "Trying to trim SRR7819995.fastq "
echo ""
fastp \
	-i /scratch/kb2784/final.proj/1_trim_fastq/1.1_gzipped_fastqs/SRR7819995.fastq \
	-o SRR7819995_out.fastq \
	--length_required 75 \
	-g \
	--html SRR7819995.fastp.html \
	--json SRR7819995.fastp.json

echo ""
echo "trimmed SRR7819995.fastq "
echo ""

# SRR7819994
echo ""
echo "Trying to trim SRR7819994.fastq "
echo ""
fastp \
	-i /scratch/kb2784/final.proj/1_trim_fastq/1.1_gzipped_fastqs/SRR7819994.fastq \
    -o SRR7819994_out.fastq \
	--length_required 75 \
	-g \
	--html SRR7819994.fastp.html \
	--json SRR7819994.fastp.json      
echo ""
echo "trimmed SRR7819994.fastq "  
echo ""
# SRR7819993
echo ""
echo "Trying to trim SRR7819993.fastq "
echo ""
fastp \
	-i /scratch/kb2784/final.proj/1_trim_fastq/1.1_gzipped_fastqs/SRR7819993.fastq \
    -o SRR7819993_out.fastq \
	--length_required 75 \
	-g \
	--html SRR7819993.fastp.html \
	--json SRR7819993.fastp.json        

echo ""
echo "trimmed SRR7819993.fastq "
echo ""

# SRR7819992
echo ""
echo "Trying to trim SRR7819992.fastq "
echo ""
fastp \
	-i /scratch/kb2784/final.proj/1_trim_fastq/1.1_gzipped_fastqs/SRR7819992.fastq \
    -o SRR7819992_out.fastq \
	--length_required 75 \
	-g \
	--html SRR7819992.fastp.html \
	--json SRR7819992.fastp.json
echo ""
echo "trimmed SRR7819992.fastq "
echo ""
# SRR7819991
echo ""
echo "Trying to trim SRR7819991.fastq "
echo ""
fastp \
	-i /scratch/kb2784/final.proj/1_trim_fastq/1.1_gzipped_fastqs/SRR7819991.fastq \
    -o SRR7819991_out.fastq \
	--length_required 75 \
	-g \
	--html SRR7819991.fastp.html \
	--json SRR7819991.fastp.json

echo ""
echo "trimmed SRR7819991.fastq "
echo ""

# SRR7819990
echo ""
echo "Trying to trim SRR7819990.fastq "
echo ""
fastp \
	-i /scratch/kb2784/final.proj/1_trim_fastq/1.1_gzipped_fastqs/SRR7819990.fastq \
    -o SRR7819990_out.fastq \
	--length_required 75 \
	-g \
	--html SRR7819990.fastp.html \
	--json SRR7819990.fastp.json
echo ""
echo "trimmed SRR7819990.fastq "
echo ""
# purge modules
module purge