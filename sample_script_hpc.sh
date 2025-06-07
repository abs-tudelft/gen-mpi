#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=2
#SBATCH -t 8:0:0
#SBATCH -p thin
#SBATCH --output=workflow.out


module load 2021
module load impi/2021.2.0-intel-compilers-2021.2.0

srun -n 1 -N 1 seqkit split2 --threads=128 -1 HG002/HG002_1.fastq -2 HG002/HG002_2.fastq -p 4 -O HG002/parts -f &

time srun -n 4 -N 2 bwa mem -t 64 -R '@RG\tID:sample_lane\tSM:sample\tPL:illumina\tLB:sample\tPU:lane' /gatk/Homo_sapiens_assembly38.fasta HG002/HG002_1.fastq HG002/HG002_2.fastq > HG002/HG002.sam

time srun -n 8 -N 2 queue_short gatk/Homo_sapiens_assembly38.fasta HG002/output/ tools_path/ <threads> gatk
