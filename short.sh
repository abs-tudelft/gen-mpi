#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=2
#SBATCH -t 8:0:0
#SBATCH -p thin
#SBATCH --output=workflow.out


module load 2021
module load impi/2021.2.0-intel-compilers-2021.2.0

srun -n 1 -N 1 /home/tahmad/tahmad/seqkit split2 --threads=128 -1 /scratch-shared/tahmad/bio_data/long/HG002/HG002.fastq -p 4 -O /scratch-shared/tahmad/bio_data/long/HG002/parts -f &

time srun -n 4 -N 2 /home/tahmad/hawk/new/posix/bwa-chrms/bwa mem -t 64 -R '@RG\tID:sample_lane\tSM:sample\tPL:illumina\tLB:sample\tPU:lane' /scratch-shared/tahmad/bio_data/short/HG002/gatk/Homo_sapiens_assembly38.fasta /scratch-shared/tahmad/bio_data/short/HG002/HG002_1.fastq /scratch-shared/tahmad/bio_data/short/HG002/HG002_2.fastq > /scratch-shared/tahmad/bio_data/short/HG002/HG002.sam

time srun -n 8 -N 2 ~/hawk/short/queue_short /scratch-shared/tahmad/bio_data/short/HG002/gatk/Homo_sapiens_assembly38.fasta /scratch-shared/tahmad/bio_data/short/HG002/output/ /home/tahmad/tahmad/tools/ 25 gatk
