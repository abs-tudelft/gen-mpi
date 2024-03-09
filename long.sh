#!/bin/bash
#SBATCH -N 1
#SBATCH -t 24:0:0
#SBATCH -p thin
#SBATCH --tasks-per-node=1
#SBATCH --exclusive
#SBATCH --output=workflow.out


module load 2021
module load impi/2021.2.0-intel-compilers-2021.2.0

/home/tahmad/tahmad/seqkit split2 --threads=128 -1 /scratch-shared/tahmad/bio_data/long/HG002/HG002.fastq -p 8 -O /scratch-shared/tahmad/bio_data/long/HG002/parts -f

time srun -N 2 -n 2 /home/tahmad/hawk/new/pc/minimap2-chrms/minimap2 -t 128 -ax map-ont -R '@RG\tID:None\tSM:sample\tPL:Pacbio\tLB:sample\tPU:lane' /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/HG002.fastq > /scratch-shared/tahmad/bio_data/long/HG002/ONT/out.sam

time srun -N 1 -n 1 ~/hawk/long/queue_long /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/output/ /home/tahmad/tahmad/tools/ 1
