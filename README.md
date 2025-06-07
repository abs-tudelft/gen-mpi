## GenMPI
### MPI-based Cluster Scalable Variant Calling forShort/Long Reads Sequencing Data

GenMPI, an MPI based scalable method for both widely used short and long reads aligners, BWA-MEM and Minimap2, respectively. One of the main goals of GenMPI is to ensure 100% identical variant output compared to the single node baseline. In addition, GenMPI provides a flexible architecture which can be used to integrate a variety of alignment and variant calling tools. 

#### Overall architecture:

<img width="1373" alt="genmpi_architecture" src="genmpi_architecture.jpg">

## Prerequisite:
Compiler
```
MPI runtime (impi/2021.2.0-intel-compilers-2021.2.0)
```
Tools
```
seqkit
samtools
sambamba
vcftools
whatshap
GATK
```
Aligners
```
bwa
minimap2
```
To use bwa-mem aligner, make sure you have above mentioned MPI complier installed: 
```cd bwa; make```

To use minimap2 aligner, make sure you have above mentioned MPI complier installed:
```cd minimap2; make```

All other 

Variant Callers
```
DeepVariant
Clair3
Octopus
```

## Usage:
### For short-reads 
```
#split FASTQ
/home/tahmad/tahmad/seqkit split2 --threads=128 -1 /scratch-shared/tahmad/bio_data/long/HG002/HG002.fastq -p 8 -O /scratch-shared/tahmad/bio_data/long/HG002/parts -f

#run MPI bwa aligner implementation on cluster
time srun -N 2 -n 2 /home/tahmad/hawk/new/pc/minimap2-chrms/minimap2 -t 128 -ax map-ont -R '@RG\tID:None\tSM:sample\tPL:Pacbio\tLB:sample\tPU:lane' /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/HG002.fastq > /scratch-shared/tahmad/bio_data/long/HG002/ONT/out.sam

#run parallel MPI variant calling on cluster
time srun -N 1 -n 1 ~/hawk/long/queue_long /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/output/ /home/tahmad/tahmad/tools/ 1
```
### For long-reads
```
#split FASTQ
/home/tahmad/tahmad/seqkit split2 --threads=128 -1 /scratch-shared/tahmad/bio_data/long/HG002/HG002.fastq -p 8 -O /scratch-shared/tahmad/bio_data/long/HG002/parts -f

#run MPI minimap2 aligner implementation on cluster
time srun -N 2 -n 2 /home/tahmad/hawk/new/pc/minimap2-chrms/minimap2 -t 128 -ax map-ont -R '@RG\tID:None\tSM:sample\tPL:Pacbio\tLB:sample\tPU:lane' /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/HG002.fastq > /scratch-shared/tahmad/bio_data/long/HG002/ONT/out.sam

#run parallel MPI variant calling on cluster
time srun -N 1 -n 1 ~/hawk/long/queue_long /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/output/ /home/tahmad/tahmad/tools/ 1
```
Tanveer Ahmad et al., "GenMPI: Cluster Scalable Variant Calling for Short/Long Reads Sequencing Data", available at: [biorxiv](https://www.biorxiv.org/content/10.1101/2022.04.01.486779v1.full)
