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
Make sure you have installed following tools:
```
seqkit
samtools
sambamba
vcftools
```
Aligners
```
bwa
minimap2
```
To use bwa-mem aligner, make sure you have above mentioned MPI complier installed: 
```
cd bwa; make
```

To use minimap2 aligner, make sure you have above mentioned MPI complier installed:
```
cd minimap2; make
```

All other Variant Callers and phasing tools to be installed as required:
```
DeepVariant
Clair3
Octopus
GATK
Whatshap
```

## Usage:
### For short-reads 
```
#split FASTQ
seqkit split2 --threads=128 -1 HG002.fastq -p 8 -O HG002/parts -f

#run MPI bwa aligner implementation on cluster
time srun -N 2 -n 2 minimap2 -t 128 -ax map-ont -R '@RG\tID:None\tSM:sample\tPL:Pacbio\tLB:sample\tPU:lane' reference/GRCh38_no_alt_analysis_set.fasta ONT/HG002.fastq > ONT/out.sam

#Compile the simple `queue_short.c` script which runs the whole variant calling worlflow
gcc queue_short.c -o queue_short

#run parallel MPI variant calling on cluster
time srun -N 1 -n 1 queue_short reference/GRCh38_no_alt_analysis_set.fasta ONT/output/ tools_path/ 1
```
### For long-reads
```
#split FASTQ
seqkit split2 --threads=128 -1 HG002.fastq -p 8 -O HG002/parts -f

#run MPI minimap2 aligner implementation on cluster
time srun -N 2 -n 2 minimap2 -t 128 -ax map-ont -R '@RG\tID:None\tSM:sample\tPL:Pacbio\tLB:sample\tPU:lane' reference/GRCh38_no_alt_analysis_set.fasta ONT/HG002.fastq > ONT/out.sam

#Compile the simple `queue_long.c` script which runs the whole variant calling worlflow
gcc queue_long.c -o queue_long

#run parallel MPI variant calling on cluster
time srun -N 1 -n 1 queue_long reference/GRCh38_no_alt_analysis_set.fasta NT/output/ tools_path/ 1
```
Tanveer Ahmad et al., "GenMPI: Cluster Scalable Variant Calling for Short/Long Reads Sequencing Data", available at: [biorxiv](https://www.biorxiv.org/content/10.1101/2022.04.01.486779v1.full)
