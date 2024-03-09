## GenMPI
### MPI-based Cluster Scalable Variant Calling forShort/Long Reads Sequencing Data

GenMPI, an MPI based scalable method for both widely used short and long reads aligners, BWA-MEM and Minimap2, respectively. One of the main goals of GenMPI is to ensure 100% identical variant output compared to the single node baseline. In addition, GenMPI provides a flexible architecture which can be used to integrate a variety of alignment and variant calling tools. 

#### Overall architecture:

<img width="1373" alt="genmpi_architecture" src="genmpi_architecture.jpg">

## Usage:
### For short-reads 
```
/home/tahmad/tahmad/seqkit split2 --threads=128 -1 /scratch-shared/tahmad/bio_data/long/HG002/HG002.fastq -p 8 -O /scratch-shared/tahmad/bio_data/long/HG002/parts -f
time srun -N 2 -n 2 /home/tahmad/hawk/new/pc/minimap2-chrms/minimap2 -t 128 -ax map-ont -R '@RG\tID:None\tSM:sample\tPL:Pacbio\tLB:sample\tPU:lane' /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/HG002.fastq > /scratch-shared/tahmad/bio_data/long/HG002/ONT/out.sam
time srun -N 1 -n 1 ~/hawk/long/queue_long /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/long/HG002/ONT/output/ /home/tahmad/tahmad/tools/ 1
```
### For long-reads
