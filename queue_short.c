// mpiicc queue.c  utils.c -o queue
// srun -n 8 ./queue $REF  $OUT_DIR  $TOOLS
// srun -n 2 ./queue /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/ERR194147/ERR194147/output/ /home/tahmad/tahmad/tools/ 24 

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <assert.h>
#include <dirent.h>
#include <sys/sysinfo.h>

const char *CHRMS[25]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"};

int nprocs, rank;

int main(int argc, char *argv[]) {

  MPI_File fh;
  int err;
  MPI_Status status;

  int NUM_FILES=atoi(argv[4]); 

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  fprintf(stderr,"Hello world from processor %s, rank %d, file count: %d \n", processor_name, rank, NUM_FILES);

  // create the window 
  int *baseptr;
  MPI_Info info;
  MPI_Win win;
  MPI_Info_create(&info);
  MPI_Info_set(info, "osc_rdma_acc_single_intrinsic", "true");
  MPI_Win_allocate(sizeof(int), 1, MPI_INFO_NULL, MPI_COMM_WORLD, &baseptr, &win);
  *baseptr = 0;
  int assert;
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Win_lock_all(0,win);

  // allocate and process files while available 
  const int one = 1;
  int file = 0;
  char cores[4]="";
  char findex[3]="";
  while (file < NUM_FILES) {
    MPI_Fetch_and_op(&one, &file, MPI_INT, 0, 0, MPI_SUM, win);
    MPI_Win_flush(0, win);

    if (file < NUM_FILES) {
      //printf("Process:%d -- File: %s\n", rank, namelist[file+2]->d_name);//process_file(file);

      fprintf(stderr,"Hello world from processor %s, rank %d, File number: %02d, Path: %s \n", processor_name, rank, file, argv[2]);

      char merge_cmd[1024], del_cmd[1024], sort_cmd[1024], index_cmd[1024], md_cmd[1024], vc_cmd[1024], bqsr_cmd[1024], applybqsr_cmd[1024];

      sprintf(merge_cmd, "cat %s*_output_%02d.sam > %s%02d.sam", argv[2], file, argv[2], file);
      fprintf(stderr,"merge_cmd: %s\n", merge_cmd); 
      system(merge_cmd);

      sprintf(del_cmd, "rm  %s*_output_%02d.sam", argv[2], file);
      fprintf(stderr,"del_cmd: %s\n", del_cmd);
      //system(del_cmd);
      
      snprintf (cores, sizeof(cores), "%d",get_nprocs_conf());
      sprintf(sort_cmd, "%ssamtools_install/samtools sort -@ %d -o %s%02d_sorted.bam %s%02d.sam", argv[3], atoi(cores)*2, argv[2], file, argv[2], file);    
      fprintf(stderr,"sort_cmd: %s\n", sort_cmd);
      system(sort_cmd);
      
      sprintf(del_cmd, "rm  %s%02d.sam", argv[2], file);
      fprintf(stderr,"del_cmd: %s\n", del_cmd);
      //system(del_cmd);

      sprintf(md_cmd, "%ssambamba markdup -t %d %s%02d_sorted.bam %s%02d.bam", argv[3], atoi(cores), argv[2], file, argv[2], file);
      fprintf(stderr,"md_cmd: %s\n", md_cmd);
      system(md_cmd);

      sprintf(index_cmd, "%ssamtools_install/samtools index -@ %d %s%02d.bam", argv[3], atoi(cores), argv[2], file);
      fprintf(stderr,"index_cmd: %s\n", index_cmd);
      system(index_cmd);

      if(argv[5]=="octopus"){
      
        sprintf(vc_cmd, "singularity exec /scratch-shared/tahmad/images/octopus_latest.sif octopus --threads %d  --reference %s  --reads %s%02d.bam -T %s --sequence-error-model PCRF.NovaSeq --forest /opt/octopus/resources/forests/germline.v0.7.4.forest -o %s%s.vcf.gz", atoi(cores), argv[1], argv[2], file, CHRMS[file], argv[2], CHRMS[file]);
        fprintf(stderr,"vc_cmd: %s\n", vc_cmd);
        system(vc_cmd);
      
      }

      else if(argv[5]=="deepvariant"){

        sprintf(vc_cmd, "singularity run -B /usr/lib/locale/:/usr/lib/locale/ /scratch-shared/tahmad/images/deepvariant_1.1.0.sif /opt/deepvariant/bin/run_deepvariant --model_type WGS --ref %s --reads %s%02d.bam --output_vcf %s%s.dv.vcf.gz --intermediate_results_dir %s%s_deepvariant --num_shards %d --regions %s", argv[1], argv[2], file, argv[2], CHRMS[file], argv[2], CHRMS[file], atoi(cores), CHRMS[file]);
        //sprintf(vc_cmd, "singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:1.2.0-gpu /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref %s --reads %s%02d.bam --output_vcf %s%s.vcf.gz --intermediate_results_dir %s%s_deepvariant --num_shards %d --regions %s", argv[1], argv[2], file, argv[2], CHRMS[file], argv[2], CHRMS[file], atoi(cores), CHRMS[file]);
        fprintf(stderr,"vc_cmd: %s\n", vc_cmd);
        system(vc_cmd);
      
      }

      else if(argv[5]=="gatk"){

        sprintf(bqsr_cmd, "java -jar %sgatk.jar BaseRecalibrator -R /scratch-shared/tahmad/bio_data/short/HG002/gatk/Homo_sapiens_assembly38.fasta -L %s  -I %s%02d.bam --known-sites  /scratch-shared/tahmad/bio_data/short/HG002/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites /scratch-shared/tahmad/bio_data/short/HG002/gatk/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /scratch-shared/tahmad/bio_data/short/HG002/gatk/Homo_sapiens_assembly38.known_indels.vcf.gz  -O %s%s.table",argv[3], CHRMS[file], argv[2], file, argv[2], CHRMS[file]);
        fprintf(stderr,"bqsr_cmd: %s\n", bqsr_cmd);
        system(bqsr_cmd);

        sprintf(applybqsr_cmd, "java -jar %sgatk.jar ApplyBQSR -R /scratch-shared/tahmad/bio_data/short/HG002/gatk/Homo_sapiens_assembly38.fasta  -L %s -I %s%02d.bam -bqsr %s%s.table -O %s%s.bam", argv[3], CHRMS[file], argv[2], file, argv[2], CHRMS[file], argv[2], CHRMS[file]);
        fprintf(stderr,"applybqsr_cmd: %s\n", applybqsr_cmd);
        system(applybqsr_cmd);

        sprintf(vc_cmd, "java -jar %sgatk.jar HaplotypeCaller -R /scratch-shared/tahmad/bio_data/short/HG002/gatk/Homo_sapiens_assembly38.fasta -L %s -I %s%s.bam -O %s%s.vcf.gz", argv[3], CHRMS[file], argv[2], CHRMS[file], argv[2], CHRMS[file]);
        fprintf(stderr,"vc_cmd: %s\n", vc_cmd);
        system(vc_cmd);
      
      }

      free(merge_cmd);
      free(del_cmd);
      free(sort_cmd);
      free(index_cmd);
      free(md_cmd);
      free(vc_cmd);
      free(bqsr_cmd);
      free(applybqsr_cmd);

    }
  }
 
  MPI_Win_unlock_all(win);
  MPI_Win_free(&win); 

  MPI_Finalize();

  return 0; 
}
