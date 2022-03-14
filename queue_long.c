// mpiicc queue.c  utils.c -o queue
// srun -n 8 ./queue $REF  $OUT_DIR  $TOOLS
// srun -n 2 ./queue /scratch-shared/tahmad/bio_data/reference/GRCh38_no_alt_analysis_set.fasta /scratch-shared/tahmad/bio_data/ERR194147/ERR194147/output/ /home/tahmad/tahmad/tools/ 24 

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <assert.h>

//#include "utils.h"

#include <dirent.h>
#include <sys/sysinfo.h>

//CHRMS=['1', '110', '120', '130', '140', '150','160', '170', '180', '190', '2', '210', '220', '230', '240', '250', '260', '270', '280', '290', '3', '31', '32', '33', '34','35','36','37', '4', '41', '42', '43', '44','45','46','47', '5', '51', '52', '53', '54','55','56','57', '6', '61', '62', '63', '64','65','66', '7', '71', '72', '73','74','75','76', '8', '81', '82', '83','84','85', '9', '91', '92', '93','94','95', '10', '101', '102', '103','104','105', '11', '111', '112', '113','114','115', '12', '121', '122', '123','124','125', '13', '131', '132','133','134', '14', '141', '142','143', '15', '151', '152','153', '16', '161', '162','163', '17', '171', '172', '18', '181', '182', '19', '191', '20', '201','202', '21', '211', '22', '221', '23', '231', '232', '233','234','235', '24', '241', '25']

const char *CHRMS[25]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"};

const char *regions[129]={
"chr1:1-24895641","chr1:24895642-49791283","chr1:49791284-74686925","chr1:74686926-99582567","chr1:99582568-124478209","chr1:124478210-149373851","chr1:149373852-174269493","chr1:174269494-199165135","chr1:199165136-224060777","chr1:224060778-248956422",
"chr2:1-24219351","chr2:24219352-48438703","chr2:48438704-72658055","chr2:72658056-96877407","chr2:96877408-121096759","chr2:121096760-145316111","chr2:145316112-169535463","chr2:169535464-193754815","chr2:193754816-217974167","chr2:217974168-242193529",
"chr3:1-24786943","chr3:24786944-49573887","chr3:49573888-74360831","chr3:74360832-99147775","chr3:99147776-123934719","chr3:123934720-148721663","chr3:148721664-173508607","chr3:173508608-198295558",
"chr4:1-23776818","chr4:23776819-47553637","chr4:47553638-71330456","chr4:71330457-95107275","chr4:95107276-118884094","chr4:118884095-142660913","chr4:142660914-166437732","chr4:166437733-190214555",
"chr5:1-22692281","chr5:22692282-45384563","chr5:45384564-68076845","chr5:68076846-90769127","chr5:90769128-113461409","chr5:113461410-136153691","chr5:136153692-158845973","chr5:158845974-181538259",
"chr6:1-24400853","chr6:24400854-48801707","chr6:48801708-73202561","chr6:73202562-97603415","chr6:97603416-122004269","chr6:122004270-146405123","chr6:146405124-170805979",
"chr7:1-22763709","chr7:22763710-45527419","chr7:45527420-68291129","chr7:68291130-91054839","chr7:91054840-113818549","chr7:113818550-136582259","chr7:136582260-159345973",
"chr8:1-24189771","chr8:24189772-48379543","chr8:48379544-72569315","chr8:72569316-96759087","chr8:96759088-120948859","chr8:120948860-145138636",
"chr9:1-23065785","chr9:23065786-46131571","chr9:46131572-69197357","chr9:69197358-92263143","chr9:92263144-115328929","chr9:115328930-138394717",
"chr10:1-22299569","chr10:22299570-44599139","chr10:44599140-66898709","chr10:66898710-89198279","chr10:89198280-111497849","chr10:111497850-133797422",
"chr11:1-22514436","chr11:22514437-45028873","chr11:45028874-67543310","chr11:67543311-90057747","chr11:90057748-112572184","chr11:112572185-135086622",
"chr12:1-22299569","chr12:22299570-44599139","chr12:44599140-66898709","chr12:66898710-89198279","chr12:89198280-111497849","chr12:111497850-133275309",
"chr13:1-22872864","chr13:22872865-45745729","chr13:45745730-68618594","chr13:68618595-91491459","chr13:91491460-114364328", 
"chr14:1-26760928","chr14:26760929-53521857","chr14:53521858-80282786","chr14:80282787-107043718",
"chr15:1-25497796","chr15:25497797-50995593","chr15:50995594-76493390","chr15:76493391-101991189",
"chr16:1-22584585","chr16:22584586-45169171","chr16:45169172-67753757","chr16:67753758-90338345",
"chr17:1-27752479","chr17:27752480-55504959","chr17:55504960-83257441",
"chr18:1-26791094","chr18:26791095-53582189","chr18:53582190-80373285",
"chr19:1-29308807","chr19:29308808-58617616",
"chr20:1-21481388","chr20:21481389-42962777","chr20:42962778-64444167",
"chr21:1-23354990","chr21:23354991-46709983",
"chr22:1-25409233","chr22:25409234-50818468",
"chrX:1-26006814","chrX:26006815-52013629","chrX:52013630-78020444","chrX:78020445-104027259","chrX:104027260-130034074","chrX:130034075-156040895",
"chrY:1-28613706","chrY:28613707-57227415",
"chrM:1-16569"};

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

      char merge_cmd[1024], del_cmd[1024], sort_cmd[1024], index_cmd[1024], vc_cmd[1024], dv_cmd1[1024], index_cmd1[1024], wh_cmd[1024], wh_cmd1[1024], tx_cmd[1024];

      sprintf(merge_cmd, "cat %s*_output_%02d.sam > %s%02d.sam", argv[2], file, argv[2], file);
      fprintf(stderr,"merge_cmd: %s\n", merge_cmd); 
      system(merge_cmd);

      sprintf(del_cmd, "rm  %s*_output_%02d.sam", argv[2], file);
      fprintf(stderr,"del_cmd: %s\n", del_cmd);
      //system(del_cmd);
      
      snprintf (cores, sizeof(cores), "%d",get_nprocs_conf());
      sprintf(sort_cmd, "%ssamtools_install/samtools sort -@ %d -o %s%02d.bam %s%02d.sam", argv[3], atoi(cores)*2, argv[2], file, argv[2], file);    
      fprintf(stderr,"sort_cmd: %s\n", sort_cmd);
      system(sort_cmd);
      
      sprintf(del_cmd, "rm  %s%02d.sam", argv[2], file);
      fprintf(stderr,"del_cmd: %s\n", del_cmd);
      //system(del_cmd);

      sprintf(index_cmd, "%ssamtools_install/samtools index -@ %d %s%02d.bam", argv[3], atoi(cores), argv[2], file);
      fprintf(stderr,"index_cmd: %s\n", index_cmd);
      system(index_cmd);

      sprintf(vc_cmd, "singularity exec /scratch-shared/tahmad/images/clair3_latest.sif /opt/bin/run_clair3.sh --bam_fn=%s%02d.bam --ref_fn=%s --threads=%d --platform=ont --model_path=/opt/models/ont --output=%sclair3/%s --ctg_name=%s", argv[2], file, argv[1], atoi(cores), argv[2], CHRMS[file], CHRMS[file]);
      fprintf(stderr,"vc_cmd: %s\n", vc_cmd);
      system(vc_cmd);

      /*

      sprintf(vc_cmd, "singularity run -B /usr/lib/locale/:/usr/lib/locale/ /scratch-shared/tahmad/images/deepvariant_1.1.0.sif /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref %s --reads %s%02d.bam --output_vcf %s%s.vcf.gz --intermediate_results_dir %s%s_deepvariant --num_shards %d --regions %s", argv[1], argv[2], file, argv[2], CHRMS[file], argv[2], CHRMS[file], atoi(cores), CHRMS[file]);
      //sprintf(vc_cmd, "singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:1.2.0-gpu /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref %s --reads %s%02d.bam --output_vcf %s%s.vcf.gz --intermediate_results_dir %s%s_deepvariant --num_shards %d --regions %s", argv[1], argv[2], file, argv[2], CHRMS[file], argv[2], CHRMS[file], atoi(cores), CHRMS[file]);
      fprintf(stderr,"vc_cmd: %s\n", vc_cmd);
      system(vc_cmd);

      sprintf(wh_cmd, "singularity exec /home/tahmad/tahmad/singularity/dv.simg /home/tahmad/.local/bin/whatshap phase --output %s%s.phased.vcf.gz --reference %s --chromosome %s %s%s.vcf.gz %s%02d.bam", argv[2], CHRMS[file], argv[1], CHRMS[file], argv[2], CHRMS[file], argv[2], file);
      fprintf(stderr,"wh_cmd: %s\n", wh_cmd);
      system(wh_cmd);

      sprintf(tx_cmd, "singularity exec /home/tahmad/tahmad/singularity/dv.simg tabix -p vcf %s%s.phased.vcf.gz", argv[2], CHRMS[file]);
      fprintf(stderr,"tx_cmd: %s\n", tx_cmd);
      system(tx_cmd);
     
      sprintf(wh_cmd1, "singularity exec /home/tahmad/tahmad/singularity/dv.simg /home/tahmad/.local/bin/whatshap haplotag --output %s%s.haplotagged.bam --reference %s %s%s.phased.vcf.gz %s%02d.bam", argv[2], CHRMS[file], argv[1], argv[2], CHRMS[file], argv[2], file);
      fprintf(stderr,"wh_cmd1: %s\n", wh_cmd1);
      system(wh_cmd1);

      sprintf(index_cmd1, "%ssamtools_install/samtools index -@ %d %s%s.haplotagged.bam", argv[3], atoi(cores), argv[2], CHRMS[file]);
      fprintf(stderr,"index_cmd1: %s\n", index_cmd1);
      system(index_cmd1);

      sprintf(dv_cmd1, "singularity run -B /usr/lib/locale/:/usr/lib/locale/ /scratch-shared/tahmad/images/deepvariant_1.1.0.sif /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref %s --reads %s%s.haplotagged.bam --use_hp_information --output_vcf %s%s.wh.vcf.gz --intermediate_results_dir %s%s_deepvariant1 --num_shards %d --regions %s", argv[1], argv[2], CHRMS[file], argv[2], CHRMS[file], argv[2], CHRMS[file], atoi(cores), CHRMS[file]);
      fprintf(stderr,"dv_cmd1: %s\n", dv_cmd1);
      system(dv_cmd1);
      */

    }
  }
 
  MPI_Win_unlock_all(win);
  MPI_Win_free(&win); 

  MPI_Finalize();

  return 0; 
}
