#cd /data1/abdelwao/maxim/splicing_pipeline/leafcutter

#Make sure to change the WORK DIRECTORY in the 

CASE_GROUP="ZRSR2_mut" \
CTRL_GROUP="normal_hematopoietic,normal_tissue" \
WORK_DIR="/data1/abdelwao/shared/splicing_analysis/leafcutter/ZRSR2_vs_normal" \
bash 00_prep_bam_lists.sh \
  --comparison custom \
  --outdir /data1/abdelwao/maxim/splicing_pipeline/leafcutter/ZRSR2_vs_normal/bam_lists

#Use above commands to generate bam lists for each comparison. Then run the following command to convert bams to junction files for LeafCutter input
#sbatch --export=COHORT=srsf2  --array=1-N_CASE  01_bam_to_junc.sh
#sbatch --export=COHORT=normal --array=1-N_CTRL  01_bam_to_junc.sh


sbatch --export=COMPARISON=ZRSR2_vs_normal 02_cluster_introns.sh

COMPARISON=ZRSR2_vs_normal bash 02b_build_groups_file.sh

sbatch --export=COMPARISON=ZRSR2_vs_normal 03_differential_splicing.sh
