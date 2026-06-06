bash submit_gatk_vaf_jobs.sh \
      --manifest           /data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv \
      --bed                /data1/abdelwao/maxim/splicing_pipeline/gatk/sf_genes_exons_merged.bed \
      --ref                /data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/GRCh38.primary_assembly.genome.fa \
      --funcotator-sources /data1/abdelwao/maxim/annotations/Homo_sapiens/funcotator_dataSources.v1.8.hg38.20230908s \
      --outdir             /data1/abdelwao/maxim/splicing_pipeline/gatk/output \
      --concurrency 1166

python3 mutation_calling/build_sample_metadata.py \
        --manifest  /data1/abdelwao/maxim/splicing_pipeline/metadata/sample_manifest.tsv \
        --metadata  output/sf_mutation_metadata.tsv \
        --variants  output/sf_all_variants.tsv \
        --outfile   output/sf_final_metadata.tsv

python mutation_calling/call_p95_r102del.py \
    --metadata  output/sf_final_metadata.tsv \
    --ref        /data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/GRCh38.primary_assembly.genome.fa \
    --gtf       /data1/abdelwao/maxim/annotations/Homo_sapiens/GENCODE/primary/gencode.v49.primary_assembly.annotation.bed \
    --outdir    ./del_calls/

python3 mutation_calling/add_p95_r102del.py \
    --metadata   output/sf_final_metadata.tsv \
    --indel-dir  ./del_calls \
    --outfile    sf_final_metadata_with_del.tsv

python make_ggsashimi_bamlist.py \
        --metadata  /data1/abdelwao/maxim/splicing_pipeline/gatk/sf_final_metadata_with_del.tsv \
        --outdir    ./plots 