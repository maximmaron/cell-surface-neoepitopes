# Generate QC summary with correct parsing
echo -e "sample\tcohort\tmapping_pct\tjunctions" \
    > /data1/abdelwao/maxim/splicing_pipeline/qc/pass1_qc_summary_v2.tsv

grep -h "Sample:\|Cohort:\|Uniquely mapped reads %\|Number of splices: Total" \
    logs/pass1_*.out 2>/dev/null | \
    awk -F'|' '
        /Sample:/{split($0,a," "); sample=a[2]}
        /Cohort:/{split($0,a," "); cohort=a[2]}
        /Uniquely mapped reads %/{gsub(/ |%/,"",$2); mapping=$2}
        /Number of splices: Total/{
            gsub(/ /,"",$2); junc=$2
            print sample"\t"cohort"\t"mapping"\t"junc
        }
    ' >> /data1/abdelwao/maxim/splicing_pipeline/qc/pass1_qc_summary_v2.tsv

# Flag suspicious samples with corrected data
echo ""
echo "=== Flagged samples ==="
awk -F'\t' 'NR>1 && $3!="" && $4!=""{
    sample=$1; cohort=$2; mapping=$3+0; junc=$4+0
    flag=""
    if(cohort=="beat" && mapping<80) flag=flag"LOW_MAPPING;"
    if(cohort=="tcga" && mapping<65) flag=flag"LOW_MAPPING;"
    if(cohort~/leucegene|pellagatti|madan|maiga/ && mapping<75) flag=flag"LOW_MAPPING;"
    if(cohort=="beat" && junc<100000) flag=flag"LOW_JUNCTIONS;"
    if(cohort=="tcga" && junc<100000) flag=flag"LOW_JUNCTIONS;"
    if(cohort~/leucegene|pellagatti/ && junc<50000) flag=flag"LOW_JUNCTIONS;"
    if(flag!="") print sample"\t"cohort"\t"mapping"%\t"junc"\t"flag
}' /data1/abdelwao/maxim/splicing_pipeline/qc/pass1_qc_summary_v2.tsv | sort -k3 -n

echo ""
echo "=== Summary by cohort ==="
awk -F'\t' 'NR>1 && $3!=""{
    cohort=$2; mapping=$3+0; junc=$4+0
    sum_m[cohort]+=mapping; sum_j[cohort]+=junc; count[cohort]++
    if(mapping<min_m[cohort] || min_m[cohort]=="") min_m[cohort]=mapping
    if(mapping>max_m[cohort]) max_m[cohort]=mapping
} END{
    for(c in count)
        printf "%s\tn=%d\tmean_map=%.1f%%\tmin_map=%.1f%%\tmax_map=%.1f%%\tmean_junc=%d\n",
            c, count[c], sum_m[c]/count[c], min_m[c], max_m[c], sum_j[c]/count[c]
}' /data1/abdelwao/maxim/splicing_pipeline/qc/pass1_qc_summary.tsv | sort