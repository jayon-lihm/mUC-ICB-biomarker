## This script is to run picard's "collectHsMetrics"
## using Agilent SureSelect Human All Exon V4 coordinate (S03723314)
## to compute the median depth 

## bash picard_collectHSMetrics.sh <tumor_bam> <normal_bam> <ref_fasta> <output_dir>

target_interval="./wes_data/reference/S03723314_Covered.bed.2.interval_list"

## picard version: 2.11

###############################################
######### USER DEFINED INPUT/OUTPUT ###########
tumor_bam=$1
normal_bam=$2
ref_fasta=$3 #"path to b37(hg19) reference fasta"
output_dir=$4
################################################


## Tumor
java -Xmx{params.mem}g -jar /picard-2.11/picard.jar CollectHsMetrics \
             I=$normal_bam  \
             O=$output_dir/normal_picard_hs_metrics.txt \
             R=$ref_fasta \
             BAIT_INTERVALS=$target_interval \
             TARGET_INTERVALS=$target_interval \
             PER_TARGET_COVERAGE=$output_dir/normal_picard_per_target_coverage.txt

## Normal
java -Xmx{params.mem}g -jar /picard-2.11/picard.jar CollectHsMetrics \
             I=$tumor_bam  \
             O=$output_dir/tumor_picard_hs_metrics.txt \
             R=$ref_fasta \
             BAIT_INTERVALS=$target_interval \
             TARGET_INTERVALS=$target_interval \
             PER_TARGET_COVERAGE=$output_dir/tumor_picard_per_target_coverage.txt

