## This script is to extract mutations from vcf of regions covered at pre-defined coverage 
## 10x for tumor, 7x for normal

## RUN: bash compute_10x_regions.sh <tumor_bam> <normal_bam> <vcf_file> <outdir>

## samtools and bedtools required
bedtools="path-to-bedtools"
samtools="path-to-samtools"
exon_regions="./wes_data/reference/NCBI_refseq_genes_exon_HG19_sorted.bed" ## output from previous step

###############################################
######### USER DEFINED INPUT/OUTPUT ###########
## DEFINE INPUT bam file paths
tumor_bam=$1
normal_bam=$2
vcf_file=$3

## DEFINE output bed files
outidr=$4
tumor_bam_10xbed="$outdir/tumor_bam_10x.bed"
normal_bam_7xbed="$outdir/normal_bam_7x.bed"
overlap_bed="$outdir/tumor_normal_region_overlap.bed"

## DEFINE VCF files
out_vcf="$outdir/mutations_in_overlap_regions.vcf"
###############################################

## 1. compute depth per position, 
## select depth >= 10 (tumor) >=7 (normal) 
## and intersect with exonic region

## 1) tumor
$samtools depth $tumor_bam | awk '$3>=10' | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2}' | $bedtools intersect -a stdin -b $exon_regions | $bedtools merge -i stdin > $tumor_bam_10xbed

## 2) normal
$samtools depth $normal_bam | awk '$3>=7' | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2}' | $bedtools intersect -a stdin -b $exon_regions | $bedtools merge -i stdin > $normal_bam_7xbed

## find overlapping regions
$bedtools intersect -a $tumor_bam_10xbed -b $normal_bam_7xbed > $overlap_bed

## 2. extract mutations from regions where both criteria meet
pass_target_region=$overlap_bed

cat $vcf_file | grep -v '^#' | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2, $3, $8}' | $bedtools intersect -a stdin -b $pass_target_region -wa > $out_vcf

