# Response Adapted Treatment with Nivolumab+Ipilimumab in mUC


Supporting bioinformatics code for Gupta et al (2026, submitted)

```
Response Adapted Immunotherapy with Nivolumab and Salvage Ipilimumab in Metastatic Urothelial Carcinoma: A Single Arm, Phase II Trial with Biomarker Analyses

Aditi Gupta1*, Estelle (Ning) Yao1,2*, Karissa Whiting3*, Niamh M. Keegan1,4, Jayon Lihm3,5,
Tatiana Shcheglova6, Mark J. Bluth7, Hikmat Al-Ahmadie8, Joseph Schmalz9, Liwei Jia10,
Jessica M. Clement11, Suresh Nair11, Adam Hagymasi6, Mohsen Abu-Akeel1, Brooke E. Kania1,
Ashley Regazzi1, Asia S. McCoy1, Etay Ziv1, Colleen Maher1, Jedd D. Wolchok12, Taha
Merghoub13, Gopa Iyer1,4, Dean F. Bajorin1,4, Jonathan E. Rosenberg1,4, Ion Mandoiu3, Yuval
Elhanati3, Katherine S. Panageas3, Benjamin Greenbaum3,14, Margaret K. Callahan4,8*, Samuel
A. Funt1,4#


*All authors contributed equally to this work

#Corresponding Author
```


## Steps

1. Run Picard's CollectHSmetrics to get QC measures

Picard image: `singularity pull docker://levim/picard:2.11`

```
bash ./bioinf_scripts/picard_collectHSMetrics.sh <tumor_bam> <normal_bam> <ref_fasta> <output_dir>
``` 

2. Converts NCBI RefSeq exon coordinates into a sorted bedfile

```
Rscript ./bioinf_scripts/NCBI_refseq_genes_hg19_to_bedfile.R
```

3. Extract regions from tumor and normal regions that are covered at 10X for tumor and 7X for normal

```
bash compute_10x_regions.sh <tumor_bam> <normal_bam> <vcf_file> <outdir>
```

Samtools and Bedtools paths need to be specified within the script.

Bedtools: https://bedtools.readthedocs.io/en/latest/

Samtools: https://www.htslib.org/

4. Filter VCF based on variant allele frequency (VAF)

```
Rscript ./bioinf_scripts/Filtering_VCF.R
```

5. Generates QC Figures in "Supplemental Figure S1".

```
Rscript ./analysis_scripts/WES_QC.R
```


