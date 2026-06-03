# Code repository for Gupta et al (2026, unver review)

### Response Adapted Treatment with Nivolumab+Ipilimumab in mUC

This repository provides codes for analyses and figures included in Gupta et al (2026, under review).


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

## Download

Download "data" folder from:

https://mskcc.box.com/s/v35jtbe6tkeswz66x8z9tafwp2g1yij3

## WES Steps

WES analysis is performed by Jayon Lihm (lihmj@mskcc.org).

### Data processing

1. Run Picard's CollectHSmetrics to get QC measures

Picard image: `singularity pull docker://levim/picard:2.11`

```
bash ./src/picard_collectHSMetrics.sh <tumor_bam> <normal_bam> <ref_fasta> <output_dir>
``` 

2. Converts NCBI RefSeq exon coordinates into a sorted bedfile

```
Rscript ./src/NCBI_refseq_genes_hg19_to_bedfile.R
```

3. Extract regions from tumor and normal regions that are covered at 10X for tumor and 7X for normal

```
bash ./bionf_scripts/compute_10x_regions.sh <tumor_bam> <normal_bam> <vcf_file> <outdir>
```

Samtools and Bedtools paths need to be specified within the script.

Bedtools: https://bedtools.readthedocs.io/en/latest/

Samtools: https://www.htslib.org/

4. Filter VCF based on variant allele frequency (VAF)

```
Rscript ./src/Filtering_VCF.R
```

### Supplemental Figure S1

Generates QC Figures in "Supplemental Figure S1".

```
Rscript ./analysis_scripts/WES_QC.R
```

### Supplemental Figure S4

Generates TMB comparison to CheckMate 275 TMB: "Supplemental Figure S4"

```
Rscript ./analysis_scripts/TMB_comparison.R
```

### Supplemental Figure S5

Generates Oncoprint
```
Rscript ./analysis_scripts/Oncoprint.R
```

## RNA-seq Steps

RNA-seq analysis is performed by Ning (Estelle) Yao. (yaon@mskcc.org)

Conda environment for the analysis can be set up with `setup.sh` script.

```
conda activate mUC_ICB_RNA

```

### Figure 2A + Supp Figure S6

`./analysis_scripts/Gene_Expression_Analysis.ipynb` 

DESEQ outputs need to be fed into the Jupyter notebook. The deseq output files used for generating figures are in `./data/rna_seq_data/deseq_output`.


