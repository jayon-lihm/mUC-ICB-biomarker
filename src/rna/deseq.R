## Author: Estelle Yao (yaon@mskcc.org)

library("DESeq2")
library("edgeR")
library("BiocParallel")
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("ggrepel")

args = commandArgs(trailingOnly=TRUE)
run_name = args[1] ## full, combofull
data_prefix = args[2] ## ./data/rna_seq_data/

## input files: pc_count.csv, sample_info_[run name].csv
## output files: deseq_result_pc_[run name].csv, dds_pc_normalize_[run name].csv

#### protein coding gene
protein_coding = fread(paste(data_prefix, "qc_pc_count_full.tsv.gz", sep="/"))
gene = protein_coding$gene.symbol
rownames(protein_coding) = protein_coding$gene.symbol
protein_coding = protein_coding[,2:ncol(protein_coding)]

#### Clinical information
cr_nr_sample = fread(paste(data_prefix, "/sample_info_", run_name, ".csv", sep=""), header = FALSE)
cr_nr_sample$response = factor(cr_nr_sample$V2, c("Non-responder",'Responder'))

#### Protein coding genes Size factor 
dds_pc <- DESeqDataSetFromMatrix(countData = protein_coding[, .SD, .SDcols =cr_nr_sample$V1],
                                  colData = cr_nr_sample,
                                  design = ~ response)
register(MulticoreParam(10))

dds_pc <- DESeq(dds_pc,betaPrior=TRUE,parallel=TRUE, BPPARAM=MulticoreParam(10))

#### run once 
counts_dds = counts(dds_pc, normalized=TRUE)
#counts_dds$ID = gene

### Results 
res <- results(dds_pc,
               parallel=TRUE, 
               BPPARAM=MulticoreParam(10))
res_df = as.data.frame(res)
res_df$ID = gene

fwrite(as.data.frame(res_df), file=paste0(data_prefix, "/deseq_output/deseq_result_pc_", run_name,".csv"))
fwrite(as.data.frame(counts_dds), file=paste0(data_prefix, "/deseq_output/dds_pc_normalize_", run_name,".csv"))

