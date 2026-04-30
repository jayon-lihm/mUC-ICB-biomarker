## oncoprint
## "final_vcf_df" needed from "Filtering_VCF.R"
## Supplemental Figure S5

## Gene list collected from:
## cBioportal: https://www.cbioportal.org/study/summary?id=blca_msk_tcga_2020
## TP53, *TTN, KMT2D, ARID1A, KDM6A, PIK3CA, MUC16, TERT, RB1, KMT2C (>15%)
## TCGA bladder
## TP53, MUC16, KMT2D, KDM6A, ARID1A, PIK3CA, KMT2C, RB1, *FAT4, *CSMD3 (>15%)

source("./bioinf_scripts/src/oncoprint_functions.R")
source("./bioinf_scripts/src/oncoprint_graphic_setting.R")

TMB_df <- read.table("./wes_data/summary_data/TMB_QC_Response.txt",
header=T, as.is=T, sep="\t")

freq_genes <- c("TP53", "TTN", "KMT2D", "ARID1A", "KDM6A", "PIK3CA",
"MUC16", "TERT", "RB1", "KMT2C", "FAT4", "CSMD3")

tier1_samples <- TMB_df$patient_name[TMB_df$QC_pass=="PASS_tier1"]

## oncoprint input
onco_input <- convert_to_oncoprint_mat(final_vcf_df, gene_list = freq_genes,
                         nonsyn_types = nonsyn_include)

colnames(onco_input) <- gsub("_T", "", colnames(onco_input))
colnames(onco_input) <- gsub("15_126_", "Patient", colnames(onco_input))

onco_input <- onco_input[, tier1_samples]

## NEED TO ADD PARTS for BOR columns mono and combo
onco_clinical_input <- data.frame(TMB_df[, c("patient_name", "n_missense", "mono_best_response",
                                                   "combo_best_response")])
rownames(onco_clinical_input) <- onco_clinical_input$patient_name
onco_clinical_input <- onco_clinical_input[ colnames(onco_input),]
onco_clinical_input$mono_best_response <- factor(onco_clinical_input$mono_best_response, levels=c("CR", "PR", "SD", "PD", "NE"))
onco_clinical_input$combo_best_response <- factor(onco_clinical_input$combo_best_response, levels=c("CR", "PR", "SD", "PD", "NE"))

## reorder based on clinical response
reordered_onco_clinical_input <- onco_clinical_input[order(onco_clinical_input$mono_best_response), ]

## Output oncoprint input files
#write.table(onco_input, "./FreqGenes_Oncoprint_15-126_2025_03_12.tsv", 
#            col.names=T, row.names=T, quote=F, sep="\t")
#write.table(reordered_onco_clinical_input, "./TMB_Oncoprint_15-126_2025_03_12.tsv",
#            col.names=T, row.names=T, quote=F, sep="\t")

## reorder genomic input as well
reordered_input <- onco_input[, reordered_onco_clinical_input$patient_name]

pdf("./wes_data/summary_data/Supplemental_Figure_S5.pdf", width=13, height = 8)

oncoPrint(reordered_input,
          alter_fun = alter_fun,
          column_order = reordered_onco_clinical_input$patient_name,
          col = col,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             TMB = anno_barplot(reordered_onco_clinical_input$n_missense),
                                             Mono = reordered_onco_clinical_input$mono_best_response,
                                             Combo = reordered_onco_clinical_input$combo_best_response,
                                             col = list( TMB = "grey",
                                                         Mono =  col_response,
                                                         Combo = col_response)
                                             ),
          column_title = "Freq. mutated genes in BLCA"
        )

dev.off()


