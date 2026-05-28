### This script contains workflow to read  VCF files and filter based on 

source("./src/functions.R")
filter_conditions <- c("t_tot_thres"=10, "t_alt_thres"=0, "t_vaf_thres"=0.02,
                       "n_tot_thres"=10, "n_vaf_thres"=0.01)

###### USER DEFINED PATH #########
unftd_vcf_list <- ""
##################################

unftd_vcf_path <- "../server_data/annotations"
unftd_vcf_list <- Sys.glob( paste(unftd_vcf_path, "/*/snpeff/*_T_ann.vcf", sep="") )
unftd_all_vcf_df <- create_vcf_df(unftd_vcf_list, file_pattern = "_ann.vcf")
unftd_all_vcf_df$sample_name <- gsub("-", "_", unftd_all_vcf_df$sample_name)

final_vcf_df <- unftd_all_vcf_df %>%
  mutate(normal_alt = normal_total - normal_ref,
         tumor_alt = tumor_total - tumor_ref) %>%
  filter(
    normal_total >= filter_conditions["n_tot_thres"] & 
      normal_vaf < filter_conditions["n_vaf_thres"] &
      tumor_total >= filter_conditions["t_tot_thres"] & 
      tumor_alt >= filter_conditions["t_alt_thres"] & 
      tumor_vaf >= filter_conditions["t_vaf_thres"] )

### OPTIONAL STEP: Output final_vcf_df into VCF files
#write_vcf_from_vcf_df(final_vcf_df, file_name_suffix = "filtered_new", output_dir = "./data/new_filtered_vcf",
#                      create_sample_dir=TRUE)

### Summary numbers used for QC
## count number of missense mutations
TMB_missense_df <- final_vcf_df %>%
  filter(grepl("missense", mut_type)) %>%
  group_by(sample_name) %>%
  summarize(n_missense=n())

## count number of mutations (any type)
TMB_all_df <- final_vcf_df %>%
  group_by(sample_name) %>%
  summarize(n_all_muts=n())

## count number of nonsynonymous mutations
## nonsyn types defined in "nonsyn_include"
TMB_nonsyn_df <- final_vcf_df %>%
    filter(grepl(paste(nonsyn_include, collapse="|"), mut_type)) %>%
    group_by(sample_name) %>%
    summarize(n_nonsyn_muts=n())

df_list <- list(TMB_all_df, TMB_nonsyn_df, TMB_missense_df)
TMB_merged_df <- Reduce(function(x, y) merge(x, y, by = "sample_name", all = TRUE), df_list)

TMB_merged_df[is.na(TMB_merged_df)] <- 0

