

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap) ## for oncoprint
## https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html

nonsyn_include <- c("splice", "missense_variant", "frameshift_variant", "stop_gained", "start_lost", "stop_lost",
                    "disruptive_inframe_deletion", "disruptive_inframe_insertion")


filter_nonsyn_mut_types <- function(vcf_df, nonsyn_types){
  a <- vcf_df %>%
    dplyr::filter( grepl(paste(nonsyn_types, collapse="|"), mut_type ))
  return(a)
}

filter_genes <- function(vcf_df, gene_list, exclude_intergenic=TRUE){
  ## by default: exclude intergenic regions
  
  gene_list2 <- paste("\\<", gene_list, "\\>", sep="")
  if (exclude_intergenic){
    ftd <- vcf_df %>%
      dplyr::filter( !grepl("intergenic", mut_type)) %>%
      dplyr::filter( grepl(paste(gene_list2, collapse="|"), gene_name ))
  } else {
    ftd <- vcf_df %>%
      filter( grepl(paste(gene_list2, collapse="|"), gene_name))
  }
  return(ftd)
}

long_to_wide_by_gene <- function(vcf_df, nonsyn_types=NULL, recoding=TRUE){
  ## column: gene
  ## row: sample
  ## value: mut types (separated by ";")
  
  ## convert format
  gene_simple_mut_df <- collect_muts_by_gene(vcf_df, recoding=recoding)
  
  ## spread it to matrix
  out <- gene_simple_mut_df %>% 
    dplyr::select(gene_name, mut_types, sample_name) %>%
    spread(gene_name, mut_types) %>%
    replace(is.na(.), " ")
  return(out)
}

collect_muts_by_gene <- function(vcf_df, recoding=TRUE){
  ## This function collects muts by gene
  ## so that each sample can have one row per gene
  ## Multiple mutations are collapsed by ";"
  ## Order of selection: 
  ## frameshift > missense > stop_lost
  ## > stop_gain, start_lost 
  ## > disruptive_inframe > splice
  
  if (recoding){
    vcf_df$mut_type <- sapply(vcf_df$mut_type, "recoding_mutation_type", USE.NAMES = FALSE)
  }
  
  out <- vcf_df %>%
    dplyr::select(gene_name, mut_type, sample_name) %>%
    group_by(sample_name, gene_name) %>%
    summarize(mut_types=paste( unique(mut_type), collapse=";"))
  
  return(out)
}

recoding_mutation_type <- function(mutation_type){
  recoded_mutation_type <- NULL
  ## Order of selection: 
  ## frameshift > missense > stop_lost
  ## > stop_gain, start_lost 
  ## > disruptive_inframe > splice
  
  if (grepl("frameshift", mutation_type)){
    recoded_mutation_type <- "frameshift"
  } else if ( grepl("missense", mutation_type)){
    recoded_mutation_type <- "missense"
  } else if ( grepl("stop_lost", mutation_type)){
    recoded_mutation_type <- "stop_lost"
  } else if ( grepl("stop_gain", mutation_type)){
    recoded_mutation_type <- "stop_gain"
  } else {
    recoded_mutation_type <- "other"
  }
  
}

convert_to_oncoprint_mat <- function(vcf_df, gene_list, nonsyn_types=NULL, recoding=TRUE){
  ## [Expected output]
  ## row: gene name
  ## column: sample name
  ## 1) filter for nonsynonymous mutations
  ## 2) select genes
  ## 3) convert to matrix
  ## 4) fill in empty samples (samples without any mutation under given conditions)
  sample_list <- unique(vcf_df$sample_name)
  
  ## filter for mutation types
  if ( ! is.null(nonsyn_types) ){
    selected_vcf_df <- filter_nonsyn_mut_types(vcf_df, nonsyn_types)
  } else {
    selected_vcf_df <- vcf_df
  }
  
  ## filter for genes
  gene_vcf_df <- filter_genes(selected_vcf_df, gene_list)
  
  ## spread to matrix
  gene_mut_mat <- long_to_wide_by_gene(gene_vcf_df, nonsyn_types = nonsyn_types,
                                       recoding=recoding)
  sample_row_names <- gene_mut_mat$sample_name
  gene_mut_mat <- as.matrix(gene_mut_mat[, -1])
  rownames(gene_mut_mat) <- sample_row_names
  
  ## add empty samples
  samples_with_no_muts <- setdiff(sample_list, sample_row_names)
  empty_mat <- matrix(" ", nrow=length(samples_with_no_muts), ncol=ncol(gene_mut_mat))
  rownames(empty_mat) <- samples_with_no_muts
  colnames(empty_mat) <- colnames(gene_mut_mat) ## gene
  gene_mut_mat <- rbind(gene_mut_mat, empty_mat)
  
  return(t(gene_mut_mat))
}

