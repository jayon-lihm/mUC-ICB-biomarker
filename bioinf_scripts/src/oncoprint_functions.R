

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap) ## for oncoprint
## https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html


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

