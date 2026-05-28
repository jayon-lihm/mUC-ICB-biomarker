## This script contains a set of functions used in downstream analyses
library("dplyr")
library("ggplot2")
library("tidyr")
library("gridExtra")


nonsyn_include <- c("splice", "missense_variant", "frameshift_variant", "stop_gained", "start_lost", "stop_lost",
                    "disruptive_inframe_deletion", "disruptive_inframe_insertion")

vcf_header <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")

get_sample_name_from_vcf_file <- function(vcf_file_name, file_pattern){
  return(gsub(file_pattern, "", basename(vcf_file_name)))
}

read_vcf_file <- function(vcf_file){
  if ( endsWith(vcf_file, ".gz") ){
    vcf_df <- read.table( gzfile(vcf_file), header=F, as.is=T, sep="\t")
  } else {
    vcf_df <- read.table(vcf_file, header=F, as.is=T, sep="\t")
  }
  colnames(vcf_df) <- get_vcf_header(vcf_file)
  return(vcf_df)
}

get_vcf_header <- function(vcf_file){
  first_line <- get_first_line(vcf_file)
  header <- strsplit(first_line, "\t")[[1]]
  header[1] <- gsub("#", "", header[1])
  return(header)
}

get_first_line <- function(file_name){
  if ( endsWith(file_name, ".gz") ){
    cmd <- paste("gzcat", file_name, " | grep CHROM")
  } else {
    cmd <- paste("grep CHROM", file_name)
  }
  first_line <- system(cmd, intern = TRUE)
  return(first_line)
}

select_ann_field <- function(ann_field){
  ann_list <- strsplit(gsub("ANN=", "", ann_field), ",")[[1]]
  mutType_list <- sapply(ann_list, "get_mutation_type", USE.NAMES = FALSE)
  if (any(grepl("missense", mutType_list))){
    selected_ann_i <- grep("missense", mutType_list)
  } else {
    selected_ann_i <- 1
  }
  return(ann_list[selected_ann_i])
}

extract_ann_field <- function(info_column){
  info_list <- strsplit(info_column, ";")[[1]]
  ann_field <- info_list[which(startsWith(info_list, "ANN=") |
                               substring(info_list, 1, 2) %in% c("A|", "C|", "G|", "T|") )]
  return(select_ann_field(ann_field) )
}

get_mutation_type <- function(ann){
  return(strsplit(ann, "\\|")[[1]][2])
}

get_gene_name <- function(ann){
  return(strsplit(ann, "\\|")[[1]][4])
}

get_gene_ensembl_ID <- function(ann){
  return(strsplit(ann, "\\|")[[1]][5])
}

get_mutation_region <- function(ann){
  return(strsplit(ann, "\\|")[[1]][6])
}

get_transcript_ensembl_ID <- function(ann){
  return(strsplit(ann, "\\|")[[1]][7])
}

get_amino_acid_change <- function(ann){
  return(strsplit(ann, "\\|")[[1]][11])
}

get_total_cov <- function(x){
  total_cov <- strsplit(x, ":")[[1]][1]
  if ( total_cov == "."){
    total_cov <- 0
  } else {
    total_cov <- as.numeric(total_cov)
  }
  return( total_cov )
}

get_ref_cov <- function(x){
  ref_cov <- strsplit(x, ":")[[1]][2]
  if ( ref_cov == "."){
    ref_cov <- 0
  } else {
    ref_cov <- as.numeric(ref_cov)
  }
  return( ref_cov )
}

parse_vcf_df <- function(vcf_df, sample_name=NULL, vcf_format="Greenbaum"){
  ## parse info field
  x <- sapply(vcf_df$INFO, "extract_ann_field", USE.NAMES = FALSE)
  
  mut_type <- sapply(x, "get_mutation_type", USE.NAMES = FALSE)
  gene_name <- sapply(x, "get_gene_name", USE.NAMES = FALSE)
  gene_id <- sapply(x, "get_gene_ensembl_ID", USE.NAMES = FALSE)
  mut_region <- sapply(x, "get_mutation_region", USE.NAMES = FALSE)
  transcript_id <- sapply(x, "get_transcript_ensembl_ID", USE.NAMES = FALSE)
  AA_change <- sapply(x, "get_amino_acid_change", USE.NAMES = FALSE)
  
  ## parse normal and tumor coverage
  if ( toupper(vcf_format) == "GREENBAUM"){
    normal_total <- sapply(vcf_df$NORMAL, "get_total_cov", USE.NAMES = FALSE)
    normal_ref <- sapply(vcf_df$NORMAL, "get_ref_cov", USE.NAMES = FALSE)
    normal_alt <- normal_total - normal_ref
    tumor_total <- sapply(vcf_df$TUMOR, "get_total_cov", USE.NAMES = FALSE)
    tumor_ref <- sapply(vcf_df$TUMOR, "get_ref_cov", USE.NAMES = FALSE)
    tumor_alt <- tumor_total - tumor_ref
  }
  
  new_df <- cbind(vcf_df[, c("CHROM", "POS", "ID", "REF", "ALT")], 
                  data.frame(mut_type, gene_name, gene_id, mut_region, transcript_id, AA_change, 
                           normal_total, normal_ref, normal_alt, tumor_total, tumor_ref, tumor_alt))
  
  if (!is.null(sample_name)){
    new_df$sample_name <- sample_name
  }
  return(new_df)
}


create_vcf_df <- function(vcf_list, file_pattern){
  all_vcf_df <- NULL
  for (v in vcf_list){
    raw_vcf_df <- read_vcf_file(v)
    vcf_df <- parse_vcf_df(raw_vcf_df, sample_name=get_sample_name_from_vcf_file(v, file_pattern=file_pattern))
    all_vcf_df <- rbind(all_vcf_df, vcf_df)
  }
  
  ## compute VAF
  all_vcf_df$normal_vaf <- (all_vcf_df$normal_total - all_vcf_df$normal_ref) / all_vcf_df$normal_total
  all_vcf_df$tumor_vaf <- (all_vcf_df$tumor_total - all_vcf_df$tumor_ref) / all_vcf_df$tumor_total
  
  return(all_vcf_df)
}

write_vcf_from_vcf_df <- function(vcf_df, file_name_suffix="filtered", output_dir="./",
                                  create_sample_dir=FALSE){
  sample_list <- unique(vcf_df$sample_name)
  select_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "normal_total", "normal_ref", "tumor_total", "tumor_ref")
  
  for (s in sample_list){
    ## If create_sample_dir=TRUE, create sample_name directory and output 
    if (create_sample_dir){
      new_output_dir <- paste(output_dir, "/", s, sep="")
    } else {
      new_output_dir <- output_dir
    }
    
    if (! dir.exists(new_output_dir)){
      dir.create(new_output_dir, recursive=TRUE)
    }
    
    output_file <- paste(new_output_dir, "/", s, "_", file_name_suffix, ".vcf", sep="")
    
    this_sample_vcf_df <- vcf_df[which(vcf_df$sample_name == s), select_columns]
    ## make VCF folumns
    this_sample_vcf_df$QUAL <- "."
    this_sample_vcf_df$FILTER <- "PASS"
    this_sample_vcf_df$INFO <- "INFO"
    this_sample_vcf_df$FORMAT <- "DP:AP"
    this_sample_vcf_df$NORMAL <- paste(this_sample_vcf_df$normal_total, this_sample_vcf_df$normal_ref, sep=":")
    this_sample_vcf_df$TUMOR <- paste(this_sample_vcf_df$tumor_total, this_sample_vcf_df$tumor_ref, sep=":")
    
    write_header(output_file, text_vcf_header)
    write.table(this_sample_vcf_df[, vcf_header_columns], output_file,
                col.names=F, row.names=F, quote=F, sep="\t",
                append = TRUE)
  }
}

