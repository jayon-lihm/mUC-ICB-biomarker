## This script is to convert NCBI refseq hg19's exon coordinates
## into a sorted bed file format.

ncbi_refseq_genes <- read.table("./wes_data/reference/NCBI_refseq_genes_HG19.tsv.gz", header=F, as.is=T, sep="\t")
exon_start_end <- ncbi_refseq_genes[, c(3, 10, 11)]
colnames(exon_start_end) <- c("chrom", "exon_starts", "exon_ends")

chrom_list <- c(paste("chr", 1:22, sep=""), "chrX", "chrY")

out <- NULL
for (chr in chrom_list){
    ## process by chromosome
    exon_start_end_chrom <- exon_start_end[which(exon_start_end$chrom == chr),]
    cat(chr, "\n")
    cat(dim(exon_start_end_chrom), "\n")

    ## Origianl table contains comma separated multiple start and end positions (exonStarts and exonEnds)
    ## Separate by comma (unpack) and rbind them into one dataframe
    exon_start_end_unpacked <- apply(exon_start_end_chrom, 1, function(x){
        start <- unlist(strsplit(x[2], ","))
        end <- unlist(strsplit(x[3], ","))
        chrom <- x[1]
        return(data.frame(chrom=chrom, start=start, end=end)) })

    cat(length(exon_start_end_unpacked), "\n")
    ## rbind
    exon_start_end_unpacked2 <- Reduce(function(...) rbind(...), exon_start_end_unpacked)
    ## remove duplicates if any
    exon_start_end_unpacked3 <- exon_start_end_unpacked2[!duplicated(exon_start_end_unpacked2), ]
    ## Re-order by coordinate
    exon_start_end_unpacked3 <- exon_start_end_unpacked3[order(exon_start_end_unpacked3$start, exon_start_end_unpacked3$end),]

    cat(dim(exon_start_end_unpacked3), "\n")
    out <- rbind(out, exon_start_end_unpacked3)
}

## write it
write.table(out, "./NCBI_refseq_genes_exon_HG19_sorted.bed", col.names=F, row.names=F, quote=F, sep="\t")
