library("dplyr")
library("ggplot2")

TMB_df <- read.table("./wes_data/summary_data/TMB_QC_table.txt",
header=T, as.is=T, sep="\t")

####### QC_pass column assignment ######
#TMB_df$QC_pass <- "FAIL:low_coverage"
#TMB_df$QC_pass[TMB_df$PCT_TARGET_BASES_10X.normal > 85 &
#            TMB_df$PCT_TARGET_BASES_10X.tumor > 85 ] <- "PASS_tier1"
#TMB_df$QC_pass[TMB_df$PCT_TARGET_BASES_10X.normal <= 85 &
#            TMB_df$PCT_TARGET_BASES_10X.tumor <= 85 &
#            TMB_df$captured_exon_size > 25*10^6 ] <- "PASS_tier2"
########################### #############   

TMB_df$scaled_TMB_missense <- TMB_df$n_missense / TMB_df$captured_exon_size * 10^6


pdf("./wes_data/summary_data/Supplemental_Figure_S1.pdf")
TMB_df %>%
  ggplot(aes(x=captured_exon_size, y=n_missense, color=QC_pass)) +
  geom_point() +
  ylab("# missense mutations") + xlab("Captured Exon size (bp)") +
  ggtitle("# missense mutations vs % target regions captured at 10X")

TMB_df %>%
  ggplot(aes(x=PCT_TARGET_BASES_10X.tumor, y=captured_exon_size, color=QC_pass)) +
  geom_point() +
  xlab("% Target bases covered >= 10X (tumor)") +
  ylab("Captured Exon size (bp)") +
  ggtitle("Captured exon size vs % target regions captured at 10X")

TMB_df %>%
  ggplot(aes(x=captured_exon_size, y=n_missense)) +
  geom_point(aes(color=QC_pass)) +
  ggtitle("# nonsyn mutations detected in exonic regions")

TMB_df %>%
  ggplot(aes(x=captured_exon_size, y=scaled_TMB_missense)) +
  geom_point(aes(color=QC_pass)) +
  ggtitle("Missense TMB/Mb vs captured exon size")

dev.off()

