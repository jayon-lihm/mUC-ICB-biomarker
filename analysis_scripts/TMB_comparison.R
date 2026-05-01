
## Create dataframe and plots comparing to CM275 tertile lines
## Generates Supplemental Figure S4

source("./src/functions.R") ## to load libraries

TMB_df <- read.table("./wes_data/summary_data/TMB_QC_Response.txt",
header=T, as.is=T, sep="\t")

line_data <- TMB_df %>%
  mutate(QC_class = case_when(
    QC_pass == "FAIL:low_coverage" ~ "FAIL",
    grepl("PASS", QC_pass) ~ "PASS",
    TRUE ~ "."
  )) %>%
  mutate(QC_class2 = factor(QC_class, 
                           levels=c("FAIL", "PASS"))) %>%
  mutate(xstart = as.numeric(QC_class2) - 0.4,
         xend = as.numeric(QC_class2) + 0.4) %>%
  select(xstart, xend) %>%
  distinct() %>%
  arrange(xstart)

## CM275 reported tertiles and median as in 
## Galsky et al, Clinical Cancer Research 2020
## Link: https://aacrjournals.org/clincancerres/article/26/19/5120/82701/Nivolumab-in-Patients-with-Advanced-Platinum
line_data$firstTertile <- 85
line_data$median <- 113
line_data$secondTertile <- 169

line_data <- tidyr::gather(line_data, "location", "values", -xstart, -xend)


## Plotting
g <- TMB_df %>%
  mutate(QC_class = case_when(
    QC_pass == "FAIL:low_coverage" ~ "FAIL",
    grepl("PASS", QC_pass) ~ "PASS",
    TRUE ~ "."
  )) %>%
  mutate(QC_class2 = factor(QC_class, 
                            levels=c("FAIL", "PASS"))) %>%
  ggplot() +
  geom_boxplot(aes(x=QC_class2, y=n_missense)) +
  geom_segment(data=line_data,
               aes(x = xstart, 
                   xend = xend,
                   y = values,
                   yend = values,
                   color=location),
               size = 0.8,
               linetype="dotted") +
  scale_color_manual(
    name = "CM275 statistics",
    labels = c("1st tertile(85)", "median(113)", "2nd tertile(169)"),
    values = c("blue", "green", "red")
  ) +
  theme_bw() +
  xlab("QC class") +
  ylab("number of missense mutations") +
  ggtitle("Number of missense mutations distribution",
          subtitle = "PASS and FAIL classes") +
  theme(legend.box.background = element_rect(colour = "black"))

pdf("./wes_data/summary_data/Supplemental_Figure_S4.pdf")
grid.arrange(g)
dev.off()

