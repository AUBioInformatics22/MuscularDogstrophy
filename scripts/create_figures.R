### Generate plots and accompanying tables for metrics at various points in data processing.
### Metrics include: coverage of raw whole genome, aligned chrX, chrX with duplicates marked,
### SNPs, and INDELs; percent mapped; and percent duplicates.

library(viridis)
library(knitr)

setwd("~/GitHub/MuscularDogstrophy")

# load data from csv (generated from stat_csv.sh)
metrics_data=read.csv(file = "data/csv/all_metrics.csv", header = TRUE,
                      colClasses = c("factor", "numeric", "numeric", "numeric",
                      "numeric", "numeric", "numeric", "numeric", "numeric"))

# define text for coverage legends
legend_text <- c("Raw Whole Genome",
                 "Estimated Raw Chromosome X",
                 "Aligned Chromosome X",
                 "Marked Chromosome X",
                 "Variants Chromosome X")

# define text for percent mapped legend
pm_legend_text <- c("Whole Genome",
                    "X Chromosome")

# define max y values for coverage plots
max_value <- max(metrics_data[c("raw", "est_raw_chrX", "aligned_chrX",
                                "marked_chrX")], na.rm = TRUE)
max_value_all <- max(metrics_data[c("raw", "est_raw_chrX", "aligned_chrX",
                                "marked_chrX", "variant_chrX")], na.rm = TRUE)


# define coordinates for labels depending on how many categories are being graphed
label_coords1 <- barplot(raw ~ id, data = metrics_data, space = 0.5)
label_coords2 <- barplot(cbind(raw, aligned_chrX) ~ id, data = metrics_data,
                         beside = TRUE, space = c(0,0.5))  
label_coords3 <- barplot(cbind(raw, aligned_chrX, marked_chrX) ~ id,
                         data = metrics_data, beside = TRUE, space = c(0,0.5))
label_coords4 <- barplot(cbind(raw, aligned_chrX, marked_chrX, variant_chrX) ~ id,
                         data = metrics_data, beside = TRUE, space = c(0,0.5))

##### TABLES
# generate markdown tables of data to accompany plots
# Note: Concatenates new tables to existing ones.

#tables <- file("analysis/0_figures/data_tables.md", open = "at")
#cat(c("__Coverage__", kable(metrics_data[,c("id", "raw")],
#                            format = "pipe", align = "c",
#                            col.names = c("Sample ID", legend_text[1])),"\n"),
#    file = tables, sep = "\n")

#cat(c("__Coverage__", kable(metrics_data[,c("id", "raw", "aligned_chrX")],
#                            format = "pipe", align = "c",
#                            col.names = c("Sample ID",
#                                          legend_text[c(1,3)])), "\n"),
#    file = tables, sep = "\n")

#cat(c("__Percent Mapped__", kable(metrics_data[,c("id", "percent_mapped_wg", "percent_mapped_chrX")],
#                                  format = "pipe", align = "c",
#                                  col.names = c("Sample ID", pm_legend_text[])), "\n"),
#    file = tables, sep = "\n")

#cat(c("__Percent Duplicates__", kable(metrics_data[,c("id", "percent_duplicates")],
#                                      format = "pipe", align = "c",
#                                      col.names = c("Sample ID", "Percent Duplicates")), "\n"),
#    file = tables, sep = "\n")
#close(tables)

##### COVERAGE PLOTS
# generate coverage plots for each step

## RAW
png(filename="analysis/0_figures/1_coverage.png", width = 700, height = 400)
barplot(
  raw ~ id,
  data = metrics_data,
  space = 0.5,
  beside = TRUE,
  ylim = c(0,max_value*1.3),
  xlab = "Sequence",
  ylab = "Coverage",
  main = "Raw Whole Genome Coverage",
  col = viridis(5)[4])
text(x = label_coords1, y = metrics_data$raw + 1.5, labels = metrics_data$raw)
dev.off()

## RAW VS ALIGNED
png(filename="analysis/0_figures/2_coverage.png", width = 700, height = 400)
barplot(cbind(raw, aligned_chrX) ~ id,
        data = metrics_data,
        space = c(0,0.5),
        beside = TRUE,
        ylim = c(0,max_value*1.3),
        xlab = "Sequence",
        ylab = "Coverage",
        main = "Raw vs. Aligned Coverage",
        legend.text = FALSE,
        #args.legend = ,
        col = viridis(5)[4:3])
legend("topright", inset = c(0,-0.1), legend = legend_text[c(1,3)],
       cex = 0.8, fill = viridis(5)[4:2], bty = "n", xpd = TRUE)
text(x = t(label_coords2), y = cbind(metrics_data$raw,
                                     metrics_data$aligned_chrX) + 1.5,
     labels = cbind(round(metrics_data$raw, digits = 3),
                    round(metrics_data$aligned_chrX,  digits = 3)), cex = 0.8)
dev.off()

## RAW VS ALIGNED VS MARKED
png(filename="analysis/0_figures/3_coverage.png", width = 700, height = 400)
barplot(cbind(raw, aligned_chrX, marked_chrX) ~ id,
        data = metrics_data,
        space = c(0,0.5),
        beside = TRUE,
        ylim = c(0,max_value*1.3),
        xlab = "Sequence",
        ylab = "Coverage",
        main = "Coverage at Different Stages of Processing",
        legend.text = FALSE,
        #args.legend = ,
        col = viridis(5)[4:2])
legend("topright", inset = c(0,-0.1), legend = legend_text[c(1,3,4)], cex = 0.8,
       fill = viridis(5)[4:2], bty = "n", xpd = TRUE)
text(x = t(label_coords3), y = cbind(metrics_data$raw,
                                     metrics_data$aligned_chrX,
                                     metrics_data$marked_chrX) + 1.5,
     labels = cbind(round(metrics_data$raw, digits = 2),
                    round(metrics_data$aligned_chrX,  digits = 2),
                    round(metrics_data$marked_chrX, digits = 2)), cex = 0.8)
dev.off()

## RAW VS ALIGNED VS MARKED VS VARIANTS
png(filename="analysis/0_figures/4_coverage.png", width = 700, height = 400)
barplot(cbind(raw, aligned_chrX, marked_chrX, var_SNP, var_INDEL) ~ id,
        data = metrics_data,
        space = c(0,0.5),
        beside = TRUE,
        ylim = c(0,max_value_all*1.3),
        xlab = "Sequence",
        ylab = "Coverage",
        main = "Coverage at Different Stages of Processing",
        legend.text = FALSE,
        #args.legend = ,
        col = viridis(6)[5:1])
legend("topright", inset = c(0,-0.07), legend = legend_text[c(1,3,4,5,6)],
       cex = 0.8, fill = viridis(6)[5:1], bty = "n", xpd = TRUE)
text(x = t(label_coords5), y = cbind(metrics_data$raw,
                                     metrics_data$aligned_chrX,
                                     metrics_data$marked_chrX,
                                     metrics_data$var_SNP,
                                     metrics_data$var_INDEL) + 2,
     labels = cbind(round(metrics_data$raw, digits = 2),
                    round(metrics_data$aligned_chrX,  digits = 2),
                    round(metrics_data$marked_chrX, digits = 2),
                    round(metrics_data$var_SNP, digits = 2),
                    round(metrics_data$var_INDEL, digits = 2)), cex = 0.7)
dev.off()



##### PERCENT PLOTS
# generate plots of percentage metrics

## PERCENT MAPPED
png(filename="analysis/0_figures/percent_mapped.png", width = 500, height = 300)
barplot(cbind(percent_mapped_wg, percent_mapped_chrX) ~ id,
        data = metrics_data,
        space = c(0,0.5),
        beside = TRUE,
        ylim = c(60,102),
        main = "Percent of Sequences Mapped by Sample",
        xlab = "Sequence",
        ylab = "Percent Mapped",
        xpd = FALSE,
        col = magma(7)[c(5,6)])
legend("topright", inset = c(0,-0.1), legend = pm_legend_text[c(1:2)], cex = 0.8,
       fill = magma(7)[c(5,6)], bty = "n", xpd = TRUE)
text(x = t(label_coords2), y = cbind(metrics_data$percent_mapped_wg,
                                 metrics_data$percent_mapped_chrX) + 2,
     labels = cbind(metrics_data$percent_mapped_wg,
                    metrics_data$percent_mapped_chrX), cex = 0.8)
dev.off()

## PERCENT DUPLICATES
png(filename="analysis/0_figures/percent_duplicates.png", width = 500, height = 300)
barplot(percent_duplicates ~ id,
        data = metrics_data,
        space = 0.5,
        ylim = c(0,20),
        main = "Percent Duplicates by Sample",
        xlab = "Sequence",
        ylab = "Percent Duplicates",
        col = magma(8)[5])
text(x = label_coords1, y = metrics_data$percent_duplicates + 1,
     labels = metrics_data$percent_duplicates, cex = 0.8)
dev.off()
