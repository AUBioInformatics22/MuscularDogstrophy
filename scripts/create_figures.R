library(viridis)
library(knitr)

setwd("~/GitHub/MuscularDogstrophy")

#v_colors = viridis(4)

metrics_data=read.csv(file = "data/csv/all_metrics.csv", header = TRUE,
                      colClasses = c("factor", "numeric", "numeric", "numeric",
                      "numeric", "numeric", "numeric", "numeric", "numeric"))

legend_text <- c("Raw Whole Genome",
                 "Estimated Raw Chromosome X",
                 "Aligned Chromosome X",
                 "Marked Chromosome X",
                 "Variants Chromosome X")

pm_legend_text <- c("Whole Genome",
                    "X Chromosome")

max_value <- max(metrics_data[c("raw", "est_raw_chrX", "aligned_chrX",
                                "marked_chrX", "variant_chrX")], na.rm = TRUE)



label_coords1 <- barplot(raw ~ id, data = metrics_data, space = 0.5)
label_coords2 <- barplot(cbind(raw, aligned_chrX) ~ id, data = metrics_data,
                         beside = TRUE, space = c(0,0.5))  
label_coords3 <- barplot(cbind(raw, aligned_chrX, marked_chrX) ~ id,
                         data = metrics_data, beside = TRUE, space = c(0,0.5))
label_coords4 <- barplot(cbind(raw, aligned_chrX, marked_chrX, variant_chrX) ~ id,
                         data = metrics_data, beside = TRUE, space = c(0,0.5))

#pd_labels <- barplot(percent_duplicates ~ id, data = metrics_data, ylim = c(0,100))

#label_coord <- barplot(raw ~ id, data = metrics_data, space = 0.5)



#pm_labels <- barplot(cbind(percent_mapped_wg, percent_mapped_chrX) ~ id,
#                     data = metrics_data, space = c(0,0.5), beside = TRUE)

##### TABLES

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
text(x = label_coords1, y = metrics_data$raw + 2, labels = metrics_data$raw)
dev.off()

#png(filename="analysis/0_figures/1_coverage_all.png", width = 500, height = 300)
#barplot(
#  raw ~ rownames(coverage_data),
#  data = coverage_data,
#  space = 0.5,
#  beside = TRUE,
#  ylim = c(0,max_value*1.3),
#  xlab = "Sequence",
#  ylab = "Coverage",
#  main = "Whole Genome Coverage of Raw Sequences",
#  col = v_colors[3])
#dev.off()


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
                                     metrics_data$aligned_chrX) + 2,
     labels = cbind(round(metrics_data$raw, digits = 3),
                    round(metrics_data$aligned_chrX,  digits = 3)), cex = 0.8)
dev.off()

png(filename="analysis/0_figures/3_coverage.png", width = 700, height = 400)
barplot(cbind(raw, aligned_chrX, marked_chrX) ~ id,
        data = metrics_data,
        space = c(0,0.5),
        beside = TRUE,
        ylim = c(0,max_value*1.3),
        xlab = "Sequence",
        ylab = "Coverage",
        main = "Raw vs. Aligned Coverage",
        legend.text = FALSE,
        #args.legend = ,
        col = viridis(5)[4:2])
legend("topright", inset = c(0,-0.1), legend = legend_text[c(1,3,4)], cex = 0.8,
       fill = viridis(5)[4:2], bty = "n", xpd = TRUE)
text(x = t(label_coords3), y = cbind(metrics_data$raw,
                                     metrics_data$aligned_chrX,
                                     metrics_data$marked_chrX) + 2,
     labels = cbind(round(metrics_data$raw, digits = 2),
                    round(metrics_data$aligned_chrX,  digits = 2),
                    round(metrics_data$marked_chrX, digits = 2)), cex = 0.8)
dev.off()

png(filename="analysis/0_figures/4_coverage.png", width = 700, height = 400)
barplot(cbind(raw, aligned_chrX, marked_chrX, variant_chrX) ~ id,
        data = metrics_data,
        space = c(0,0.5),
        beside = TRUE,
        ylim = c(0,max_value*1.3),
        xlab = "Sequence",
        ylab = "Coverage",
        main = "Raw vs. Aligned Coverage",
        legend.text = FALSE,
        #args.legend = ,
        col = viridis(5)[4:1])
legend("topright", inset = c(0,-0.1), legend = legend_text[c(1,3,4,5)],
       cex = 0.8, fill = viridis(5)[4:1], bty = "n", xpd = TRUE)
text(x = t(label_coords4), y = cbind(metrics_data$raw,
                                     metrics_data$aligned_chrX,
                                     metrics_data$marked_chrX,
                                     metrics_data$variant_chrX) + 2,
     labels = cbind(round(metrics_data$raw, digits = 2),
                    round(metrics_data$aligned_chrX,  digits = 2),
                    round(metrics_data$marked_chrX, digits = 2)), cex = 0.8)
dev.off()


##### PERCENT PLOTS

## PERCENT MAPPED
png(filename="analysis/0_figures/percent_mapped.png", width = 500, height = 300)
barplot(cbind(percent_mapped_wg, percent_mapped_chrX) ~ id,
        data = metrics_data,
        space = c(0,0.5),
        beside = TRUE,
        ylim = c(60,102),
        main = "Percent of Sequences Mapped in Samples",
        xlab = "Sequence",
        ylab = "Percent Mapped",
        xpd = FALSE,
        col = magma(8)[c(4,6)])
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
        col = magma(8)[7])
text(x = label_coords1, y = metrics_data$percent_duplicates + 1,
     labels = metrics_data$percent_duplicates, cex = 0.8)
dev.off()



##### HISTOGRAMS

#png(filename="percent_mapped_plot.png", width = 500, height = 300)
# hist(
#   metrics_data$percent_mapped/100,
#   #ylim = c(0,length(mapped_data$percent_mapped)),
#   main = "Percent of Sequences Mapped in Samples",
#   xlab = "Percent Mapped",
#   ylab = "Frequency",
#   col = viridis(5)[3])
# #dev.off()
#
# #png(filename="percent_duplicates_plot.png", width = 500, height = 300)
# #hist(
# #  mapped_data$percent_mapped/100,
# #ylim = c(0,length(mapped_data$percent_mapped)),
# #  main = "Percent of Duplicate  in Samples",
# #  xlab = "Percent Mapped",
# #  ylab = "Frequency",
# #  col = viridis(5)[3])
# #dev.off()

