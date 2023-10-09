###########################################
#                                         #
#  PLOT FREQUENCIES OF PS-SPECIFIC CELLS  #
#                                         #
###########################################

rm(list=ls()) #cleans the environment

setwd('')
getwd()

library(ggplot2)
library(pheatmap)
library(viridisLite)
library(viridis)
library(RColorBrewer)
library(reshape)
library(xlsx)
library(cowplot)
library(grid.arrange)
library(plot_grid)
library(plot_layout)
library(ggarrange)
library(devtools)
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)

sink('session info.txt')
sessionInfo()
sink()

df1 <- read.csv(".csv", sep=";")
df1

df2 <- data.frame(t(df1[-1]))  #transpose columns and rows and make new column names
colnames(df2) <- df1[, 1]
names(df2)[names(df2) == '0.Negative'] <- 'TotalBcells'   #Change the name of 0.Negative, as the dot leads to error
df2 <- cbind(Filename= rownames(df2), df2)  #Make the index into a new column with Donor numbers
rownames(df2) <- 1:nrow(df2)
df2

df3 <- filter(df2, Batch != "2022")  #Delete Batch 2022
df3 <- filter(df3, Timepoint != "wk3") 
df3 <- filter(df3, Timepoint != "M18") 
df3 <- filter(df3, Donor != "XX" | Batch != "1") #Remove batch repeats from XX
df3 <- filter(df3, Donor != "XX" | Batch != "3") 
df3 <- filter(df3, Filename != "Batch3_D6rep1_wk0") #remove batch repeat for D6
df3

df3xreactive <- df3
ps_column_indices <- grep("^PS", colnames(df3xreactive))
columns_to_sum <- c(ps_column_indices, which(colnames(df3xreactive) == "Crossreactive"))
df3xreactive$TotalPScells <- rowSums(df3xreactive[, columns_to_sum], na.rm = TRUE)
df3xreactive$PercNonXreactive <- ((df3xreactive$Crossreactive/df3xreactive$TotalPScells)*100)
 
p3 <- ggplot(df3xreactive, mapping = aes(x = "1", y = PercNonXreactive)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(size = 1) +
  scale_y_log10() +
  ylab("% crossreactive cells") +
  scale_color_manual(values = c("gray75")) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 18))  # Set the y-axis limits from 1 to 14
p3

average_perc <- median(df3xreactive$PercNonXreactive)
average_perc2 <- 100 - average_perc
range_perc <- range(df3xreactive$PercNonXreactive)
range_perc2 <- 100 - range_perc

print(paste("Average PercNonXreactive:", average_perc))
print(paste("Range of PercNonXreactive:", range_perc[1], "-", range_perc[2]))

pdf(file = "FreqCrossreactives.pdf", width = 1.4, height = 2)   #function to save PDF
p3 
dev.off()

