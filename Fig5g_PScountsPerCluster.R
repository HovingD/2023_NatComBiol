###########################################
#                                         #
#  PLOT FREQUENCIES OF PS-SPECIFIC CELLS  #
#                                         #
###########################################

rm(list=ls()) #cleans the environment

getwd()

library(ggplot2)
library(ggpubr)
library(plyr)             # install.packages('plyr') always before dplyr, otherwise warning
library(dplyr)            # install.packages('dplyr')
library(tidyr)            # install.packages('tidyr')
library(magrittr)         # install.packages('magrittr')
library(stringr)   

sink('Script 4 session info.txt')
sessionInfo()
sink()

#USE data_perc from bootstrap
data_perc2 <- data_perc

#Cleanup
data_perc2$type <- gsub("PS1a","PSIa",as.character(data_perc2$type))
data_perc2$type <- gsub("PS1b","PSIb",as.character(data_perc2$type))
data_perc2$type <- gsub("PS2","PSII",as.character(data_perc2$type))
data_perc2$type <- gsub("PS3","PSIII",as.character(data_perc2$type))
data_perc2$type <- gsub("PS5","PSV",as.character(data_perc2$type))
data_perc2$type <- factor(data_perc2$type, levels = c("PSIa","PSIb","PSII","PSIII","PSV"))
data_perc2$timepoint <- gsub("EU","Dutch",as.character(data_perc2$timepoint))
data_perc2$timepoint  <- gsub("SA","South African",as.character(data_perc2$timepoint))

p <- ggplot(data_perc2, mapping = aes(x = type, y = value, fill = timepoint)) +
  facet_wrap(~cluster, ncol = 1, scales = "free") + 
  geom_boxplot(aes_string(fill = "timepoint"), outlier.shape = NA, alpha = 0.5) +
  geom_point(shape = 16, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             aes_string(fill = "timepoint")) +
  #geom_vline(xintercept=c(13.5), linetype="dashed", size = 0.5) +
  ylab("% Total PS-specific B cells per cluster") +
  scale_fill_manual(values=c("gray75", "gray20")) +
  #stat_compare_means(method = "wilcox.test", aes_string(group = "PrePostPCV13"), bracket.size = 0.6, label = "p.signif", size = 4) +   #do ANOVA between SA and EU. with label = "p.format" or "p.signif" you can see p-value/stars
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20),
        axis.line=element_line()) + scale_y_continuous()

pdf(file = ".pdf", width = 6, height = 60)   #function to save PDF
p     
dev.off()
