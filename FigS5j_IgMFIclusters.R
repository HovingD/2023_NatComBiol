####################################
#                                  #
#  PLOT UMAP OF PS-SPECIFIC CELLS  #
#                                  #
####################################

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

#load different files and merge
dfA <- read.csv("GBS_MFIs_mayority.csv", sep=",")
dfA$PS <- as.character(dfA$PS) 
dfA$PS[dfA$PS == '0'] <- 'PSnegBcells'
dfA$PS[dfA$PS == '1'] <- 'Crossreactive'
dfA$PS[dfA$PS == '2'] <- 'PSIa'
dfA$PS[dfA$PS == '3'] <- 'PSIb'
dfA$PS[dfA$PS == '4'] <- 'PSII'
dfA$PS[dfA$PS == '5'] <- 'PSIII'
dfA$PS[dfA$PS == '6'] <- 'PSV'
dfA$PS[dfA$PS == '7'] <- 'empty'

dfB <- read.csv("GBS_MFIs_PS13_41_42.csv", sep=",")
dfB$PS <- as.character(dfB$PS) 
dfB$PS[dfB$PS == '0'] <- 'PSnegBcells'
dfB$PS[dfB$PS == '1'] <- 'Crossreactive'
dfB$PS[dfB$PS == '2'] <- 'PSIa'
dfB$PS[dfB$PS == '3'] <- 'PSII'
dfB$PS[dfB$PS == '4'] <- 'PSIII'
dfB$PS[dfB$PS == '5'] <- 'PSV'

dfC <- read.csv("GBS_MFIs_Donor22.csv", sep=",")
dfC$PS <- as.character(dfC$PS) 
dfC$PS[dfC$PS == '0'] <- 'PSnegBcells'
dfC$PS[dfC$PS == '1'] <- 'Crossreactive'
dfC$PS[dfC$PS == '2'] <- 'PSIa'
dfC$PS[dfC$PS == '3'] <- 'PSIb'
dfC$PS[dfC$PS == '4'] <- 'PSII'
dfC$PS[dfC$PS == '5'] <- 'PSV'
dfB$PS[dfB$PS == '6'] <- 'empty'

dfD <- read.csv("GBS_MFIs_Donor39.csv", sep=",")
dfD$PS <- as.character(dfD$PS) 
dfD$PS[dfD$PS == '0'] <- 'PSnegBcells'
dfD$PS[dfD$PS == '1'] <- 'Crossreactive'
dfD$PS[dfD$PS == '2'] <- 'PSIa'
dfD$PS[dfD$PS == '3'] <- 'PSII'
dfD$PS[dfD$PS == '4'] <- 'PSV'

df0 <- do.call("rbind", list(dfA, dfB, dfC, dfD))

#Check if binding went correct
ncol(dfA)
ncol(dfB)
ncol(dfC)
ncol(dfD)
ncol(df0)

nrow(dfA)+nrow(dfB)+nrow(dfC)+nrow(dfD)   
nrow(df0)

#exclude samples
df0 <- df0[df0$OmiqFileIndex != 'D43specific_cells_allCSV.csv',]
df0 <- df0[df0$OmiqFileIndex != 'D18specific_cells_allCSV.csv',]

#Change names of columns
colnames(df0)[colnames(df0) == "IgD.BUV563.Norm"] ="IgD"
colnames(df0)[colnames(df0) == "IgM.BV570.Norm"] ="IgM"
colnames(df0)[colnames(df0) == "IgG.PE.CF594.Norm"] ="IgG"
colnames(df0)[colnames(df0) == "IgA.APC.Vio770.Norm"] ="IgA"

df0Long <- df0 %>% pivot_longer(cols = c("IgD", "IgM", "IgG", "IgA"),    #long df so each IG can be mapped for each PS
                            names_to='isotype',
                            values_to='MFI') %>% as.data.frame()

df0Long$isotype <- factor(df0Long$isotype, levels = c("IgD", "IgM", "IgG", "IgA"))

p2 <- ggplot(df0Long, aes(x = factor(fsom.metaclust.elbow.1), y = MFI, fill = isotype, color = isotype)) +
  geom_boxplot(outlier.shape = NA, position = "dodge") +
  ylab("MFI") +
  theme_classic() +
  scale_fill_manual(breaks = c("IgD", "IgM", 'IgG', 'IgA'),
                     values=c("black", "grey40", "grey60", "grey85"),
                     guide = guide_legend(override.aes = list(size = 3,
                                                              alpha = 1) )) +
  scale_color_manual(breaks = c("IgD", "IgM", "IgG", "IgA"),
                     values = c("black", "grey40", "grey60", "grey85")) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 20)
  ) +
  ylim(-1.5, 4.5) 
p2

pdf('.pdf', width = 10, height = 1.5)
print(p2)
dev.off()
