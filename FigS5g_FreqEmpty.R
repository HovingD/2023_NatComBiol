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
df2 <- cbind(Donor= rownames(df2), df2)  #Make the index into a new column with Donor numbers
rownames(df2) <- 1:nrow(df2)
df2

df2 <- df2[df2$Donor != 'D43',] #remove donor 43 (is repeat of donor19)
df2 <- df2[df2$Donor != 'D18',] #remove donor 18 (is non-colonized, but colonized 3 weeks prior)

class(df2$PS1a)

df2$TotalBcells <- as.numeric(df2$TotalBcells) #Make sure the vector class is numeric, not character.
df2$PS1a <- as.numeric(df2$PS1a)
df2$PS1b <- as.numeric(df2$PS1b) 
df2$PS2 <- as.numeric(df2$PS2) 
df2$PS3 <- as.numeric(df2$PS3) 
df2$PS5 <- as.numeric(df2$PS5) 
df2$PS5 <- as.numeric(df2$PS5) 
df2$Crossreactive <- as.numeric(df2$Crossreactive) 
df2$empty <- as.numeric(df2$empty) 
  
class(df2$PS1a)

df2$empty_freq <- (df2$empty/df2$TotalBcells)*100
df2

p2 <- ggplot(df2, mapping = aes(x = "empty", y = empty_freq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, alpha = 0.8, size = 3) +
  theme_classic()
p2  

pdf('FrequencyEmpty.pdf', width = 3, height = 5)
print(p2)
dev.off() #is required to finish creating the plot.