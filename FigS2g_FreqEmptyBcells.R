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
df3

colnames(df2)
PS <- c("PS1", "PS3", "PS4", "PS5", "PS6A", "PS6B", "PS7F", "PS9V", "PS14", "PS18C",  #make a vector the for loop can go through
        "PS19A", "PS19F", "PS23F", "PS15B", "Crossreactive", "empty")

names(df3)[names(df3) == "Pre/Post(wk9M4)"] <- "PrePostwk9M4"  #change column name to avoid issues

for(i in unique(PS)){       #go through all PS
  
  df3[, i] <- as.numeric(df3[, i])
  df3$TotalBcells <- as.numeric(df3$TotalBcells) #make sure its numeric
  
  Freq = ((df3[, i]/df3[,"TotalBcells"])*100) #calculate %of Ig+ in all B cells
  
  new_col_name = paste0("Freq_", i) #name the vector of the newly calculated %
  
  df3[,new_col_name] = Freq #add it to the df
  
} 

df4 <- df3
df4$Timepoint <- factor(df4$Timepoint, levels=c("wk0", "wk3", "wk9", "M4", "M18"))  #reorder the groups
df4$Totalempty <- "empty" 

p2 <- ggplot(df4, mapping = aes(x = Totalempty, y = Freq_empty)) +
  geom_boxplot() +
  geom_jitter(shape = 16, alpha = 0.8, size = 3) +
  theme_classic()
p2  

pdf('FrequencyEmptyCombined.pdf', width = 3, height = 5)
print(p2)
dev.off() #is required to finish creating the plot.
