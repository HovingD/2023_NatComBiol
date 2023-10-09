###########################################
#                                         #
#  PLOT FREQUENCIES OF PS-SPECIFIC CELLS  #
#                                         #
###########################################

rm(list=ls()) #cleans the environment

setwd('')
getwd()

library(ggplot2)
library(dplyr)
library(tidyr)
library(janitor)

install.packages("janitor")

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
library(stringr)

sink('session info.txt')
sessionInfo()
sink()

df1Bcells <- read.csv(".csv", sep=",")
df2Bcells <- data.frame(t(df1Bcells[-1]))  #transpose columns and rows and make new column names
colnames(df2Bcells) <- df1Bcells[, 1]
names(df2Bcells)[names(df2Bcells) == '0.Negative'] <- 'Totalcells'   #Change the name of 0.Negative, as the dot leads to error
df2Bcells <- cbind(Filename= rownames(df2Bcells), df2Bcells)  #Make the index into a new column with Donor numbers
rownames(df2Bcells) <- 1:nrow(df2Bcells)
df2Bcells$Batch <- c(3,10,1)
df2Bcells$Celltype <- "Bcells"
df2Bcells

df1TNKcells <- read.csv(".csv", sep=",")
df2TNKcells <- data.frame(t(df1TNKcells[-1]))  #transpose columns and rows and make new column names
colnames(df2TNKcells) <- df1TNKcells[, 1]
names(df2TNKcells)[names(df2TNKcells) == '0.Negative'] <- 'Totalcells'   #Change the name of 0.Negative, as the dot leads to error
df2TNKcells <- cbind(Filename= rownames(df2TNKcells), df2TNKcells)  #Make the index into a new column with Donor numbers
rownames(df2TNKcells) <- 1:nrow(df2TNKcells)
df2TNKcells$Batch <- c(3,10,1)
df2TNKcells$Celltype <- "TNKcells"
df2TNKcells

df1Monocytes <- read.csv(".csv", sep=",")
df2Monocytes <- data.frame(t(df1Monocytes[-1]))  #transpose columns and rows and make new column names
colnames(df2Monocytes) <- df1Monocytes[, 1]
names(df2Monocytes)[names(df2Monocytes) == '0.Negative'] <- 'Totalcells'   #Change the name of 0.Negative, as the dot leads to error
df2Monocytes <- cbind(Filename= rownames(df2Monocytes), df2Monocytes)  #Make the index into a new column with Donor numbers
rownames(df2Monocytes) <- 1:nrow(df2Monocytes)
df2Monocytes$Batch <- c(3,10,1)
df2Monocytes$Celltype <- "Monocytes"
df2Monocytes

compare_df_cols(df2Bcells, df2TNKcells, df2Monocytes)  #add columns that have 0 cells
df2TNKcells$PS18C <- 0
df2TNKcells$PS5 <- 0
df2TNKcells$PS9V <- 0
df2Monocytes$PS5 <- 0  

df3 <- rbind(df2Bcells, df2TNKcells, df2Monocytes)

PS <- c("PS1", "PS3", "PS4", "PS5", "PS6A", "PS6B", "PS7F", "PS9V", "PS14", "PS18C", 
        "PS19A", "PS19F", "PS23F", "PS15B", "Crossreactive", "empty")

for(e in unique(PS)){       #go through all PS
  
  df3[, e] <- as.numeric(df3[, e])
  df3$Totalcells <- as.numeric(df3$Totalcells) #make sure its numeric
  
  Freq = ((df3[, e]/df3[,"Totalcells"])*100) #calculate %of Ig+ in all B cells
  
  new_col_name = paste0("Freq_", e) #name the vector of the newly calculated %
  
  df3[,new_col_name] = Freq #add it to the df
  
} 

AllPS <- paste0("Freq_", PS)

df4 <- df3 %>% pivot_longer(cols = c(AllPS),               #long df so each IG can be mapped for each PS
                            names_to='Serotype',
                            values_to='Freqs') %>% as.data.frame()

df4$Serotype <- factor(df4$Serotype, levels=c(AllPS)) #apply a factor to reorder the X-axis

p1 <- ggplot(df4, mapping = aes_string(x = "Serotype", y = "Freqs"), fill = "Celltype") +
  geom_boxplot(aes_string(fill = "Celltype"), outlier.shape = NA, alpha = 0.5) +
  geom_point(shape = 16, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             aes_string(fill = "Celltype")) +
  geom_vline(xintercept=c(13.5), linetype="dashed", size = 1) +
  ylab("% Total cells") +
  scale_fill_manual(values=c("gray90", "gray60", "gray20")) +
  stat_compare_means(method = "anova", aes_string(group = "Celltype"), bracket.size = 0.6, label = "p.signif", size = 4) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
p1

p2 <- ggplot(df4, mapping = aes(x = Serotype, y = log10(Freqs)), fill = Celltype) +
  geom_boxplot(aes_string(fill = "Celltype"), outlier.shape = NA, alpha = 0.5) +
  geom_point(shape = 16, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             aes_string(fill = "Celltype")) +
  geom_vline(xintercept=c(13.5), linetype="dashed", size = 1) +
  ylab("Log 10 % Total cells") +
  scale_fill_manual(values=c("gray90", "gray60", "gray20")) +
  stat_compare_means(method = "anova", aes_string(group = "Celltype"), bracket.size = 0.6, label = "p.signif", size = 4) +   #do ANOVA between SA and EU. with label = "p.format" or "p.signif" you can see p-value/stars
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
p2

pdf(file = "PSspecificBTNKcells.pdf", width = 15, height = 4)   #function to save PDF
p1
p2  
dev.off()
