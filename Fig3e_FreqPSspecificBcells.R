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
df3 <- filter(df3, Donor != "XX" | Batch != "1") #Remove batch repeats from 
df3 <- filter(df3, Donor != "XX" | Batch != "3") 
df3 <- filter(df3, Filename != "Batch3_D6rep1_wk0") #remove batch repeat for D6
df3

df3$LongCross <- "Cross-sectional"
df3$LongCross[df3$Donor == "SJ"] <- "Longitudinal"

colnames(df2)
PS <- c("PS1", "PS3", "PS4", "PS5", "PS6A", "PS6B", "PS7F", "PS9V", "PS14", "PS18C",  #make a vector the for loop can go through
        "PS19A", "PS19F", "PS23F", "PS15B", "Crossreactive", "empty")

names(df3)[names(df3) == "Pre/Post(wk9M4)"] <- "PrePostwk9M4"  #change column name to avoid issues
df3$PrePostwk9M4 <- gsub("Pre-PCV13","No PCV13", as.character(df3$PrePostwk9M4))

for(i in unique(PS)){       #go through all PS
  
  df3[, i] <- as.numeric(df3[, i])
  df3$TotalBcells <- as.numeric(df3$TotalBcells) #make sure its numeric
  
  Freq = ((df3[, i]/df3[,"TotalBcells"])*100) #calculate %of Ig+ in all B cells
  
  new_col_name = paste0("Freq_", i) #name the vector of the newly calculated %
  
  df3[,new_col_name] = Freq #add it to the df
  
} 

AllPS <- paste0("Freq_", PS)

df4 <- df3 %>% pivot_longer(cols = c(AllPS),               #long df so each IG can be mapped for each PS
                            names_to='Serotype',
                            values_to='Freqs') %>% as.data.frame()

df4$Serotype <- factor(df4$Serotype, levels=c(AllPS))  #reorder the PS
df4$PrePostwk9M4 <- factor(df4$PrePostwk9M4, levels=c("No PCV13", "Post-PCV13"))  #reorder the groups

p <- ggplot(df4, mapping = aes(x = Serotype, y = Freqs, fill = PrePostwk9M4)) +
  geom_boxplot(aes_string(fill = "PrePostwk9M4"), outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             aes_string(fill = "PrePostwk9M4", shape = "LongCross")) +
  geom_vline(xintercept=c(13.5), linetype="dashed", size = 1) +
  scale_y_log10() +
  ylab("Log 10 % Total B cells") +
  scale_fill_manual(values=c("gray75", "gray20")) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes_string(group = "PrePostwk9M4"), bracket.size = 0.6, label = "p.signif", size = 4) +   #nonpaired wilcons = mann whitney. do ANOVA between SA and EU. with label = "p.format" or "p.signif" you can see p-value/stars
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20)) 
p

pdf(file = ".pdf", width = 15, height = 3.5)   #function to save PDF
p      
dev.off()

#extract P values:
p <- ggplot(df4, mapping = aes(x = Serotype, y = log10(Freqs), fill = PrePostwk9M4)) +
  geom_boxplot(aes_string(fill = "PrePostwk9M4"), outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             aes_string(fill = "PrePostwk9M4", shape = "LongCross")) +
  geom_vline(xintercept=c(13.5), linetype="dashed", size = 1) +
  ylab("Log 10 % Total B cells") +
  scale_fill_manual(values=c("gray75", "gray20")) +
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes_string(group = "PrePostwk9M4"), bracket.size = 0.6, label = "p.format", size = 4) +   #nonpaired wilcons = mann whitney. do ANOVA between SA and EU. with label = "p.format" or "p.signif" you can see p-value/stars
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20)) 
p

PvaluePS <- PS <- c("PS1", "PS3", "PS4", "PS5", "PS6A", "PS6B", "PS7F", "PS9V", "PS14", "PS18C",  
                    "PS19A", "PS19F", "PS23F")
Pvalues <- c(0.0190, 0.7619, 0.0190, 0.0381, 0.1143, 0.0667, 0.0381, 0.0190,
             0.3524, 0.0095, 0.4762, 0.1143, 0.5476)

PvaluesAdjusted <- p.adjust(Pvalues, method='BH') 


