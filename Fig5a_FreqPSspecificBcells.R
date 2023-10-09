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

df1 <- read.csv("", sep=";")
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

df2$PS1a_freq <- (df2$PS1a/df2$TotalBcells)*100 #add frequencies
df2$PS1b_freq <- (df2$PS1b/df2$TotalBcells)*100
df2$PS2_freq <- (df2$PS2/df2$TotalBcells)*100
df2$PS3_freq <- (df2$PS3/df2$TotalBcells)*100
df2$PS5_freq <- (df2$PS5/df2$TotalBcells)*100
df2$Crossreactive_freq <- (df2$Crossreactive/df2$TotalBcells)*100
df2$empty_freq <- (df2$empty/df2$TotalBcells)*100
df2

#Make a long dataframe:
df3 <- df2 %>% pivot_longer(cols=c("PS1a_freq", "PS1b_freq", "PS2_freq", "PS3_freq", "PS5_freq", "Crossreactive_freq", "empty_freq"),
                            names_to='Serotype',
                            values_to='Serotype_Freq') %>% as.data.frame()

head(df3)
  
#
## 
###
####
##### SA vs EU
####
###
##
#


##################################################################
# Plot together in one single plot - Without FMT and with T-TEST #
##################################################################

df4 <- df3
df4 <- df4[df4$Donor != 'FMT.SA2',] #remove FMT
df4 <- df4[df4$Donor != 'FMT.EU2',] #remove FMT

df4$SA_EU <- factor(df4$SA_EU, levels = c("EU", "SA"))

p4 <- ggplot(df4, aes(x = factor(Serotype, level=c("PS1a_freq", "PS1b_freq", "PS2_freq", "PS3_freq", "PS5_freq",    #reorganize the columns with factor
                                                    "Crossreactive_freq", "empty_freq")), y = log10(Serotype_Freq), fill = SA_EU)) +    #with group=SA_EU you can group boxes again
  geom_point(aes(fill = SA_EU), size = 2, position = position_jitterdodge()) + #needed to separate the points according to the group. "shape = 21" can be added for diff shapes of dots
  geom_boxplot (alpha = 0.7, outlier.shape = NA) +   #remove the outliers shown by boxplot, as you already show it in geom_point
  geom_vline(xintercept=c(5.5), linetype="dashed", size = 0.5) +
  scale_fill_manual(values=c("gray75", "gray20")) +
  scale_x_discrete(labels=c("PS1a_freq" = "PS1a", "PS1b_freq" = "PS1b", "PS2_freq" = "PS2", "PS3_freq" = "PS3",
                            "PS5_freq" = "PS5", "Crossreactive_freq" = "Crossreactive", "empty_freq" = "empty")) +
  labs(y= "Log10 % of total B cells", x = "PS-specific B cells") + 
  stat_compare_means(method = "wilcox.test", aes(group = SA_EU), bracket.size = 0.6, label = "p.signif", size = 4) +   #do T test between SA and EU. with label = "p.format" or "p.signif" you can see p-value/stars
  theme_classic() +
  theme(text = element_text(size = 15)) 
p4

pdf('./GBSFreq_SAvsEU.pdf', width = 7, height = 3)
print(p4)
dev.off() #is required to finish creating the plot.
