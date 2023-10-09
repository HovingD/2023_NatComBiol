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

df0 <- read.csv(".csv", sep=",")
df0

df1 <- df0
df1[is.na(df1)] <- 0

#transpose columns and rows and make new column names:
df2 <- data.frame(t(df1[-1]))  
colnames(df2) <- df1[, 1]

#Change the name of 0.Negative, as the dot leads to error
names(df2)[names(df2) == '0.Negative'] <- 'TotalCells'  

#Make the index into a new column with Donor numbers
df2 <- cbind(Donor= rownames(df2), df2)  
rownames(df2) <- 1:nrow(df2)
df2

#make new columns for group and celltype:
df2$group <- do.call(c, lapply(strsplit(df2$Donor, c('_')), function(i)i[1]))
df2$celltype <- do.call(c, lapply(strsplit(df2$Donor, c('_')), function(i)i[2]))

#remove numbers from group
df2$group <- gsub('[[:digit:]]+', '', df2$group)

#change names:
df2$celltype[df2$celltype  == "T.cells"] <- "Tcells"

#make numeric
class(df2$PS1a)
df2$TotalCells <- as.numeric(df2$TotalCells) #Make sure the vector class is numeric, not character.
df2$PS1a <- as.numeric(df2$PS1a)
df2$PS1b <- as.numeric(df2$PS1b) 
df2$PS2 <- as.numeric(df2$PS2) 
df2$PS3 <- as.numeric(df2$PS3) 
df2$PS5 <- as.numeric(df2$PS5) 
df2$PS5 <- as.numeric(df2$PS5) 
df2$Crossreactive <- as.numeric(df2$Crossreactive) 
#df2$empty <- as.numeric(df2$empty) 
class(df2$PS1a)

df2$PS1a_freq <- (df2$PS1a/df2$TotalCells)*100 #add frequencies
df2$PS1b_freq <- (df2$PS1b/df2$TotalCells)*100
df2$PS2_freq <- (df2$PS2/df2$TotalCells)*100
df2$PS3_freq <- (df2$PS3/df2$TotalCells)*100
df2$PS5_freq <- (df2$PS5/df2$TotalCells)*100
df2$Crossreactive_freq <- (df2$Crossreactive/df2$TotalCells)*100
#df2$empty_freq <- (df2$empty/df2$TotalBcells)*100
df2

#Make a long dataframe, so that you can plot everything together in only one piece of code:
df3 <- df2 %>% pivot_longer(cols=c("PS1a_freq", "PS1b_freq", "PS2_freq", "PS3_freq", "PS5_freq", "Crossreactive_freq"),
                            names_to='Serotype',
                            values_to='Serotype_Freq') %>% as.data.frame()
head(df3)

#remove "_freq"
df3$Serotype <- gsub('_freq', '', df3$Serotype)

#Change PS Arabic for Roman:
df3$Serotype[df3$Serotype == "PS1a"] <- "PSIa"
df3$Serotype[df3$Serotype == "PS1b"] <- "PSIb"
df3$Serotype[df3$Serotype == "PS2"] <- "PSII"
df3$Serotype[df3$Serotype == "PS3"] <- "PSIII"
df3$Serotype[df3$Serotype == "PS5"] <- "PSV"
  
#make new group:
df3$group2 <- paste(df3$celltype,df3$group)

#reorder:
df3$Serotype <- factor(df3$Serotype, levels = c("PSIa", "PSIb", "PSII", "PSIII", "PSV", "Crossreactive"))
df3$celltype <- factor(df3$celltype, levels = c("Bcells", "Monocytes", "Tcells"))
df3$group2 <- factor(df3$group2, levels = c("Bcells NoBlock","Bcells SiglecBlock", "Monocytes NoBlock", 
                                            "Monocytes SiglecBlock","Tcells NoBlock", "Tcells SiglecBlock"))

#plot:
plot1 <- ggplot(df3, aes(x = Serotype, y = Serotype_Freq), group = group2, color = factor(celltype)) +
  geom_boxplot(aes(fill = group2), alpha = 0.7, outlier.shape = NA) +
  geom_point (aes(fill = group2, shape=factor(celltype)), size=2, position = position_jitterdodge()) +
  #stat_compare_means(method = "T.test", aes(group = group2), bracket.size = 20, label = "p.signif", size = 4) +
  theme_classic() +
  geom_vline(xintercept=c(5.5), linetype="dashed", size = 0.5) +
  scale_fill_manual(values=c("gray75", "gray20","gray75", "gray20","gray75", "gray20")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  #rotate x-axis labels
  labs(y= "% per total celltype", colour = "group", shape = "celltype")
plot1

pdf('FreqBTMcells2.pdf', width = 12, height = 3)
print(plot1)
dev.off() #is required to finish creating the plot.



