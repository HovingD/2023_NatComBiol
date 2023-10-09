###########################################
#                                         #
#  PLOT FREQUENCIES OF PS-SPECIFIC CELLS  #
#                                         #
###########################################

rm(list=ls()) #cleans the environment

#setwd('../')
#setwd("")
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

sink('20230624 Freq comp session info.txt')
sessionInfo()
sink()

df1 <- read.csv(".csv", sep = ",")
df1

df2 <- data.frame(t(df1[-1]))  #transpose columns and rows and make new column names
colnames(df2) <- df1[, 1]
names(df2)[names(df2) == '0.Negative'] <- 'TotalBcells'   #Change the name of 0.Negative, as the dot leads to error
df2 <- cbind(Donor= rownames(df2), df2)  #Make the index into a new column with Donor numbers
rownames(df2) <- 1:nrow(df2)
df2

#make donor number column:
df2$DonorNR <- gsub("Donor(\\d+).*", "Donor\\1", df2$Donor)
df2$Stainmix <- gsub("Donor\\d+(Mix\\d+)", "\\1", df2$Donor)
new_order <- c(1, 9, 10, 2:(ncol(df2)-2))
df2 <- df2[, new_order]
df2

class(df2$PS19F)
df2[, 4:ncol(df2)] <- lapply(df2[, 4:ncol(df2)], as.numeric) #Make sure the vector class is numeric, not character.
class(df2$PS19F)

#make percentages:
for (i in 5:ncol(df2)) {
  col_name <- paste0(names(df2)[i], "_percentage")
  df2[[col_name]] <- df2[[i]] / df2[[4]] * 100
}
df2

#Make a long dataframe, so that you can plot everything together in only one piece of code:
df3 <- df2 %>%
  pivot_longer(cols = contains("percentage"),
               names_to = "serotype",
               values_to = "frequency")  %>% as.data.frame()
head(df3)
df3  

#remove the word "percentage":
df3$serotype <- gsub("_percentage", "", df3$serotype)

#Ensure order of PS
PS <- c("PS7F", "PS19F", "PS15B", "Crossreactive", "empty", "BV711BV785")
df3$serotype <- factor(df3$serotype, levels=c(PS))
Stainmixorder <- c("Mix1", "Mix2", "Mix3", "FMM")
df3$Stainmix <- factor(df3$Stainmix, levels=c(Stainmixorder))
serotypes_to_keep <- c("PS7F", "PS19F", "PS15B")

df4 <- subset(df3, serotype %in% serotypes_to_keep & Donor != "FMM")
df4$Stainmix <- gsub("Mix", "Fluorchrome combination ", df4$Stainmix)

p2 <- ggplot(df4, aes(x = serotype, y = frequency, group = interaction(serotype, DonorNR), color = DonorNR)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_y_log10() +
  geom_point(aes(shape = Stainmix), position = position_dodge(width = 0.8), size = 3) + 
  theme_classic() +
  labs(y = "% total B cells", x = "serotype")
p2

pdf(file = ".pdf", width = 6, height = 3)   #function to save PDF
p2
dev.off()
