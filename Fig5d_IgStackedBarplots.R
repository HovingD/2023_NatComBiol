###################################################
#                                                 #
#  Stacked Barplots OF PS-SPECIFIC CELLS ISOTYPES #
#                                                 #
###################################################

rm(list=ls()) #cleans the environment

library(ggplot2)
library(plyr)             # install.packages('plyr') always before dplyr, otherwise warning
library(dplyr)            # install.packages('dplyr')
library(tidyr)            # install.packages('tidyr')
library(magrittr)         # install.packages('magrittr')
library(stringr)          # install.packages('stringr')

sink('Script3 session info.txt')
sessionInfo()
sink()

############################
#CROSS-SECTIONAL ANALYSIS: #
############################

df0 <- read.csv(".csv", sep=",")

#analysis:
df1 <- df0 %>% gather('key', 'value', -c(file, Col_NonCol_EU, Colonization, Donor, SA_EU, 
                                         SerotypePS1a, SerotypePS1b, SerotypePS2, SerotypePS3,
                                         SerotypePS5, Serotype_specific_Colonization)) %>%
  mutate(type=str_extract(key, '^PS[[:digit:]]+.'), group=str_extract(key, '\\..+\\.')) %>%
  ddply(~ type + file, mutate, total=sum(value)) %>%
  mutate(total=ifelse(total < 5, NA, total)) %>%
  mutate(percent=100 * (value / total)) %>%
  ddply(~ SA_EU + type + group, mutate, mean_percent=mean(percent, na.rm=TRUE))

#clean-up:
df1$type <- gsub("\\.","",as.character(df1$type))
df1$group <- gsub("\\.","",as.character(df1$group))
df1$type <- gsub("PS1a","PSIa",as.character(df1$type))
df1$type <- gsub("PS1b","PSIb",as.character(df1$type))
df1$type <- gsub("PS2","PSII",as.character(df1$type))
df1$type <- gsub("PS3","PSIII",as.character(df1$type))
df1$type <- gsub("PS5","PSV",as.character(df1$type))
df1$type <- factor(df1$type , levels = c("PSIa","PSIb","PSII","PSIII","PSV"))
df1$SA_EU <- gsub("EU","Dutch",as.character(df1$SA_EU))
df1$SA_EU <- gsub("SA","South African",as.character(df1$SA_EU))
df1$group <- factor(df1$group , levels = c("IgDIgM","IgD", "IgM",  "IgA", "IgG"))
#dfLong1 <- filter(dfLong1, type != "PS6B")  #remove 6B and 23F, not enough cells
#dfLong1 <- filter(dfLong1, type != "PS23F") 

#plot:
p <- ggplot(df1, aes(x=SA_EU, y=mean_percent, fill=as.factor(group))) + 
  geom_bar(stat='identity', position='fill') + 
  facet_wrap(~type, ncol = 6) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Mean Abundance Filtered 5 cells')
p

#Save:
pdf('.pdf', width = 5, height = 4)
print(p)
dev.off() #is required to finish creating the plot.
