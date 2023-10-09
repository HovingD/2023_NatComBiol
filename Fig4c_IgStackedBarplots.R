###########################################
#                                         #
#  Stacked Barplots OF PS-SPECIFIC CELLS  #
#                                         #
###########################################

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

#########################
#LONGITUDINAL ANALYSIS: #
########################

dfLong0 <- read.csv("./.csv", sep=",")

#analysis:
dfLong1 <- dfLong0 %>% gather('key', 'value', -c(file, donor, PrePostPCV13, timepoint)) %>%
  mutate(type=str_extract(key, '^PS[[:digit:]]+.'), group=str_extract(key, '\\..+\\.')) %>%
  ddply(~ type + file, mutate, total=sum(value)) %>%
  mutate(total=ifelse(total < 5, NA, total)) %>%
  mutate(percent=100 * (value / total)) %>%
  ddply(~ timepoint + type + group, mutate, mean_percent=mean(percent, na.rm=TRUE))

#clean-up:
dfLong1$type <- gsub("\\.","",as.character(dfLong1$type))
dfLong1$group <- gsub("\\.","",as.character(dfLong1$group))

dfLong1$type <- factor(dfLong1$type , levels = c("PS1","PS3","PS4","PS5","PS6A","PS6B","PS7F","PS9V","PS14",
                                            "PS18C","PS19A","PS19F", "PS23F","PS15B"))
dfLong1$timepoint <- gsub("pre","wk0",as.character(dfLong1$timepoint))
dfLong1$timepoint <- factor(dfLong1$timepoint , levels = c("wk0", "wk3", "wk9", "M18"))
dfLong1$group <- factor(dfLong1$group , levels = c("IgDIgM","IgD", "IgM",  "IgA", "IgG"))
dfLong1 <- filter(dfLong1, type != "PS6B")  #remove 6B and 23F, not enough cells
dfLong1 <- filter(dfLong1, type != "PS23F") 

#plot:
pLong <- ggplot(dfLong1, aes(x=timepoint, y=mean_percent, fill=as.factor(group))) + 
  geom_bar(stat='identity', position='fill') + 
  facet_wrap(~type, ncol = 6) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Mean Abundance Filtered 5 cells')
pLong

#Save:
pdf('./.pdf', width = 7, height = 4)
print(pLong)
dev.off() 