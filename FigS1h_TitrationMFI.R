
rm(list=ls()) #cleans the environment

library(ggplot2)

setwd("")

df1 <- read.csv(".csv", sep=";")

df1$Conc <- as.character(df1$Conc)
class(df1$Conc)

df1$Temp <- factor(df1$Temp, levels = c("RT", "Ice"))

p <- ggplot(df1, aes(x=Conc, y=MFI, fill=Temp)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("MFI PS7F-SA-BV785") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 11000)) +
  scale_fill_manual(values = c("#d3d3d3","#545454")) +
  ggtitle('PBMC MFI PS7F-SA-BV785')
p

pdf('.pdf', width = 9, height = 3)
print(p)
dev.off() 
