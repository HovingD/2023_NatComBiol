#Competition ELISA analysis

library(ggplot2)
library(ggrepel)

data <- read.delim('Competition_all.txt')

data$label <- data$PS.1
data$label[data$Condition != 'auto'] <- NA
data$Condition <- gsub('CDAP', 'auto_biotin', data$Condition)
ggplot(data, aes(x=Condition, y=Average, group = PS.1, colour = PS.1)) + 
  geom_point() + 
  geom_line() +
  geom_hline(yintercept = 0, colour = 'red', linetype = 'dashed') +
  theme_bw() + 
  ylab('Percent inhibition') + 
  geom_text_repel(aes(label = label), nudge_x = -0.1, min.segment.length = 10) +
  ggtitle('Competition ELISA')
ggsave('Competition ELISA.pdf', width = 6, height = 5)  
