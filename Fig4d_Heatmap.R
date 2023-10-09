#############
# HEATMAP   #
#############

#Information for PS gates:
# 0 = B cells
# 1 = Cross-reactive
# 2 = PS1
# 3 = PS14
# 4 = PS15B
# 5 = PS18C
# 6 = PS19A
# 7 = PS19F
# 8 = PS23F
# 9 = PS3
# 10 = PS4
# 11 = PS5
# 12 = PS6A
# 13 = PS6B
# 14 = PS7F
# 15 = PS9V
# 16 = empty

#clean environment
rm(list=ls())

setwd('R:\\')

library("CATALYST")
library("RColorBrewer")
library("grid")
source("plotMultiHeatmap_customized.R")
library('pheatmap')
library('ComplexHeatmap')
library(reshape2)
library(stringi)

sink('session info.txt')
sessionInfo()
sink()

###
## Data
#

# meta 
meta <- read.csv2("metadata.csv", na.strings=c("","NA"))
meta$Filename <- gsub(pattern = ".csv", replacement = ".fcs", x = meta$Filename)
meta$Donor <- factor(trimws(meta$Donor), levels = unique(meta$Donor))
meta$Pre.Post.PCV13 <- factor(meta$Pre.Post.PCV13, levels = unique(meta$Pre.Post.PCV13))
meta$Pre.Post.wk9M4. <- factor(meta$Pre.Post.wk9M4., levels = unique(meta$Pre.Post.wk9M4.))
meta$Timepoint <- factor(meta$Timepoint, levels = unique(meta$Timepoint))
meta$Control <- factor(meta$Control, levels = unique(meta$Control))
meta$Batch <- factor(meta$Batch, levels = unique(meta$Batch))
meta$Filename  <- factor(meta$Filename, 
                         levels = meta$Filename[order(meta$Donor)])


# panel
panel <- read.csv2("panel.csv", na.strings=c("","NA"))

# fcs files
fcsFiles <- list.files(path = getwd(), pattern = ".fcs", full.names = TRUE)

# construct SingleCellExperiment
sce <- CATALYST::prepData(x = fcsFiles, md = meta, panel = panel, features = panel$channel_name,
                          panel_cols = list(channel = "channel_name",
                                            antigen = "marker_label",
                                            class = "marker_class"),
                          
                          md_cols = list(file = "Filename", id = "Donor", 
                                         factors = c("Pre.Post.PCV13","Pre.Post.wk9M4.", "Timepoint", "Control", 
                                                     "Batch")),
                          transform = FALSE,
                          FACS = TRUE)
names(sce@assays) <- "exprs"

#don't do this here as then it will not create the heatmaps for the phenotype using these cells!!!!!!
#if you want to remove them from the heatmap, it's best to remove them from the heatmap at the end
#sce <- sce[,!sce@assays@data$exprs["PS",] %in% c(0, 1, 7)] #use this to remove 0, 1 and 7


### 
## Import Set clustering as metadata
#

clusters <- assay(sce, i = "exprs")

##get expression of markers
markers <- data.frame(t(unlist(clusters[c(25,1:23),]))) #exclude PS (24) and fsom-metaclust-elbow-1 (25)
colnames(markers)[1] <- 'fs'
expression <- split(markers, f = markers$fs)
quant <- lapply(expression, function(x){
    a <- x
    a <- a[,2:ncol(a)]
    #apply(a, 2, quantile, .5)
    apply(a, 2, median)
})
median_expr <- do.call(rbind, quant)

cluster <- clusters["fsom-metaclust-elbow-1",]
tst <- factor(unlist(clusters["fsom-metaclust-elbow-1",]))
sce$cluster_id <- tst

tst2 <- factor(unlist(clusters["PS",]))
sce$cell_id <- tst2

#make the split per PS type
split_data <- split(sce@colData, sce@colData$cell_id)

#get counts from the non-PS
negcounts <- table(split_data[[1]]$cluster_id)

#get the cluster frequencies of PS negative B cells:
negfreqs <- negcounts/sum(negcounts)*100

#make empty matrices
abundance2 <-  abundance1 <-  abundance <- 
    matrix(ncol=length(split_data), nrow=length(unique(cluster)))

#go per PS type to get the abundance out
for(i in 1:length(split_data)){
    counts <- table(split_data[[i]]$cluster_id, split_data[[i]]$sample_id)
    
    #remove empty donors
    counts <- counts[,colSums(counts)>0]
    
    #make mean abundance
    props <- apply(counts, 2, FUN = function(x) x/sum(x)*100)
    props[!is.finite(props)] <- NA
    colSums(props, na.rm = T)
    abundance[,i] <- apply(props, 1, FUN=mean, na.rm=T)
    
    #make total abundance and store in abundance 1
    total <- rowSums(counts)
    abundance1[,i] <- total / sum(total) * 100
    
    #make mean abundance if there are at least 5 cells per person and store in abundance 2
    #if there are not enough cells for at least 1 people, take all
    counts2 <- counts[,colSums(counts) > 4 ]
    if(ncol(counts2) > 0){
    props2 <- apply(counts2, 2, FUN = function(x) x/sum(x)*100)
    colSums(props2, na.rm = T)
    abundance2[,i] <- apply(props2, 1, FUN=mean, na.rm=T)
    } else { abundance2[,i] <- total / sum(total) * 100  }
}

colnames(abundance) <- colnames(abundance1) <- colnames(abundance2) <- names(split_data)
rownames(abundance) <- rownames(abundance1) <- rownames(abundance2) <- rownames(props)

#remove here the empty, xreactive and total b cells
abundance <- abundance[,-c(1,2,17)]
abundance1 <- abundance1[,-c(1,2,17)]
abundance2 <- abundance2[,-c(1,2,17)]

#make heatmaps
row_dend = hclust(dist(median_expr)) #make the ordering of clusters
hm_expr <- Heatmap(median_expr, 
         col=colorRampPalette(c("white", "darkgreen"))(50), 
         cluster_rows = row_dend, 
         border_gp = gpar(col = "black", lwd = 1))

#Change names of PS gates to PS
colnames(abundance) <- stri_replace_all_regex(colnames(abundance), pattern = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"), 
                                              replacement = c("PS1",	"PS14",	"PS15B",	"PS18C",	"PS19A",	"PS19F",	"PS23F",	"PS3",	"PS4",	"PS5",	"PS6A",	"PS6B",	"PS7F",	"PS9V"))
colnames(abundance1) <- stri_replace_all_regex(colnames(abundance1), pattern = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"), 
                                               replacement = c("PS1",	"PS14",	"PS15B",	"PS18C",	"PS19A",	"PS19F",	"PS23F",	"PS3",	"PS4",	"PS5",	"PS6A",	"PS6B",	"PS7F",	"PS9V"))
colnames(abundance2) <- stri_replace_all_regex(colnames(abundance2), pattern = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15"),
                                               replacement = c("PS1",	"PS14",	"PS15B",	"PS18C",	"PS19A",	"PS19F",	"PS23F",	"PS3",	"PS4",	"PS5",	"PS6A",	"PS6B",	"PS7F",	"PS9V"))

col.order <- c("PS1",	"PS3",	"PS4",	"PS5",	"PS6A",	"PS6B",	"PS7F",	"PS9V", "PS14",	"PS18C",	"PS19A",	"PS19F",	"PS23F", "PS15B") #change order of PS
abundance <- abundance[ , col.order]
abundance1 <- abundance1[ , col.order]
abundance2 <- abundance2[ , col.order]

#make heatmaps for the different abundance options:
hm <- Heatmap(abundance, 
                cluster_columns = F, 
                col=colorRampPalette(c("white", "darkgreen"))(50), 
                cluster_rows = row_dend, 
               border_gp = gpar(col = "black", lwd = 1))

hm1 <- Heatmap(abundance1, 
               cluster_columns = F, 
               col=colorRampPalette(c("white", "darkgreen"))(50), 
               cluster_rows = row_dend, 
               border_gp = gpar(col = "black", lwd = 1))

hm2 <- Heatmap(abundance2, 
               cluster_columns = F, 
               col=colorRampPalette(c("white", "darkgreen"))(50), 
               cluster_rows = row_dend, 
               border_gp = gpar(col = "black", lwd = 1))

#combine heatmaps
ht_list <- hm_expr + 
    rowAnnotation(percentage = anno_barplot(as.numeric(negfreqs), width = unit(3, "cm"))) +
    hm  


ht_list1 <- hm_expr + 
    rowAnnotation(percentage = anno_barplot(as.numeric(negfreqs), width = unit(3, "cm"))) +
    hm1

ht_list2 <- hm_expr + 
    rowAnnotation(percentage = anno_barplot(as.numeric(negfreqs), width = unit(3, "cm"))) +
    hm2

#generate files
pdf('.pdf', width = 8, height = 5)
draw(ht_list, column_title = "Mean abundance", column_title_gp = gpar(fontsize = 16))
draw(ht_list1, column_title = "Total abundance", column_title_gp = gpar(fontsize = 16))
draw(ht_list2, column_title = "Filtered mean abundance", column_title_gp = gpar(fontsize = 16))
dev.off() 
