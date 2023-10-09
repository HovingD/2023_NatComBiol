###
## Import libraries
#


#For PS gates:
# 0 = B cells
# 1 = Cross-reactive
# 2 = PS1a
# 3 = PS1b
# 4 = PS2
# 5 = PS3
# 6 = PS5
# 7 = empty

#clean environment
rm(list=ls())

setwd('')

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
meta <- read.csv2("metadataGBS.csv", na.strings=c("","NA"))
meta$Filename <- gsub(pattern = ".csv", replacement = ".fcs", x = meta$Filename)
meta$OmiqID <- factor(meta$OmiqID, levels = unique(meta$OmiqID))
meta$Batch <- factor(meta$Batch, levels = unique(meta$Batch))
meta$Col_NonCol_EU <- factor(meta$Col_NonCol_EU, levels = unique(meta$Col_NonCol_EU))
meta$Colonization <- factor(meta$Colonization, levels = unique(meta$Colonization))
meta$ConcatBatch100K_Control <- factor(meta$ConcatBatch100K_Control, levels = unique(meta$ConcatBatch100K_Control))
meta$ConcatBatch10K_Control <- factor(meta$ConcatBatch10K_Control, levels = unique(meta$ConcatBatch10K_Control))
meta$ConcatBatch500K_Control <- factor(meta$ConcatBatch500K_Control, levels = unique(meta$ConcatBatch500K_Control))
meta$ConcatBatch50K_Control <- factor(meta$ConcatBatch50K_Control, levels = unique(meta$ConcatBatch50K_Control))
meta$Donor <- factor(trimws(meta$Donor), levels = unique(meta$Donor))
meta$SA_EU <- factor(trimws(meta$SA_EU), levels = unique(meta$SA_EU))
meta$SampleBatch_Control <- factor(trimws(meta$SampleBatch_Control), levels = unique(meta$SampleBatch_Control))
meta$Samples <- factor(meta$Samples, levels = unique(meta$Samples))
meta$SerotypePS1a <- factor(meta$SerotypePS1a, levels = unique(meta$SerotypePS1a))
meta$SerotypePS1b <- factor(meta$SerotypePS1b, levels = unique(meta$SerotypePS1b))
meta$SerotypePS2 <- factor(meta$SerotypePS2, levels = unique(meta$SerotypePS2))
meta$SerotypePS3 <- factor(meta$SerotypePS3, levels = unique(meta$SerotypePS3))
meta$SerotypePS5 <- factor(meta$SerotypePS5, levels = unique(meta$SerotypePS5))
meta$Serotype_specific_Colonization <- factor(meta$Serotype_specific_Colonization, levels = unique(meta$Serotype_specific_Colonization))
meta$Filename  <- factor(meta$Filename, 
                         levels = meta$Filename[order(meta$Donor)])

# panel
panel <- read.csv2("panelGBS.csv", na.strings=c("","NA"))

# fcs files
fcsFiles <- list.files(path = getwd(), pattern = ".fcs", full.names = TRUE)

# construct SingleCellExperiment
sce <- CATALYST::prepData(x = fcsFiles, md = meta, panel = panel, features = panel$channel_name,
                          panel_cols = list(channel = "channel_name",
                                            antigen = "marker_label",
                                            class = "marker_class"),
                          
                          md_cols = list(file = "Filename", id = "Donor", 
                                         factors = c("OmiqID", "Batch", "Col_NonCol_EU")),
                          transform = FALSE,
                          FACS = TRUE)
names(sce@assays) <- "exprs"

### 
## Import Set clustering as metadata
#
clusters <- assay(sce, i = "exprs")

#get expression of markers
markers <- data.frame(t(unlist(clusters[c(26,1:24),])))  #exclude PS (25) and fsom-metaclust-elbow-1 (26)
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
abundance <- abundance[,-c(1,2,8)]
abundance1 <- abundance1[,-c(1,2,8)]
abundance2 <- abundance2[,-c(1,2,8)]

#make heatmaps
row_dend = hclust(dist(median_expr)) #make the ordering of clusters
hm_expr <- Heatmap(median_expr, 
         col=colorRampPalette(c("white", "darkgreen"))(50), 
         cluster_rows = row_dend, 
         border_gp = gpar(col = "black", lwd = 1))

#Change names of PS gates to PS
colnames(abundance) <- stri_replace_all_regex(colnames(abundance), pattern = c("2","3","4","5","6"), replacement = c( "PSIa", "PSIb", "PSII", "PSIII", "PSV"))
colnames(abundance1) <- stri_replace_all_regex(colnames(abundance1), pattern = c("2","3","4","5","6"), replacement = c( "PSIa", "PSIb", "PSII", "PSIII", "PSV"))
colnames(abundance2) <- stri_replace_all_regex(colnames(abundance2), pattern = c("2","3","4","5","6"), replacement = c( "PSIa", "PSIb", "PSII", "PSIII", "PSV"))

#make heatmaps for the different abundance options
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
dev.off() #is required to finish creating the plot.