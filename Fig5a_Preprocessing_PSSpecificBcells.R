#This script takes out duplicated cells, separates the crossreactive cells and downsamples the file to 10.000 non-specific cells.

rm(list=ls()) #cleans the environment

setwd('')
getwd()

#loads in required packages
library(mvtnorm)
library(usethis)
library(reshape2)
library(ggplot2)
library(flowCore)
library(cowplot)
library(plyr)
library(devtools) 

#Save package version info:
sink('session info.txt')
sessionInfo()
sink()

#Go to folder where files are
setwd('')
getwd()

####################################################################
files <- list.files() #make a list of all files

sample <- unique(do.call(c, lapply(strsplit(files, '_'), function(i)i[1])))
sample3 <- sample

counter = 0 #set a counter that is used to check if this is the first sample used. (?)

#cycle through each of the samples
for(vol in sample3){
  
  file <- list.files(pattern=vol)
  
  #read in data and fix column names
  data <- read.FCS(file[1], transformation  = F, truncate_max_range = F, ignore.text.offset = TRUE)
  dat2 <- data.frame(data@exprs) #make a data frame of the data: cols = markers and rows = cells
  
  #get the marker names only out and remove the fluorochromes
  desc <- data@parameters@data$desc
  desc[is.na(desc)] <- data@parameters@data$name[is.na(data@parameters@data$desc)]
  colnames(dat2) <- desc #save the markers as data frame column names
  dat2 <- dat2[,colnames(dat2) != 'Event #'] #remove the event parameter if it's there
  
  
  dat2$ID <- paste0(dat2$IgM, '_', dat2$IgG) #create an ID of the cells based on their IgM, IgG (and IgD if needed?) values
  print(paste0(vol, ' has no duplicated B cells: '))
  print(sum(duplicated(dat2$ID)) == 0) #needs to be true so that each cell is unique. If this is FALSE then markers should be added or there are duplicates
  dat2 <- dat2[!duplicated(dat2$ID),] #remove duplicate rows (cell duplicates)
  
  #add in new rows
  dat2$Vol <- vol
  dat2$GBS <- NA
  
  #load in the specific fcs files and add a tag with GBS specificity
  for(i in 2:length(file)){ #use for loop to look at each specific exported 
    data <- read.FCS(file[i], transformation  = F, truncate_max_range = F, ignore.text.offset = TRUE)
    dat3 <- data.frame(data@exprs) #get the data frame for these specific cells and save as dat3
    
    if(nrow(dat3)>0){ #checks if there are any cells in there
      
      #fix the names again as above
      desc <- data@parameters@data$desc
      desc[is.na(desc)] <- data@parameters@data$name[is.na(data@parameters@data$desc)]
      colnames(dat3) <- desc
      dat3 <- dat3[,colnames(dat3) != 'Event #']
      
      dat3 <- dat3[!duplicated(dat3),] #remove duplicate rows
      dat3$ID <- paste0(dat3$IgM, '_', dat3$IgG) #check for ID
               
      if(sum(dat2$ID %in% dat3$ID) == nrow(dat3)){ #check that all the cells are in the parent fcs file
        
        #in the cells with already a PS assigned, add the specific cells onto them
        dat2$GBS[dat2$ID %in% dat3$ID & !is.na(dat2$GBS)] <- paste(dat2$GBS[dat2$ID %in% dat3$ID & !is.na(dat2$GBS)], 
                                                                   strsplit(file[i], '_')[[1]][length(strsplit(file[i], '_')[[1]])], 
                                                                   sep = '_') #add on if co-pos
        
        #in the cells with not yet a PS assigned, name them as the specific cells
        dat2$GBS[dat2$ID %in% dat3$ID & is.na(dat2$GBS)] <- strsplit(file[i], '_')[[1]][length(strsplit(file[i], '_')[[1]])] #assign the names if not existing yet
        
        dat2$GBS <- gsub('.fcs', '', dat2$GBS)
      } else{print('not all cells from the specific population are found in the B cells total')}
    }
  }
  table(dat2$GBS) #check if the correct files are added
  
  #fix the names of the PS
  names <- file[2:7] #Make sure to adjust this to total amount of files (Bcells, empty and 5 GBS PS. We now exclude Bcells)
  names <- do.call(c, lapply(strsplit(names, '_'), function(i) i[length(strsplit(names, '_')[[1]])]))
  names <- gsub('.fcs', '', names)  
  
  #check for crossreactives
  dat2$PS <- dat2$GBS
  dat2$PS[is.na(dat2$GBS)] <- '0.Negative'
  xreactive <- which(!dat2$GBS %in% names & !is.na(dat2$GBS))
  dat2$PS[xreactive] <- 'Crossreactive'
  table(dat2$PS)
  
  #add 1 to the counter and check if it's the first file
  counter = counter +1
  if(counter == 1){
    
    counts <- data.frame(table(dat2$PS))
    colnames(counts) <- c('Population', vol)
    
    #take only 10,000 of the negative cells and all of the positive ones
    down = 10000
    dat4 <- dat2[sample(which(is.na(dat2$GBS)), down),] #Get a random cell (sample function), do this 10000 times.
    dat_all <- rbind(dat4, dat2[!is.na(dat2$GBS),])
    
    
    if(sum(counts[,2]) == nrow(dat2)){print(paste0(vol, ' ok'))}
    
  } else {
    
    #make the counts for the new file and merge them with the other counts
    counts2 <-  data.frame(table(dat2$PS))
    colnames(counts2) <- c('Population', vol)
    counts <-  merge(counts2, counts, all = T, by = 'Population')
    
    down = 10000
    dat4 <- dat2[sample(which(is.na(dat2$GBS)), down),]
    dat5 <- rbind(dat4, dat2[!is.na(dat2$GBS),])
         
    #paste the different samples together
    dat_all <- rbind(dat_all, dat5)
    
    if(sum(counts2[,2]) == nrow(dat2)){print(paste0(vol, ' ok'))}
    
  }
  
}

counts[is.na(counts)] <- 0 #set msissing PS to 0 

#make a new folder and set as write out
#setwd('')
dir.create('')
setwd('')
getwd()

write.table(counts, 'counts_specific_allCSV.csv', row.names = F, quote = F, sep = ',')

table(dat_all$PS)

dat_all2 <- dat_all[dat_all$Vol != 'Batch16_FMT',] #if there is a FMT sample this is removed. just change name of FMT sample

#Save all donors combined in one file:
write.table(dat_all2, 'specific_cells_allCSV.csv', row.names = F, quote = F, sep = ',') #it save the specific cells

#separate in individual files
split.data <- split(dat_all2, dat_all2$Vol) 

for(j in (1:length(split.data))) {
  
  one <- split.data[[j]]
  
  write.table(one, paste0(one$Vol[1], 'specific_cells_allCSV.csv'), sep = ',', row.names = F) 
  
}