#Load required packages
require(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(ggpubr)
library(viridis)
library(ggpattern)
library(logistf)
library(RColorBrewer)
library(reshape2)
library(ggrepel)   
library(broom)

#Set working directory and convert .txt files to .csv files
setwd("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis")
filelist = list.files(pattern = ".txt")
for (i in 1:length(filelist)){
  input<-filelist[i]
  output <- paste0(gsub("\\.txt$", "", input), ".csv")
  print(paste("Processing the file:", input))
  data = read.delim(input, header = TRUE)   
  write.table(data, file=output, sep=",", col.names=TRUE, row.names=FALSE)
}

Pfs47_1_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/pfs47_1_finalTab.csv",header=TRUE,sep=',')
Pfs47_2_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/pfs47_2_finalTab.csv",header=TRUE,sep=',')

Pfs47_1_processed <- process_data(Pfs47_1_raw, "pfs47_1")
Pfs47_2_processed <- process_data(Pfs47_2_raw, "pfs47_2")

Haplotypes_Pfs47_1 <- process_haplotypes(Pfs47_1_processed)
Haplotypes_Pfs47_2 <- process_haplotypes(Pfs47_2_processed)

Haplotypes_Pfs47_1 <- calculate_reads(Haplotypes_Pfs47_1)
Haplotypes_Pfs47_2 <- calculate_reads(Haplotypes_Pfs47_2)

Haplotypes_Pfs47_1$species <- ifelse(grepl("Mosq", Haplotypes_Pfs47_1$SampleID), "mosquito", "human")
Haplotypes_Pfs47_2$species <- ifelse(grepl("Mosq", Haplotypes_Pfs47_2$SampleID), "mosquito", "human")

Pfs47_FinalHaplotypes <- merge(Haplotypes_Pfs47_1, 
                               Haplotypes_Pfs47_2, 
                               by = "SampleID", 
                               all = TRUE)
Pfs47_FinalHaplotypes$species.x <- NULL
Pfs47_FinalHaplotypes$species.y <- NULL
colnames(Pfs47_FinalHaplotypes) <- c("SampleID", "Pfs47_1_hap", "Pfs47_1_reads", "Pfs47_1_MOI", "Pfs47_1_TotalReads", "Pfs47_1_Perc", "Pfs47_2_hap", "Pfs47_2_reads", "Pfs47_2_MOI", "Pfs47_2_TotalReads", "Pfs47_2_Perc")

# Split by host type and process SampleID
Pfs47_FinalHaplotypes_Human <- process_human(filter(Pfs47_FinalHaplotypes, !grepl("Mosq", SampleID)))
Pfs47_FinalHaplotypes_Mosq <- process_mosquito(filter(Pfs47_FinalHaplotypes, grepl("Mosq", SampleID)))
FINAL_Pfs47<-bind_rows(Pfs47_FinalHaplotypes_Mosq,Pfs47_FinalHaplotypes_Human)
FINAL_Pfs47$Species <- ifelse(is.na(FINAL_Pfs47$Mosquito), "human", "mosquito")


# Melt data for further analysis
melt_data <- function(data) {
  data %>%
    separate_rows(Haplotype, Reads, Percentage, sep = " ")
}

# Melt and combine data
Haplotypes_Pfs47_1_Molten <- melt_data(Haplotypes_Pfs47_1)
Haplotypes_Pfs47_1_Molten$Timepoint <- gsub("_Mosq\\d+", "", Haplotypes_Pfs47_1_Molten$SampleID)
Haplotypes_Pfs47_2_Molten <- melt_data(Haplotypes_Pfs47_2)
Haplotypes_Pfs47_2_Molten$Timepoint <- gsub("_Mosq\\d+", "", Haplotypes_Pfs47_2_Molten$SampleID)

# Only keep rows of timepoints at which there are both human and mosquito samples
Matching_Haplotypes_Pfs47_1_Molten <- Haplotypes_Pfs47_1_Molten %>%
  group_by(Timepoint) %>%
  filter("human" %in% species & "mosquito" %in% species) %>%
  ungroup()
Matching_Haplotypes_Pfs47_2_Molten <- Haplotypes_Pfs47_2_Molten %>%
  group_by(Timepoint) %>%
  filter("human" %in% species & "mosquito" %in% species) %>%
  ungroup()

# Function to determine matched or unmatched haplotypes in each species
determine_sample_type <- function(group) {
  if("human" %in% group$species && "mosquito" %in% group$species) {
    return("matching")
  } else if("human" %in% group$species) {
    return("human_only")
  } else if("mosquito" %in% group$species) {
    return("mosquito_only")
  }
}

Matching_Haplotypes_Pfs47_1_Molten <- Matching_Haplotypes_Pfs47_1_Molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
  ungroup()
Matching_Haplotypes_Pfs47_2_Molten <- Matching_Haplotypes_Pfs47_2_Molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
  ungroup()


