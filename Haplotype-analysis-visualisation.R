---
title: "_Plasmodium falciparum_ transmissibility in asymptomatic gametocyte carriers in Mali "
author: "Leen Vanheer"

---
  
#Load libraries
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
library(dunn.test)
library(ape)
library(pegas)
library(Biostrings)
library(lubridate)
require(dplyr)

# set seed for reproducibility
set.seed(0)

# determine color palettes
pal1<-c("mosquito_only"="#f8997c","human_only"= "#a891cf", "matching" = "lightgrey")
pal2<-c("lightgrey", "grey48")
pal3<-c("monoclonal"="#9AC77B","polyclonal"='#7EC2DE')
pal4 <- c("#8dd3c6","#fccde5", "#f98073", "#80b0d3")
pal5 <- c("#4c9184", "#bf6b96", "#a8392d","#3373a1")
pal6 <- c("#ebccb5","#A2A77F", "#8CBCB9", "#DB8A74", "#972D07","#C19AAF", "#D1AC00")

################################################### Functions ###################################################
process_data <- function(data, type) {
  # Common preprocessing
  data <- data %>%
    filter(grepl(type, Haplotype)) %>%
    mutate(SampleID = gsub("(.*)_.*", "\\1", SampleID),
           Reads = as.numeric(Reads)) %>%
    group_by(SampleID, SampleName) %>%
    dplyr::summarise(Haplotype = paste(Haplotype, collapse = " "),
              Reads = paste(Reads, collapse = " ")) %>%
    ungroup()
  
  # Initialize an empty data frame for results
  results <- data.frame(SampleID = character(),
                        Matching_Values = character(),
                        Average_Reads = numeric())
  
  # Loop through unique samples
  for (sample in unique(data$SampleID)) {
    subset_data <- filter(data, SampleID == sample)
    
    if (nrow(subset_data) == 2) {
      # Extract haplotypes and reads for both reps
      haplotypes_rep1 <- unlist(strsplit(filter(subset_data, SampleName == "Rep1")$Haplotype, " "))
      haplotypes_rep2 <- unlist(strsplit(filter(subset_data, SampleName == "Rep2")$Haplotype, " "))
      reads_rep1 <- as.numeric(unlist(strsplit(filter(subset_data, SampleName == "Rep1")$Reads, " ")))
      reads_rep2 <- as.numeric(unlist(strsplit(filter(subset_data, SampleName == "Rep2")$Reads, " ")))
      
      # Identify matching haplotypes and calculate average reads
      matching_haplotypes <- intersect(haplotypes_rep1, haplotypes_rep2)
      avg_reads <- sapply(matching_haplotypes, function(haplotype) {
        mean(c(reads_rep1[which(haplotypes_rep1 == haplotype)], reads_rep2[which(haplotypes_rep2 == haplotype)]))
      })
      
      # Append to results data frame
      results <- rbind(results, data.frame(SampleID = rep(sample, length(matching_haplotypes)),
                                           Matching_Values = matching_haplotypes,
                                           Average_Reads = avg_reads))
    }
  }
  
  return(results)
}


# Function to process haplotype data and calculate MOI
process_haplotypes <- function(data) {
  data %>%
    group_by(SampleID) %>%
    dplyr::summarise(Haplotype = paste(Matching_Values, collapse = " "),
              Reads = paste(Average_Reads, collapse = " ")) %>%
    mutate(MOI = str_count(Haplotype, " ") + 1)
}

# Function to calculate total reads and percentage
calculate_reads <- function(data) {
  data %>%
    mutate(TotalReads = sapply(strsplit(Reads, " "), function(x) sum(as.numeric(x))),
           Percentage = sapply(Reads, function(x) {
             reads_list <- as.numeric(strsplit(x, " ")[[1]])
             if (sum(reads_list) == 0) {
               perc_list <- rep(0, length(reads_list))
             } else {
               perc_list <- round((reads_list / sum(reads_list)) * 100, 2)
             }
             return(paste(perc_list, collapse = " "))
           }))
}

# Function to split and process human haplotypes for further data analysis
process_human <- function(data) {
  data_split <- data %>%
    dplyr::mutate(Individual = str_extract(SampleID, "^[^_]+"),
           Day = str_extract(SampleID, "_([^_]+)$"),
           Day = str_replace_all(Day, "_", ""), 
           Day = as.numeric(str_sub(Day, 4)),
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)),
           Timepoint = SampleID,
           MOI_Combined = pmax(trap_MOI,csp_MOI, na.rm = TRUE)) %>%
    dplyr::select(Individual, Day, Timepoint, everything())
  
  return(data_split)
}

# Function to split and process mosquito haplotypes for further data analysis
process_mosquito <- function(data) {
  data_split <- data %>%
    dplyr::mutate(Individual = str_extract(SampleID, "^[^_]+"),
           Day = str_extract(SampleID, "_([^_]*)_"),  # Capture value between underscores
           Mosquito = str_extract(SampleID, "_([^_]+)$"),
           Day = str_replace_all(Day, "_", ""), 
           Timepoint = paste0(Individual,"_", Day),
           Mosquito = str_replace_all(Mosquito, "_", ""),
           Day = as.numeric(str_sub(Day, 4)),
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)),
           MOI_Combined = pmax(trap_MOI,csp_MOI, na.rm = TRUE)) %>%
    
    dplyr::select(Individual, Day, Mosquito, Timepoint, everything())
  
  return(data_split)
}

# Function to melt data for further analysis
melt_data <- function(data) {
  data %>%
    separate_rows(Haplotype, Reads, Percentage, sep = " ")%>%
    dplyr::mutate(Timepoint = gsub("_Mosq\\d+", "", SampleID),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
    left_join(csp_trap_finalhaplotypes[, c("SampleID", "Individual")], by = "SampleID")
}

# Function to determine matched or unmatched haplotypes in each species
determine_sample_type <- function(group) {
  type <- if("human" %in% group$species && "mosquito" %in% group$species) {
    "matching"
  } else if("human" %in% group$species) {
    "human_only"
  } else if("mosquito" %in% group$species) {
    "mosquito_only"
  } else {
    NA_character_  # Handle cases where none of the conditions are met
  }
  
  # Convert to factor with specified levels
  factor(type, levels = c("human_only", "mosquito_only", "matching"))
}

# Define a function to perform the matching and transformation
process_data2 <- function(data) {
  data %>%
    group_by(Timepoint) %>%
    dplyr::filter("human" %in% species & "mosquito" %in% species) %>%
    ungroup() %>%
    group_by(Timepoint, Haplotype) %>%
    dplyr::mutate(Comparison = determine_sample_type(cur_data()),
           Day = str_extract(Timepoint, "_([^_]+)$"),  
           Day = str_replace_all(Day, "_", ""), 
           Day = as.numeric(str_sub(Day, 4)),
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
    ungroup()
}

process_data3 <- function(data) {
  data %>%
    group_by(Day, Individual) %>%
    dplyr::filter("human" %in% species & "mosquito" %in% species) %>%
    ungroup() %>%
    group_by(Day, Individual, Haplotype) %>%
    dplyr::mutate(Comparison = determine_sample_type(cur_data())) %>%
    ungroup()
}

calculate_haplotype_fractions <- function(data, haplotype_column) {
  # Step 1: Split haplotypes and expand the DataFrame
  expanded_df <- data %>%
    separate_rows(!!sym(haplotype_column), sep = " ") %>%
    filter(!!sym(haplotype_column) != "")
  
  # Step 2: Count unique samples for each haplotype and species
  haplotype_counts <- expanded_df %>%
    group_by(!!sym(haplotype_column), Species) %>%
    dplyr::summarise(count = n_distinct(SampleID), .groups = 'drop')
  
  # Step 3: Calculate the percentage for each species
  total_samples_human <- n_distinct(data$SampleID[data$Species == "human"])
  total_samples_mosquito <- n_distinct(data$SampleID[data$Species == "mosquito"])
  
  haplotype_counts %>%
    mutate(frac_sample = case_when(
      Species == "human" ~ (count / total_samples_human) * 100,
      Species == "mosquito" ~ (count / total_samples_mosquito) * 100
    )) %>%
    dplyr::select(haplotype = !!sym(haplotype_column), Species, frac_sample)
}

# Function to extract sequence by haplotype and write to file
write_haplotype_sequences <- function(sample_id, haplotype, day, species, sequence, file_path) {
  # Construct the FASTA header
  header <- paste0(">", sample_id, "_", haplotype, "_", day, "_",species)
  # Combine header and sequence
  fasta_entry <- paste(header, sequence, sep="\n")
  # Write to file
  cat(fasta_entry, file=file_path, append=TRUE, sep="\n")
}

################################################### Data processing ###################################################
######## Set working directory and convert HaplotypR output .txt files to .csv files
setwd("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis")
filelist = list.files(pattern = ".txt")
for (i in 1:length(filelist)){
  input<-filelist[i]
  output <- paste0(gsub("\\.txt$", "", input), ".csv")
  print(paste("Processing the file:", input))
  data = read.delim(input, header = TRUE)   
  write.table(data, file=output, sep=",", col.names=TRUE, row.names=FALSE)
}


######## Load and process the data
trap_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/TRAP_finalTab.csv",header=TRUE,sep=',')
csp_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/CSP_finalTab.csv",header=TRUE,sep=',')
clinicaldata<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Paper_AmpSeq/ClinicalData.csv",header=TRUE,sep=',')
oocystdata<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Paper_AmpSeq/Oocyst_data.csv",header=TRUE,sep=',')
TRAP_haplotypes <- readDNAStringSet("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork/trap_HaplotypeSeq_merge.fasta", format="fasta")
CSP_haplotypes <- readDNAStringSet("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork/csp_HaplotypeSeq_merge.fasta", format="fasta")
labstrain_ratios_data<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Paper_AmpSeq/Labstrain_Ratios_Data.csv",header=TRUE,sep=',')
variants_Exp51<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/variants_Exp51.csv")
variants_Exp52<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/variants_Exp52.csv")
coverage_Exp51<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/coverage_report_Exp51.csv")
coverage_Exp52<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/coverage_report_Exp52.csv")
missing_positions_Exp51<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/missing_positions_report_Exp51.csv")
missing_positions_Exp52<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/missing_positions_report_Exp52.csv")

# Initial processing
trap_processed <- process_data(trap_raw, "trap")
csp_processed <- process_data(csp_raw, "csp")

# Process Matching Haplotypes
csp_haplotypes <- process_haplotypes(csp_processed)
trap_haplotypes <- process_haplotypes(trap_processed)

# Calculate total reads and percentages + assign species
csp_haplotypes <- calculate_reads(csp_haplotypes)
trap_haplotypes <- calculate_reads(trap_haplotypes)
csp_haplotypes$species <- ifelse(grepl("Mosq", csp_haplotypes$SampleID), "mosquito", "human")
trap_haplotypes$species <- ifelse(grepl("Mosq", trap_haplotypes$SampleID), "mosquito", "human")

# Combine markers and format haplotype data
finalhaplotypes <- merge(trap_haplotypes, csp_haplotypes, by = "SampleID", all = TRUE)
finalhaplotypes$species.x <- NULL
finalhaplotypes$species.y <- NULL
colnames(finalhaplotypes) <- c("SampleID", "trap_hap", "trap_reads", "trap_MOI", "trap_TotalReads", "trap_Perc", "csp_hap", "csp_reads", "csp_MOI", "csp_TotalReads", "csp_Perc")

# Split by host type and process SampleID, then combine again
finalhaplotypes_human <- process_human(filter(finalhaplotypes, !grepl("Mosq", SampleID)))
finalhaplotypes_mosq <- process_mosquito(filter(finalhaplotypes, grepl("Mosq", SampleID)))
csp_trap_finalhaplotypes<-bind_rows(finalhaplotypes_mosq,finalhaplotypes_human)
csp_trap_finalhaplotypes$Species <- ifelse(is.na(csp_trap_finalhaplotypes$Mosquito), "human", "mosquito")
csp_trap_finalhaplotypes_human<-csp_trap_finalhaplotypes[csp_trap_finalhaplotypes$Species == "human",]
csp_trap_finalhaplotypes_mosq<-csp_trap_finalhaplotypes[csp_trap_finalhaplotypes$Species == "mosquito",]

# Melt data
csp_haplotypes_molten <- melt_data(csp_haplotypes)
trap_haplotypes_molten <- melt_data(trap_haplotypes)

# Only keep rows of timepoints at which there are both human and mosquito samples + determine matched or unmatched haplotypes in each species
csp_matching_haplotypes_molten <- process_data2(csp_haplotypes_molten)
trap_matching_haplotypes_molten <- process_data2(trap_haplotypes_molten)

# Calculate haplotype fractions
trap_frac <- calculate_haplotype_fractions(csp_trap_finalhaplotypes, "trap_hap")
csp_frac <- calculate_haplotype_fractions(csp_trap_finalhaplotypes, "csp_hap")
csp_trap_frac <- bind_rows(trap_frac, csp_frac)

# Combine clinical data and haplotype data
clinicaldata$SampleID<-paste(clinicaldata$studycode,"_Day",clinicaldata$studyvisit,sep="")
merge_clinical_haplotypes <- csp_trap_finalhaplotypes %>%
  merge(clinicaldata, by = "SampleID") %>%
  mutate(percentagemosqinfected = (mosq_pos / mosq_total) * 100)
csp_merge_clinical_haplotypes_molten <- csp_haplotypes_molten %>%
  merge(clinicaldata, by = "SampleID") %>%
  mutate(
    percentagemosqinfected = (mosq_pos / mosq_total) * 100,
    transmission = ifelse(percentagemosqinfected > 0, 1, 0)
  )
trap_merge_clinical_haplotypes_molten <- trap_haplotypes_molten %>%
  merge(clinicaldata, by = "SampleID") %>%
  mutate(
    percentagemosqinfected = (mosq_pos / mosq_total) * 100,
    transmission = ifelse(percentagemosqinfected > 0, 1, 0)
  )
csp_merge_clinical_matching_haplotypes_molten <- csp_matching_haplotypes_molten %>%
  merge(clinicaldata, by = "SampleID") %>%
  mutate(
    percentagemosqinfected = (mosq_pos / mosq_total) * 100,
    transmission = ifelse(percentagemosqinfected > 0, 1, 0)
  )
trap_merge_clinical_matching_haplotypes_molten <- trap_matching_haplotypes_molten %>%
  merge(clinicaldata, by = "SampleID") %>%
  mutate(
    percentagemosqinfected = (mosq_pos / mosq_total) * 100,
    transmission = ifelse(percentagemosqinfected > 0, 1, 0)
  )

# Combine oocyst data and haplotype data
oocystdata$SampleID<-paste(oocystdata$studycode,"_Day",oocystdata$studyvisit,sep="")
transformed_oocystdata <- data.frame()
for (i in 1:nrow(oocystdata)) {
  row <- oocystdata[i,]
  sample_id <- as.character(row$SampleID)
  
  # Iterate over each pair of mosqnbr# and oocynbr# columns
  for (j in 1:75) { # Adjust the range according to the number of mosqnbr# and oocynbr# columns
    mosq_col <- paste0("mosqnbr", j)
    oocy_col <- paste0("oocynbr", j)
    
    if (!is.na(row[[mosq_col]]) && row[[mosq_col]] == 1) {
      # Create a new row with SampleID_Mosq# format and the corresponding oocyst number
      new_row <- data.frame(
        SampleID = paste0(sample_id, "_Mosq", j),
        oocysts = row[[oocy_col]]
      )
      
      # Add the new row to the transformed dataframe
      transformed_oocystdata <- rbind(transformed_oocystdata, new_row)
    }
  }
}
merge_oocyst_haplotypes<-merge(csp_trap_finalhaplotypes,transformed_oocystdata,by="SampleID")
transformed_oocystdata$Individual_ID <- gsub("(.*)_(.*)", "\\1", transformed_oocystdata$SampleID)
merge_oocyst_individual_haplotypes <- merge(x = csp_trap_finalhaplotypes_human, y = transformed_oocystdata, by.x = 'SampleID', by.y = 'Individual_ID', all.y = TRUE)
merge_oocyst_individual_haplotypes <- merge_oocyst_individual_haplotypes[!is.na(merge_oocyst_individual_haplotypes$MOI_Combined), ]

# Combine oocyst, clinical and haplotype data
merge_oocyst_clinical_individual_haplotypes<-merge(merge_oocyst_individual_haplotypes,clinicaldata,by="SampleID")

# Edit drug resistance data
variants_Exp52<-variants_Exp52 %>%
  mutate(uniqueid = sub("^([^_]+_[^_]+_)(.*)", "\\2", sample_id)) %>%
  group_by(genome_pos, gene, protein_change) %>%
  mutate(species="human",
         mosqid=NA)%>%
  filter(!grepl("fs", protein_change))%>%
  ungroup() %>%
  dplyr::select(1:10, 13, 11, 12)

variants_Exp51<-variants_Exp51 %>%
  mutate(mosqid = sub("^([^_]+_[^_]+_)(.*)", "\\2", sample_id),
         uniqueid = sub("^([^_]+).*", "\\1",  mosqid)) %>%
  group_by(genome_pos, gene, protein_change) %>%
  mutate(species="mosquito")%>%
  filter(!grepl("fs", protein_change))

combined_variants<-rbind(variants_Exp52,variants_Exp51)
combined_coverage<-rbind(coverage_Exp51,coverage_Exp52)
combined_missing_positions<-rbind(missing_positions_Exp51,missing_positions_Exp52) %>%
  mutate(sample_id = Sample.Name,
         gene=Gene,
         protein_change = Variants) %>%
  dplyr::select(sample_id,gene,protein_change)
################################################### Data visualisation ############################################################

######## Figure 1 - Polyclonal infections have increased gametocyte densities and cause higher oocyst density infections ########
#Gametocytemia Monoclonal versus polyclonal
a<-merge_clinical_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined == 1, "monoclonal", "polyclonal")) %>%
  group_by(MonoPoly) %>%
  dplyr::summarise(
    Mean = mean(totalgct_ul, na.rm = TRUE),
    Min = min(totalgct_ul, na.rm = TRUE),
    Max = max(totalgct_ul, na.rm = TRUE),
    totalgct_ul = list(totalgct_ul)) %>%
  unnest(totalgct_ul) %>%
  ggplot(aes(x = as.factor(MonoPoly), y = totalgct_ul)) + 
  geom_violin(aes(fill=MonoPoly),alpha = 0.5)+
  #stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")+
  geom_jitter(position = position_jitter(width = 0.3, height = 0.2), alpha = 0.5) +
  geom_point(aes(y = Mean), color = "#f8997c", size = 2) +
  scale_y_log10() + 
  annotation_logticks(sides = "l") +
  scale_fill_manual(values=pal3)+
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "MOI human blood", y = "Gametocytes per µL") 

# Perform the Wilcoxon rank-sum test
merge_clinical_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined == 1, "monoclonal", "polyclonal"))%>%
  wilcox.test(totalgct_ul ~ MonoPoly, data = ., exact = FALSE)



#Infection rates Monoclonal versus polyclonal
b<-merge_clinical_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  group_by(MonoPoly) %>%
  dplyr::summarise(
    Median_InfRate = median(percentagemosqinfected, na.rm = TRUE),
    IQR_Lower = quantile(percentagemosqinfected, probs = 0.25, na.rm = TRUE),
    IQR_Upper = quantile(percentagemosqinfected, probs = 0.75, na.rm = TRUE),
    percentagemosqinfected = list(percentagemosqinfected)
  ) %>%
  unnest(percentagemosqinfected) %>%
  ggplot(aes(x=as.factor(MonoPoly), y=percentagemosqinfected)) + 
  geom_violin(aes(fill=MonoPoly),alpha = 0.5)+
  geom_boxplot(width=0.03,fill="grey") +
  theme_classic()+
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100))+
  #geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  #geom_point(aes(y = Median_InfRate), color = "#f8997c", size = 2) +
  theme(legend.position = "none")+
  scale_fill_manual(values=pal3)+
  labs(x= "MOI human blood",y="Mosquito infection rate (%)")

# Perform the Wilcoxon rank-sum test
merge_clinical_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined == 1, "monoclonal", "polyclonal")) %>%
  wilcox.test(percentagemosqinfected ~ MonoPoly, data = ., exact = FALSE)


#Oocyst density Monoclonal versus polyclonal
c<-merge_oocyst_individual_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  group_by(MonoPoly) %>%
  dplyr::summarise(
    Median_oocysts = median(oocysts, na.rm = TRUE),
    IQR_Lower = quantile(oocysts, probs = 0.25, na.rm = TRUE),
    IQR_Upper = quantile(oocysts, probs = 0.75, na.rm = TRUE),
    oocysts = list(oocysts)
  ) %>%
  unnest(oocysts) %>%
  ggplot(aes(x=as.factor(MonoPoly), y=oocysts)) + 
  geom_violin(aes(fill=MonoPoly),alpha = 0.5)+
  geom_boxplot(width=0.03,fill="grey") +
  scale_y_log10(breaks=c(1,5,10,50,100,200))+
  annotation_logticks(sides = "l") +
  theme_classic()+
  scale_fill_manual(values=pal3)+
  #geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  #geom_point(aes(y = Median_oocysts), color = "#f8997c", size = 2) +
  theme(legend.position = "none")+
  labs(x= "MOI human blood",y="Oocyst density")


# Perform the Wilcoxon rank-sum test
merge_oocyst_individual_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined == 1, "monoclonal", "polyclonal")) %>%
  wilcox.test(oocysts ~ MonoPoly, data = ., exact = FALSE)



#Infection rates Monoclonal versus polyclonal normalised by GC
d<-merge_clinical_haplotypes %>%
  mutate(normalized_infection_rate = (percentagemosqinfected / totalgct_ul) * 100,
         MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  ggplot(aes(x=as.factor(MonoPoly), y=normalized_infection_rate)) + 
  geom_violin(aes(fill=MonoPoly),alpha = 0.5)+
  geom_boxplot(width=0.03,fill="grey") +
  theme_classic()+
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100))+
  theme(legend.position = "none")+
  scale_fill_manual(values=pal3)+
  labs(x= "MOI human blood",y="Normalised mosquito infection rate (%)")

# Perform the Wilcoxon rank-sum test
merge_clinical_haplotypes %>%
  mutate(normalized_infection_rate = (percentagemosqinfected / totalgct_ul) * 100,
         MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  wilcox.test(normalized_infection_rate ~ MonoPoly, data = ., exact = FALSE)


#Oocyst density Monoclonal versus polyclonal normalised by GC
e<-merge_oocyst_clinical_individual_haplotypes %>%
  mutate(normalized_oocysts = (oocysts / totalgct_ul) * 100,
         MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  ggplot(aes(x=as.factor(MonoPoly), y=normalized_oocysts)) + 
  geom_violin(aes(fill=MonoPoly),alpha = 0.5)+
  geom_boxplot(width=0.03,fill="grey") +
  theme_classic()+
  scale_fill_manual(values=pal3)+
  theme(legend.position = "none")+
  labs(x= "MOI human blood",y="Normalized oocyst density (%)")

# Perform the Wilcoxon rank-sum test
merge_oocyst_clinical_individual_haplotypes %>%
  mutate(normalized_oocysts = (oocysts / totalgct_ul) * 100,
         MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  wilcox.test(normalized_oocysts ~ MonoPoly, data = ., exact = FALSE)


Figure1<-a+b+c+d+e+ plot_layout(ncol=3, widths=c(1, 1,1))

ggsave("Figure1.pdf", Figure3, width=14, height=6)
ggsave("Figure1A.pdf", a , width=4.6, height=6)

######## Figure 2A - Haplotype prevalence in human vs mosquito samples ########
csp_trap_frac %>%
  group_by(haplotype) %>%
  dplyr::summarise(frac_sample_human = sum(frac_sample[Species == "human"], na.rm = TRUE),
            frac_sample_mosquito = sum(frac_sample[Species == "mosquito"], na.rm = TRUE),
            .groups = 'drop') %>%
  mutate(marker = sub("-.*", "", haplotype),
         difference = frac_sample_human - frac_sample_mosquito,
         hap_lab = ifelse(haplotype == 'csp-22', "RTS,S/AS01\nvaccine haplotype", NA),
         hap_lab2 = ifelse(haplotype == 'trap-8', "3D7 trap", NA)) %>%
  ggplot(aes(x = frac_sample_human, y = frac_sample_mosquito, col = difference)) +
  geom_abline(intercept = 0, col = 'grey80') +
  geom_point(aes(shape = marker), size = 3) +
  scale_shape_manual(values = c(17, 19)) +
  geom_text_repel(aes(label = haplotype), color = 'black') +
  labs(x = 'Haplotype prevalence in human samples', y = 'Haplotype prevalence in mosquito samples',
       shape = 'Marker', col = 'Difference in prevalence') +
  coord_cartesian(xlim = c(0,40), ylim = c(0,40)) + 
  scale_colour_gradient2(low = "red", mid = "white", high = "darkblue", midpoint = 0) +
  theme_light()
ggsave("Haplotype_prevalence_humans_mosquitoes.pdf", width=7, height=5)




######## Figure 2B - Percentage matching clones/human only/mosquito only at each timepoint ########
a<- ggplot(csp_matching_haplotypes_molten[csp_matching_haplotypes_molten$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  xlab("Days after treatment initiation")+
  ylab("Haplotype percentage")+
  ggtitle("csp")+
  theme(legend.position="none",
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values=pal1)
b<- ggplot(trap_matching_haplotypes_molten[trap_matching_haplotypes_molten$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  #facet_wrap(~Matching_FinalHaplotypes_ALL$Host)
  xlab("Days after treatment initiation")+
  ylab("")+
  scale_fill_manual(values=pal1, name =" ")+
  ggtitle("trap")+
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))
combined_plot <- a + b + plot_layout(ncol=2, widths=c(1.1, 1))
ggsave("Fig1B_Matching_Haplotypes_Percentages.pdf", combined_plot, width=14, height=6)

#Combined
csp_data_prepared <- csp_matching_haplotypes_molten %>%
  filter(Day %in% c("0", "2", "7", "14")) %>%
  group_by(Day, Comparison) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  mutate(Source = "csp")

trap_data_prepared <- trap_matching_haplotypes_molten %>%
  filter(Day %in% c("0", "2", "7", "14")) %>%
  group_by(Day, Comparison) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  mutate(Source = "trap") 

combined_data <- bind_rows(csp_data_prepared, trap_data_prepared) %>%
  group_by(Day, Comparison) %>%
  dplyr::summarise(AvgCount = mean(Count), .groups = "drop")

ggplot(combined_data, aes(x = as.factor(Day), y = AvgCount, fill = factor(Comparison))) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  xlab("Days after treatment initiation") +
  ylab("Haplotype Percentage (%)")+
  scale_fill_manual(values = pal1) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), 
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 2), # Add border
        plot.margin = margin(5.5, 5.5, 5.5, 5.5) # Optional: Adjust plot margin
  )
ggsave("Fig2D_Matching_Haplotypes_Percentages_bothMarkers.pdf",width=7, height=5)

######## Figure 2C - Evidence of transmission of minority clones in polyclonal infections ########

#Check if matching haplotypes in multiclonal individuals more often have high or low coverage at Day0?
test<-csp_merge_clinical_matching_haplotypes_molten %>%
  filter(MOI > 2,
         species == "human")

# Filtering the human-only rows
human_only_df <- test %>%
  filter(Comparison == "human_only") %>%
  mutate(Percentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day,.groups='keep')%>%
  dplyr::summarise(TotalPercentage = sum(Percentage, na.rm = TRUE))

# Filtering the matching rows
non_human_only_df <- test %>%
  filter(Comparison != "human_only") %>%
  mutate(TotalPercentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day)

# Binding the two dataframes together
filtered_df <- bind_rows(human_only_df, non_human_only_df) %>%
  filter(Day%in% c(0,2,7,14))


ggplot(filtered_df, aes(x=Comparison,y=as.numeric(TotalPercentage))) + 
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 2) +
  theme_classic()+
  xlab("Type")+
  ylab("Percentage of total reads")+
  facet_wrap(~Day,ncol=4)+
  ggtitle("CSP")
ggsave("Percentage_Reads_Human_only_vs_matching_CSP.pdf", width=4, height=3)

######## Figure 2D - Percentage matching clones/human only/mosquito comparing day 0 mosquito clones to all human timepoints ########

#csp
csp_haplotypes_molten_filter<-csp_haplotypes_molten %>%
  filter(species=="human" | Day==0)

mosquito_day0 <- csp_haplotypes_molten_filter %>%
  filter(species == "mosquito", Day == 0)

human_timepoints_individuals <- csp_haplotypes_molten_filter %>%
  filter(species == "human") %>%
  dplyr::select(SampleID, Day, Timepoint, Individual) %>%
  distinct()

# Replicate Day 0 mosquito samples for each Day present in the human samples
expanded_mosquito_samples <- mosquito_day0 %>%
  # Perform a full join to replicate mosquito data for each human timepoint
  full_join(human_timepoints_individuals, by = "Individual") %>%
  mutate(species = "mosquito",
         SampleID = SampleID.x,
         Day = Day.y) %>% 
  # Retain necessary columns and arrange for clarity
  dplyr::select(SampleID, Haplotype, Reads, MOI, TotalReads, Percentage, species, Day, Individual) %>%
  arrange(SampleID, Individual, Day)

csp_plotme <- bind_rows(csp_haplotypes_molten_filter, expanded_mosquito_samples) %>%
  arrange(SampleID, Individual, Day) %>%
  filter(!is.na(SampleID)) %>%
  process_data3()

csp_plotme<-csp_plotme[-which(is.na(csp_plotme$Timepoint) & csp_plotme$Day == 0),]

#trap
trap_haplotypes_molten_filter<-trap_haplotypes_molten %>%
  filter(species=="human" | Day==0)

mosquito_day0 <- trap_haplotypes_molten_filter %>%
  filter(species == "mosquito", Day == 0)

human_timepoints_individuals <- trap_haplotypes_molten_filter %>%
  filter(species == "human") %>%
  dplyr::select(SampleID, Day, Timepoint, Individual) %>%
  distinct()

# Replicate Day 0 mosquito samples for each Day present in the human samples
expanded_mosquito_samples <- mosquito_day0 %>%
  # Perform a full join to replicate mosquito data for each human timepoint
  full_join(human_timepoints_individuals, by = "Individual") %>%
  mutate(species = "mosquito",
         SampleID = SampleID.x,
         Day = Day.y) %>% 
  # Retain necessary columns and arrange for clarity
  dplyr::select(SampleID, Haplotype, Reads, MOI, TotalReads, Percentage, species, Day, Individual) %>%
  arrange(SampleID, Individual, Day)

trap_plotme <- bind_rows(trap_haplotypes_molten_filter, expanded_mosquito_samples) %>%
  arrange(SampleID, Individual, Day) %>%
  filter(!is.na(SampleID)) %>%
  process_data3()

trap_plotme<-trap_plotme[-which(is.na(trap_plotme$Timepoint) & trap_plotme$Day == 0),]

#Combined
csp_data_prepared <- csp_plotme %>%
  filter(Day %in% c("0", "2", "7", "14")) %>%
  group_by(Day, Comparison) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  mutate(Source = "csp")

trap_data_prepared <- trap_plotme %>%
  filter(Day %in% c("0", "2", "7", "14")) %>%
  group_by(Day, Comparison) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  mutate(Source = "trap") 

combined_data <- bind_rows(csp_data_prepared, trap_data_prepared) %>%
  group_by(Day, Comparison) %>%
  dplyr::summarise(AvgCount = mean(Count), .groups = "drop")

ggplot(combined_data, aes(x = as.factor(Day), y = AvgCount, fill = factor(Comparison))) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  xlab("Days after treatment initiation") +
  ylab("Haplotype Percentage (%)")+
  scale_fill_manual(values = pal1) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12), 
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 2), # Add border
        plot.margin = margin(5.5, 5.5, 5.5, 5.5) # Optional: Adjust plot margin
  )
ggsave("Fig2E_Matching_Haplotypes_Percentages_Day0_bothMarkers.pdf",width=7, height=5)


######## Figure 2E - Evidence of transmission of minority clones in polyclonal infections - comparing to day 0 mosquito clones########

#Check if matching haplotypes in multiclonal individuals more often have high or low coverage at Day0?
test<-csp_plotme %>%
  filter(MOI > 1,
         species == "human")

# Filtering the human-only rows
human_only_df <- test %>%
  filter(Comparison == "human_only") %>%
  mutate(Percentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day)%>%
  dplyr::summarise(TotalPercentage = sum(Percentage, na.rm = TRUE))

# Filtering the matching rows
non_human_only_df <- test %>%
  filter(Comparison != "human_only") %>%
  mutate(TotalPercentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day)

# Binding the two dataframes together
filtered_df <- bind_rows(human_only_df, non_human_only_df)%>%
  filter(Day%in% c(0,2,7,14))


ggplot(filtered_df, aes(x=Comparison,y=as.numeric(TotalPercentage))) + 
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 2) +
  theme_classic()+
  xlab("Type")+
  ylab("Percentage of total reads")+
  facet_wrap(~Day,ncol=4)+
  ggtitle("CSP")
ggsave("Percentage_Reads_Human_only_vs_matching_ComparedtoDay0Mosq_CSP.pdf", width=4, height=3)



######## Figure 3A - Overall frequency of SNPs ########

# Create a complete list of combinations of sample_id and protein_change
complete_combinations <- expand.grid(sample_id = unique(combined_variants$sample_id), protein_change = unique(combined_variants$protein_change))

# Merge with original data
merged_data <- left_join(complete_combinations, combined_variants, by = c("sample_id", "protein_change"))

# remove missing positions
merged_data2 <- merged_data %>%
  anti_join(combined_missing_positions, by = c("sample_id", "protein_change"))

# replace NA in freq with zero
merged_data2$freq[is.na(merged_data2$freq)] <- 0

gene_protein_mapping <- merged_data2 %>%
  dplyr::select(protein_change, gene) %>%
  distinct() %>%
  na.omit()

merged_data2 <- merged_data2 %>%
  left_join(gene_protein_mapping, by = "protein_change", suffix = c("", "_ref"))

merged_data2 <- merged_data2 %>%
  mutate(gene = ifelse(is.na(gene), gene_ref, gene)) %>%
  dplyr::select(-gene_ref)

sampleid_species_mapping <- merged_data2 %>%
  dplyr::select(sample_id, species) %>%
  distinct() %>%
  na.omit()

merged_data2 <- merged_data2 %>%
  left_join(sampleid_species_mapping, by = "sample_id", suffix = c("", "_ref"))

merged_data2 <- merged_data2 %>%
  mutate(species = ifelse(is.na(species), species_ref, species)) %>%
  dplyr::select(-species_ref)

merged_data2 <- merged_data2 %>%
  filter(!protein_change %in% c("p.Gly102Gly","p.Ser436Phe","p.SerAla436AlaGly","p.Phe180Phe","p.Asp575Tyr","p.Ser436Tyr","p.Arg571Met","p.Val494Phe")) %>%
  mutate(protein_change = factor(recode(protein_change,
                                        "p.Lys76Thr" = "Lys76Thr",
                                        "p.Tyr184Phe" = "Tyr184Phe",
                                        "p.Asn86Tyr" = "Asn86Tyr",
                                        "p.Asn51Ile" = "Asn51Ile",
                                        "p.Cys59Arg" = "Cys59Arg",
                                        "p.Ser108Asn" = "Ser108Asn",
                                        "p.Ile431Val" = "Ile431Val",
                                        "p.Ser436Ala" = "Ser436Ala",
                                        "p.Ala437Gly" = "Ala437Gly",
                                        "p.Lys540Glu" = "Lys540Glu",
                                        "p.Ala581Gly" = "Ala581Gly",
                                        "p.Ala613Ser" = "Ala613Ser"),
                                 levels = c("Lys76Thr", "Asn86Tyr","Tyr184Phe", "Asn51Ile",
                                            "Cys59Arg", "Ser108Asn", "Ile431Val", "Ser436Ala", "Ser436Phe",
                                            "Ala437Gly", "Lys540Glu", "Ala581Gly", "Ala613Ser")))


summary_data <- merged_data2 %>%
  group_by(protein_change,species) %>%
  dplyr::mutate(count = n()) %>%
  ungroup()%>%
  group_by(protein_change,gene,species,count) %>%
  dplyr::summarise(
    mean_freq = mean(freq, na.rm = TRUE),
    se = sd(freq, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )%>%
  dplyr::mutate(gene= factor(gene,levels=c("CRT","MDR1","DHFR-TS","PPPK-DHPS","K13")))



ggplot(summary_data, aes(x = as.factor(protein_change), y = mean_freq,pattern = species,color =gene)) +
  geom_bar_pattern(aes(pattern_color=gene), fill= "white",position = "dodge", stat = "identity", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  geom_errorbar(aes(ymin = mean_freq - se, ymax = mean_freq + se,color = gene), 
                width = 0.2, position = position_dodge(0.9)) +
  #geom_text(aes(label = paste(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Protein change", y = "Overall frequency in entire sample set", fill = "Gene") +
  scale_pattern_manual(values = c("none", "stripe"),labels=c("Human", "Mosquito")) +
  scale_color_manual(values = pal5) +
  scale_pattern_fill_manual(values = pal5) +
  scale_pattern_color_manual(values = pal5)+
  theme_light()+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),limits=c(0,1))+
  theme(legend.position = "none",
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=14))
ggsave("Figure 3A_Overall_frequencies_bar.pdf",width=20, height=7)


# Fisher exact
fisher_test_func <- function(species, freq) {
  fisher.test(table(species, freq > 0), alternative = "two.sided")$p.value
}

# Perform Fisher exact test for each protein change
fisher_results <- merged_data2 %>%
  group_by(protein_change) %>%
  dplyr::summarise(fisher_p_value = fisher_test_func(species, freq), .groups = 'drop') %>%
  ungroup()

# Print or inspect the results
print(fisher_results)



######## Figure 3B - Pairwise comparison of SNPs ########

human_uniqueid<-unique(merged_data2$uniqueid[merged_data2$species=="human"])
mosq_uniqueid<-unique(merged_data2$uniqueid[merged_data2$species=="mosquito"])
Variants_Hum_matching<-merged_data2[merged_data2$species=="human" & merged_data2$uniqueid %in% mosq_uniqueid,]
Variants_Mosq_matching<-merged_data2[merged_data2$species=="mosquito" & merged_data2$uniqueid %in% human_uniqueid,]
Variants_Hum_matching$species<-"human"
Variants_Mosq_matching$species<-"mosquito"


hum_df <- Variants_Hum_matching
mosq_df <- Variants_Mosq_matching


# Initialize an empty dataframe to store the results
final_hum_df <- data.frame()

# Iterate over each unique mosqid
for (mosqid in unique_mosqids) {
  # Find the uniqueids corresponding to the current mosqid
  corresponding_uniqueids <- mosq_df$uniqueid[mosq_df$mosqid == mosqid]
  
  # Find the rows in hum_df that match these uniqueids
  matched_rows <- hum_df[hum_df$uniqueid %in% corresponding_uniqueids, ]
  
  # If there are matching rows, add the mosqid to these rows
  if (nrow(matched_rows) > 0) {
    matched_rows$mosqid <- mosqid  # Add the mosqid as a new column
    
    # Append the matched rows to the final dataframe
    final_hum_df <- rbind(final_hum_df, matched_rows)
  }
}


# Perform left joins to combine 'human' and 'mosquito' data on uniqueid and genome_pos
merged_df <- merge(Variants_Mosq_matching, final_hum_df, by = c('mosqid', 'uniqueid','genome_pos','protein_change','gene'), suffixes = c('_mosquito', '_human'), all = TRUE)

# Replace NA values with 0 in freq columns after merge
merged_df[is.na(merged_df)] <- 0

# Calculate the frequency differences
merged_df$freq_diff <- as.numeric(merged_df$freq_human - merged_df$freq_mosquito)

# Select necessary columns for the output dataframe
output_df <- merged_df %>%
  filter(!protein_change %in% c("p.Gly102Gly","p.Ser436Phe","p.SerAla436AlaGly","p.Phe180Phe","p.Asp575Tyr","p.Ser436Tyr","p.Arg571Met","p.Val494Phe")) %>%
  mutate(protein_change = factor(recode(protein_change,
                                        "p.Lys76Thr" = "Lys76Thr",
                                        "p.Tyr184Phe" = "Tyr184Phe",
                                        "p.Asn86Tyr" = "Asn86Tyr",
                                        "p.Asn51Ile" = "Asn51Ile",
                                        "p.Cys59Arg" = "Cys59Arg",
                                        "p.Ser108Asn" = "Ser108Asn",
                                        "p.Ile431Val" = "Ile431Val",
                                        "p.Ser436Ala" = "Ser436Ala",
                                        "p.Ala437Gly" = "Ala437Gly",
                                        "p.Lys540Glu" = "Lys540Glu",
                                        "p.Ala581Gly" = "Ala581Gly",
                                        "p.Ala613Ser" = "Ala613Ser"),
                                 levels = c("Lys76Thr", "Asn86Tyr","Tyr184Phe" ,"Asn51Ile",
                                            "Cys59Arg", "Ser108Asn", "Ile431Val", "Ser436Ala", "Ser436Phe",
                                            "Ala437Gly", "Lys540Glu", "Ala581Gly", "Ala613Ser"))) %>%
  dplyr::select(uniqueid, genome_pos, freq_diff, protein_change, gene, mosqid)

# Calculate mean and standard error for each protein_change
protein_summary <- output_df %>%
  filter(freq_diff != 0) %>%
  group_by(protein_change, gene) %>%
  dplyr::summarise(
    mean_freq_diff = mean(freq_diff, na.rm = TRUE),
    se = sd(freq_diff, na.rm = TRUE) / sqrt(n()),
    count = n()
  ) %>%
  mutate(gene = factor(gene, levels = c("CRT", "MDR1", "DHFR-TS", "PPPK-DHPS")))


# Create the bar plot with error bars
ggplot(protein_summary, aes(x = as.factor(protein_change), y = mean_freq_diff,fill=gene, color=gene)) +
  geom_bar(stat = "identity", position = position_dodge(),alpha=0.5) +
  geom_errorbar(aes(ymin = mean_freq_diff - se, ymax = mean_freq_diff + se), width = .2, position = position_dodge(.9)) +
  geom_text(aes(label = count, 
                y = ifelse(mean_freq_diff > 0, mean_freq_diff + se + 0.05, mean_freq_diff - se - 0.05)),
            vjust = 0, position = position_dodge(0.9)) +
  scale_color_manual(values = pal4) +
  scale_fill_manual(values=pal4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "\n Protein Change", y = "Mean Freq Diff")+
  theme_light()+
  theme(legend.position = "none",
        axis.text.x=element_text(size=18),
        axis.title.x=element_text(size=18),
        axis.text.y=element_text(size=14))
ggsave("Pairwise_comparison.pdf",width=20, height=7)



######## Figure S1 - Detection of minority clones ######## 
labstrain_ratios_data %>%
  mutate(Ratio = factor(Ratio,levels=c("1:500","1:200","1:100","1:75","1:50","1:20","1:10","1:5","1:1","1:0","0:1")),
         Haplotype = factor(Haplotype, levels =c("Background","3D7","HB3")),
         Marker = factor(Amplicon, levels =c("TRAP","CSP"))) %>%
  ggplot(aes(fill=Haplotype, y=Percentage, x=Ratio)) + 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  xlab("HB3:3D7 ratio")+
  ylab("Percentage of total reads")+
  facet_wrap(~Marker,ncol=2, strip.position ="top")+
  scale_fill_manual(values = c("#FDCDAC","aliceblue","#B3E2CD"))+
  coord_flip()

ggsave("FigureS1_Labstrain_ratios_horizontal.pdf",height=4,width=8)

######## Figure S2 - Parasite densities and mosquito infection rates over the course of ACT treatment ########
clinicaldata_gam<-clinicaldata %>%
  mutate(studycode =as.factor(studycode),
         MIR = mosq_pos/mosq_total) %>%
  filter(arm_num %in% c(1,3),
         !studyvisit==10)
mod2b = gam(totalgct_ul ~  s(studyvisit)+s(studycode, bs="re"), family=nb(), data=clinicaldata_gam)
summary(mod2b)
plot(mod2b)

newd = expand.grid(studyvisit = seq(0,max(clinicaldata_gam$studyvisit),by=0.1),studycode = clinicaldata_gam$studycode[1])
newd$prop=predict(mod2b,newd, type="response", exclude=c("s(studycode)"))
fam <- family(mod2b)
fam
str(fam)
ilink <- fam$linkinv
ilink
ilink <- family(mod2b)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2b, newd, se.fit = TRUE)[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=studyvisit, y= fit_resp))+
  geom_jitter(data=clinicaldata_gam,aes(x=studyvisit,y=totalgct_ul),alpha=0.5)+
  geom_line(col="darkgreen", linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr), fill="darkgreen", alpha=0.3)+
  theme_prism()+
  coord_cartesian(ylim=c(0,100))+
  scale_x_continuous(breaks=c(0,2,7,14,21,28,35,42,49))+
  ylab("Gametocyte density/µL (pcr)")+
  xlab("Duration of infection (days)")
ggsave("Gam_GC.pdf", width = 7, height = 6)

#Asexual density
mod2a = gam(ring_ul_all ~  s(studyvisit,k=9)+s(studycode, bs="re", k=9), family=nb(), data=clinicaldata_gam)
summary(mod2a)
plot(mod2a)

newd = expand.grid(studyvisit = seq(0,max(clinicaldata_gam$studyvisit),by=0.1),studycode = clinicaldata_gam$studycode[1])
newd$prop=predict(mod2a,newd, type="response", exclude=c("s(studycode)"))
fam <- family(mod2a)
fam
str(fam)
ilink <- fam$linkinv
ilink
ilink <- family(mod2a)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2a, newd, se.fit = TRUE)[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=studyvisit, y= fit_resp))+
  geom_jitter(data=clinicaldata_gam,aes(x=studyvisit,y=ring_ul_all),alpha=0.5)+
  geom_line(col="darkblue", linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr), fill="darkblue", alpha=0.3)+
  theme_prism()+
  coord_cartesian(ylim=c(0,1200))+
  scale_x_continuous(breaks=c(0,2,7,14,21,28,35,42,49))+
  ylab("Asexual parasite density/µL")+
  xlab("Duration of infection (days)")
ggsave("Gam_Asex.pdf", width = 7, height = 6)


#Mosquito infection rate
clinicaldata_gam_clean = clinicaldata_gam[!is.na(clinicaldata_gam$MIR), ]
mod2c = gam(MIR ~  s(studyvisit, k=8)+s(studycode, bs="re", k=8), family=nb(), data=clinicaldata_gam_clean)

summary(mod2c)
plot(mod2c)

newd = expand.grid(studyvisit = seq(0,max(clinicaldata_gam_clean$studyvisit),by=0.1),studycode = clinicaldata_gam_clean$studycode[1])

newd$prop=predict(mod2c,newd, type="response", exclude=c("s(studycode)"))
fam <- family(mod2c)
fam
str(fam)
ilink <- fam$linkinv
ilink
ilink <- family(mod2c)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2c, newd, se.fit = TRUE)[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=studyvisit, y= fit_resp))+
  geom_jitter(data=clinicaldata_gam_clean,aes(x=studyvisit,y=MIR),alpha=0.5)+
  geom_line(col="darkorange", linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr), fill="darkorange", alpha=0.3)+
  theme_prism()+
  coord_cartesian(ylim=c(0,1))+
  scale_x_continuous(breaks=c(0,2,7,14,21,28,35,42,49))+
  ylab("Mosquito infection rate")+
  xlab("Duration of infection (days)")
ggsave("Gam_MIR.pdf", width = 7, height = 6)

######## Figure S3 - Correlation between replicates and markers ########
trap_raw %>%
  filter(grepl("^trap-", Haplotype)) %>%
  group_by(SampleID) %>%
  dplyr::summarise(Count = n(), .groups = 'drop') %>%
  separate(SampleID, into = c("IndividualID", "Rep"), sep = "_Rep") %>%
  mutate(Rep = paste0("Rep", Rep)) %>%
  pivot_wider(names_from = Rep, values_from = Count, values_fill = list(Count = 0)) %>%
  mutate(Marker = "trap") %>%
  bind_rows(
    csp_raw %>%
      filter(grepl("^csp-", Haplotype)) %>%
      group_by(SampleID) %>%
      dplyr::summarise(Count = n(), .groups = 'drop') %>%
      separate(SampleID, into = c("IndividualID", "Rep"), sep = "_Rep") %>%
      mutate(Rep = paste0("Rep", Rep)) %>%
      pivot_wider(names_from = Rep, values_from = Count, values_fill = list(Count = 0)) %>%
      mutate(Marker = "csp")
  ) %>%
  { 
    # Calculate Pearson correlation coefficients 
    cor_trap <- cor.test(.$Rep1[.$Marker == "trap"], .$Rep2[.$Marker == "trap"], method = "pearson")$estimate
    cor_csp <- cor.test(.$Rep1[.$Marker == "csp"], .$Rep2[.$Marker == "csp"], method = "pearson")$estimate
    
    df_trap <- cor.test(.$Rep1[.$Marker == "trap"], .$Rep2[.$Marker == "trap"], method = "pearson")$parameter
    df_csp <- cor.test(.$Rep1[.$Marker == "csp"], .$Rep2[.$Marker == "csp"], method = "pearson")$parameter
    
    p_trap <- cor.test(.$Rep1[.$Marker == "trap"], .$Rep2[.$Marker == "trap"], method = "pearson")$p.value
    p_csp <- cor.test(.$Rep1[.$Marker == "csp"], .$Rep2[.$Marker == "csp"], method = "pearson")$p.value
    print(p_trap)
    print(p_csp)
    
    ggplot(., aes(x = Rep1, y = Rep2, color = Marker)) + 
      geom_jitter(alpha = 0.5) +
      geom_smooth(method = lm) +
      annotate("text", x = 8, y = 2, label = paste0("trap: r(",df_trap,")=", round(cor_trap,2),", p < 0.0001"), hjust = 0, color = "#4C67BD") +
      annotate("text", x = 8, y = 1, label = paste0("csp: r(",df_csp,")=", round(cor_csp,2),", p < 0.0001"), hjust = 0, color = "#D90368") +
      theme_classic() +
      xlab("Rep1 MOI") +
      ylab("Rep2 MOI") +
      scale_x_discrete(limits = 0:12) +
      scale_y_discrete(limits = 0:10) +
      ggtitle("Correlation between replicates MOI") +
      scale_color_manual(values = c("trap" = "#4C67BD", "csp" = "#D90368"))

  }
ggsave("FigS2_Correlation_Replicates.pdf", width = 7, height = 6)

csp_trap_finalhaplotypes %>%
  {
    # Calculate Pearson correlation coefficients within the pipe using curly braces
    cor <- cor.test(.$csp_MOI, .$trap_MOI, method = "pearson")$estimate
    df <- cor.test(.$csp_MOI, .$trap_MOI, method = "pearson")$parameter
    p <- cor.test(.$csp_MOI, .$trap_MOI, method = "pearson")$p.value
    
    print(p)
    
    ggplot(csp_trap_finalhaplotypes, aes(x=csp_MOI, y=trap_MOI)) + 
      geom_jitter()+
      theme_classic()+
      theme(legend.position = "none")+
      annotate("text", x = 6, y = 5, label = paste0("r(",df,")=", round(cor,2),", p < 0.0001"), hjust = 0) +
      xlab("csp MOI")+
      ylab("trap MOI")+
      scale_x_discrete(limits = 0:12) +
      scale_y_discrete(limits = 0:10) +
      geom_smooth(method=lm)+
      ggtitle("Correlation between csp and trap MOI")
  }
ggsave("FigS2_Correlation_csp_trap.pdf", width=7, height=6)


######## Figure S4 - Haplotype Networks ########

# Convert raw byte vectors to character strings
haplotype_sequences_TRAP <- sapply(TRAP_haplotypes, function(x) {
  # Convert DNAString to character
  as.character(x)
})
# Convert raw byte vectors to character strings
haplotype_sequences_CSP <- sapply(CSP_haplotypes, function(x) {
  # Convert DNAString to character
  as.character(x)
})


# Iterate over the dataframe and write sequences to the FASTA file
for (i in 1:nrow(trap_matching_haplotypes_molten)) {
  individual <- trap_matching_haplotypes_molten$Individual[i]
  haplotype <- trap_matching_haplotypes_molten$Haplotype[i]
  species <- trap_matching_haplotypes_molten$species[i]
  day <- trap_matching_haplotypes_molten$Day[i]
  
  # Check if the haplotype is in the list
  if (haplotype %in% names(haplotype_sequences_TRAP)) {
    sequence <- haplotype_sequences_TRAP[[haplotype]]
    sequence <- as.character(sequence) # Make sure to convert to character
  } else {
    sequence <- "Sequence not found"
  }
  
  # Write the sequence to the FASTA file
  write_haplotype_sequences(individual, haplotype, day, species, sequence, 'TRAP_plotme.fasta')
}

# Iterate over the dataframe and write sequences to the FASTA file
for (i in 1:nrow(csp_matching_haplotypes_molten)) {
  individual <- csp_matching_haplotypes_molten$Individual[i]
  haplotype <- csp_matching_haplotypes_molten$Haplotype[i]
  species <- csp_matching_haplotypes_molten$species[i]
  day <- csp_matching_haplotypes_molten$Day[i]
  
  # Check if the haplotype is in the list
  if (haplotype %in% names(haplotype_sequences_CSP)) {
    sequence <- haplotype_sequences_CSP[[haplotype]]
    sequence <- as.character(sequence) # Make sure to convert to character
  } else {
    sequence <- "Sequence not found"
  }
  
  # Write the sequence to the FASTA file
  write_haplotype_sequences(individual, haplotype, day,species, sequence, 'CSP_plotme.fasta')
}


seqs <- read.FASTA("CSP_plotme.fasta")
seqs <- read.FASTA("TRAP_plotme.fasta")
haps <- haplotype(seqs)
dist <- dist.dna(haps, "N")
net <- rmst(dist, B=100)
(sz <- summary(haps))
nt.labs <- names(seqs)

matrix <- as.matrix(seqs)
rownames(matrix) <- nt.labs
regions <- haploFreq(matrix,split="_",what=3)

#Match haplotypR haplotype names to roman numbers
first_numbers <- sapply(attributes(haps)$index, function(x) x[1])
selected_labels <- nt.labs[first_numbers]
numbers_csp <- str_extract(selected_labels, "csp-[0-9]+")
numbers_trap <- str_extract(selected_labels, "trap-[0-9]+")
roman_numbers<-rownames(as.data.frame(regions))
final_df_csp <- data.frame(RomanNumeral = roman_numbers, Haplotype = numbers_csp)
final_df_trap <- data.frame(RomanNumeral = roman_numbers, Haplotype = numbers_trap)

#col=usecol(pal_unikn_pref)

plot(net, size=sz, scale.ratio =2, cex=0.75, pie = regions, bg = pal6,threshold = 0, col.link="black")

replot()
legend("left", colnames(regions),col=pal6, pch=20, cex=0.7)

current_plot <- recordPlot()
pdf("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork/TRAP_network.pdf",width = 12,height = 10)
pdf("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork/legend.pdf",width = 8,height = 4)
replayPlot(current_plot)
dev.off()






######## Figure S5 - Percentage of polyclonal infections ########

csp_trap_finalhaplotypes %>%
  mutate(
    trap_MonoPoly = ifelse(is.na(trap_MOI), NA, ifelse(trap_MOI == 1, "Mono", "Poly")),
    csp_MonoPoly = ifelse(is.na(csp_MOI), NA, ifelse(csp_MOI == 1, "Mono", "Poly")),
    MonoPoly = ifelse(
      is.na(csp_MonoPoly) & is.na(trap_MonoPoly), NA, 
      ifelse(is.na(csp_MonoPoly), trap_MonoPoly,
             ifelse(is.na(trap_MonoPoly), csp_MonoPoly,
                    ifelse(csp_MonoPoly == "Poly" | trap_MonoPoly == "Poly", "Poly", "Mono")
             )
      )
    ),
    Species = factor(Species, levels = c("mosquito", "human"))
  ) %>%
  select(-Mosquito) %>%
  {
    counts <- . %>%
      group_by(Day, Species, MonoPoly) %>%
      summarize(Count = n(), .groups = 'drop') %>%
      group_by(Day, Species) %>%
      mutate(LabelPos = cumsum(Count) - (0.5 * Count))
    
    ggplot(.[.$Day %in% c("0", "2", "7", "14"),], aes(x = as.factor(Day), fill = factor(MonoPoly))) + 
      geom_bar(position = "fill", stat = "count") +
      geom_text(data = counts, aes(x = as.factor(Day), y = LabelPos/sum(Count), label = Count), position = position_fill(), vjust = 1.5) +
      xlab("Days after treatment initiation") +
      ylab("Percentage of samples") +
      facet_wrap(~Species, ncol = 1) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
      scale_fill_manual(values = pal2)
  }
ggsave("MonoPoly_MarkersCombined.pdf", width=6, height=6)



######## Figure S6 - Median MOI at all timepoints in both species ########

a <- csp_trap_finalhaplotypes %>%
  filter(Species == "human", Day %in% c("0", "2", "7", "14", "21", "28")) %>%
  mutate(Day = as.factor(Day)) %>%
  group_by(Day) %>%
  dplyr::summarise(
    Median_MOI = median(MOI_Combined, na.rm = TRUE),
    IQR_Lower = quantile(MOI_Combined, probs = 0.25, na.rm = TRUE),
    IQR_Upper = quantile(MOI_Combined, probs = 0.75, na.rm = TRUE),
    MOI_Combined = list(MOI_Combined)
  ) %>%
  unnest(MOI_Combined) %>%
  ggplot(aes(x = Day, y = MOI_Combined)) +
  geom_violin(fill = "grey", alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  geom_point(aes(y = Median_MOI), color = "#DF7070", size = 2) + 
  geom_line(aes(y = Median_MOI, group = 1), color = "#DF7070", size = 1) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x = "Days after treatment initiation", y = "MOI", title = "MOI in human samples")


b<-csp_trap_finalhaplotypes %>%
  filter(Species == "mosquito", Day %in% c("0", "2", "7", "14", "21", "28")) %>%
  mutate(Day = as.factor(Day)) %>%
  group_by(Day) %>%
  dplyr::summarise(
    Median_MOI = median(MOI_Combined, na.rm = TRUE),
    IQR_Lower = quantile(MOI_Combined, probs = 0.25, na.rm = TRUE),
    IQR_Upper = quantile(MOI_Combined, probs = 0.75, na.rm = TRUE),
    MOI_Combined = list(MOI_Combined)
  ) %>%
  unnest(MOI_Combined) %>%
  ggplot(aes(x = Day, y = MOI_Combined)) +
  geom_violin(fill = "grey", alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  geom_point(aes(y = Median_MOI), color = "#DF7070", size = 2) +
  geom_line(aes(y = Median_MOI, group = 1), color = "#DF7070", size = 1) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x = "Days after treatment initiation", y = "MOI", title = "MOI in mosquito samples")


combined_clonality_plot <- b + a + plot_layout(ncol=1)
ggsave("Median_MOI_species_timepoints.pdf", combined_clonality_plot, width=7, height=10)





######## Figure S7 - Gametocyte densities, mosquito infection rates and oocyst densities at each MOI ########
a<-ggplot(merge_clinical_haplotypes,aes(x=as.factor(MOI_Combined), y=totalgct_ul)) + 
  geom_boxplot(alpha = 0.5, fill="grey")+
  theme_classic()+
  #geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI human blood",y="Gametocytes per µL")+
  scale_y_log10()+
  annotation_logticks(sides = "l")
  

b<-ggplot(merge_clinical_haplotypes, aes(x=as.factor(MOI_Combined), y=percentagemosqinfected)) + 
  geom_boxplot(alpha = 0.5,fill="grey")+
  theme_classic()+
  #geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI human blood",y="Mosquito infection rate (%)")

c<-ggplot(merge_oocyst_individual_haplotypes, aes(x=as.factor(MOI_Combined), y=oocysts)) + 
  geom_boxplot(alpha = 0.5, fill="grey")+
  theme_classic()+
  #geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI human blood",y="Oocyst density")

FigureS6<-a+b+c + plot_layout(ncol=1)
ggsave("FigureS6.pdf", FigureS6, width=5, height=12)

######## Figure S8 - Baseline MOI by age group and month ########
x<-merge_clinical_haplotypes %>%
  mutate(
    age_group = case_when(
      age >= 5 & age < 10 ~ "5-10",
      age >= 10 & age < 15 ~ "10-15",
      age >= 15 ~ ">15"
    ),
    age_group = factor(age_group, levels = c("5-10", "10-15", ">15"))
  ) %>%
  filter(Day == "0") %>%
  ggplot(aes(x = age_group, y = MOI_Combined)) +
  geom_boxplot(fill="grey")+
  geom_jitter(width=0.1)+
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Age group", y = "MOI") +
  scale_y_continuous(breaks = c(1, 2, 3,4,5,6,7,8,9)) 


# Perform pairwise Wilcoxon Rank Sum test
prepared_data <- merge_clinical_haplotypes %>%
  mutate(
    age_group = case_when(
      age >= 5 & age < 10 ~ "5-10",
      age >= 10 & age < 15 ~ "10-15",
      age >= 15 ~ ">15"
    )
  ) %>%
  mutate(age_group = factor(age_group, levels = c("5-10", "10-15", ">15")))

pairwise.wilcox.test(prepared_data$MOI_Combined, prepared_data$age_group, p.adjust.method = "bonf",paired=FALSE)



y<-merge_clinical_haplotypes %>%
  mutate(
    visitdate_hb = as.Date(visitdate_hb, format="%d/%m/%Y"),
    month = format(visitdate_hb, "%m"),
    month_group = case_when(
      month =="09" ~ "Sep",
      month ==10 ~ "Oct",
      month ==11 ~ "Nov",
      month ==12 ~ "Dec",
    ),
    month_group = factor(month_group, levels = c("Sep", "Oct", "Nov","Dec"))
  ) %>%
  filter(month_group %in% c("Sep", "Oct", "Nov")) %>%
  filter(Day == "0") %>%
  ggplot(aes(x=month_group, y=MOI_Combined)) + 
  geom_boxplot(fill="grey")+
  geom_jitter(width=0.1)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = c(1, 2, 3,4,5,6,7,8,9)) +
  labs(x= "Month",y="MOI")


FigureS7<-x+y+plot_layout(ncol=1)
ggsave("FigureS7.pdf", FigureS7, width=6, height=9)
######## Figure S? - Pftrap Evidence of transmission of minority clones ######## 
#Check if matching haplotypes in multiclonal individuals more often have high or low coverage at Day0?
test<-trap_merge_clinical_matching_haplotypes_molten %>%
  filter(MOI > 1,
         species == "human")

# Filtering the human-only rows
human_only_df <- test %>%
  filter(Comparison == "human_only") %>%
  mutate(Percentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day,.groups='keep')%>%
  dplyr::summarise(TotalPercentage = sum(Percentage, na.rm = TRUE))

# Filtering the matching rows
non_human_only_df <- test %>%
  filter(Comparison != "human_only") %>%
  mutate(TotalPercentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day)

# Binding the two dataframes together
filtered_df <- bind_rows(human_only_df, non_human_only_df) %>%
  filter(Day%in% c(0,2,7,14))


ggplot(filtered_df, aes(x=Comparison,y=as.numeric(TotalPercentage))) + 
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 2) +
  theme_classic()+
  xlab("Type")+
  ylab("Percentage of total reads")+
  facet_wrap(~Day,ncol=4)+
  ggtitle("trap")
ggsave("Percentage_Reads_Human_only_vs_matching_trap.pdf", width=4, height=3)



######## Figure S9 - Haplotype count in human and mosquito hosts ########
csp_grouped_df <- csp_matching_haplotypes_molten %>%
  group_by(Haplotype, Comparison) %>%
  dplyr::summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
csp_grouped_df <- csp_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("csp-", "", Haplotype))) %>%
  arrange(Haplotype_num)

trap_grouped_df <- trap_matching_haplotypes_molten %>%
  group_by(Haplotype, Comparison) %>%
  dplyr::summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
trap_grouped_df <- trap_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("trap-", "", Haplotype))) %>%
  arrange(Haplotype_num)

# Plot
e<-ggplot(csp_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal1)+
  labs(title = "csp",
       x = "Haplotype",
       y = "Count") +
  theme_light()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
f<-ggplot(trap_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal1)+
  labs(title = "trap",
       x = "Haplotype",
       y = "Count") +
  theme_light()+
  theme(axis.text.x = element_blank())
combined_plot <- e / f 
ggsave("Haplotypes_transmission.pdf", combined_plot, width=5, height=7)


#Check if matching haplotypes in multiclonal individuals more often have high or low coverage at Day0?
test<-trap_plotme %>%
  filter(MOI > 2,
         species == "human")

# Filtering the human-only rows
human_only_df <- test %>%
  filter(Comparison == "human_only") %>%
  mutate(Percentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day,.groups='keep')%>%
  dplyr::summarise(TotalPercentage = sum(Percentage, na.rm = TRUE))

# Filtering the matching rows
non_human_only_df <- test %>%
  filter(Comparison != "human_only") %>%
  mutate(TotalPercentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison, Day)

# Binding the two dataframes together
filtered_df <- bind_rows(human_only_df, non_human_only_df)%>%
  filter(Day%in% c(0,2,7,14))


ggplot(filtered_df, aes(x=Comparison,y=as.numeric(TotalPercentage))) + 
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 2) +
  theme_classic()+
  xlab("Type")+
  ylab("Percentage of total reads")+
  facet_wrap(~Day,ncol=4)+
  ggtitle("trap")
ggsave("Percentage_Reads_Human_only_vs_matching_ComparedtoDay0Mosq_trap.pdf", width=4, height=3)














######## Figure S10 - Odds of transmission for each haplotype ########
#Pfcsp
Clin_csp <- csp_haplotypes_molten %>%
  left_join(dplyr::select(clinicaldata, SampleID, totalgct_ul, ring_ul_all), by = "SampleID")

# Normalize read percentages in human data
human_data_csp <- Clin_csp %>% 
  filter(species != 'mosquito') %>%
  mutate(Normalized_csp_Perc = as.numeric(Percentage) / as.numeric(totalgct_ul))

# Prepare mosquito data
mosquito_data_csp <- Clin_csp %>% 
  filter(species == 'mosquito')

# Merge human and mosquito data
merged_data_csp <- merge(human_data_csp, mosquito_data_csp, by = c("Timepoint", "Haplotype"), all = TRUE) %>%
  mutate(Normalized_csp_Perc = ifelse(is.na(Normalized_csp_Perc), 0, Normalized_csp_Perc),
         Percentage.y = ifelse(is.na(Percentage.y), 0, Percentage.y),
         Higher_in_Mosquito = Percentage.y > Normalized_csp_Perc)

# Calculate the odds of transmission for each haplotype
haplotype_odds <- merged_data_csp %>%
  group_by(Haplotype) %>%
  summarize(Odds_of_Transmission = mean(Higher_in_Mosquito, na.rm = TRUE), Count = n(), .groups = 'drop') %>%
  mutate(Haplotype_Num = as.numeric(gsub("[^0-9]", "", Haplotype))) %>%
  arrange(Haplotype_Num) %>%
  dplyr::select(-Haplotype_Num) %>%
  mutate(Haplotype = factor(Haplotype, levels = unique(Haplotype))) %>%
  filter(Count>3) %>%
  left_join(final_df_csp,by="Haplotype") %>%
  mutate(DisplayRoman = ifelse(Odds_of_Transmission > 0.2, as.character(RomanNumeral), NA))

y_axis_data <- data.frame(Y = seq(0, 1, by = 0.1), Odds_of_Transmission = 0)

# Create the circular plot with the y-axis circles
a<-ggplot(haplotype_odds, aes(x = Haplotype, y = Odds_of_Transmission, fill = Odds_of_Transmission)) +
  geom_hline(yintercept = seq(0, 1, by = 0.2), color = "grey", linetype = "dashed")+
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = DisplayRoman), position = position_stack(vjust = 1.1), color = "black", size = 3) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.8), color = "white", size = 2) +
  scale_fill_gradient2(low = "white", high = "lightblue4", limits = c(0, 1)) +
  coord_polar(start = 0) +
  theme_minimal() +
  labs(
    title = "Circular Visualization of Transmission Odds for Each Haplotype normalized for totalGCT",
    x = NULL,
    y = "Odds of Transmission"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )



#Pftrap

Clin_trap <- trap_haplotypes_molten %>%
  left_join(dplyr::select(clinicaldata, SampleID, totalgct_ul, ring_ul_all), by = "SampleID")

# Normalize read percentages in human data
human_data_trap <- Clin_trap %>% 
  filter(species != 'mosquito') %>%
  mutate(Normalized_trap_Perc = as.numeric(Percentage) / as.numeric(totalgct_ul))

# Prepare mosquito data
mosquito_data_trap <- Clin_trap %>% 
  filter(species == 'mosquito')

# Merge human and mosquito data
merged_data_trap <- merge(human_data_trap, mosquito_data_trap, by = c("Timepoint", "Haplotype"), all = TRUE) %>%
  mutate(Normalized_trap_Perc = ifelse(is.na(Normalized_trap_Perc), 0, Normalized_trap_Perc),
         Percentage.y = ifelse(is.na(Percentage.y), 0, Percentage.y),
         Higher_in_Mosquito = Percentage.y > Normalized_trap_Perc)

# Calculate the odds of transmission for each haplotype
haplotype_odds <- merged_data_trap %>%
  group_by(Haplotype) %>%
  dplyr::summarize(Odds_of_Transmission = mean(Higher_in_Mosquito, na.rm = TRUE), Count = n(), .groups = 'drop') %>%
  mutate(Haplotype_Num = as.numeric(gsub("[^0-9]", "", Haplotype))) %>%
  arrange(Haplotype_Num) %>%
  dplyr::select(-Haplotype_Num) %>%
  mutate(Haplotype = factor(Haplotype, levels = unique(Haplotype))) %>%
  filter(Count>3)%>%
  left_join(final_df_trap,by="Haplotype") %>%
  mutate(DisplayRoman = ifelse(Odds_of_Transmission > 0.2, as.character(RomanNumeral), NA))

y_axis_data <- data.frame(Y = seq(0, 1, by = 0.1), Odds_of_Transmission = 0)

# Create the circular plot with the y-axis circles
b<-ggplot(haplotype_odds, aes(x = Haplotype, y = Odds_of_Transmission, fill = Odds_of_Transmission)) +
  geom_hline(yintercept = seq(0, 1, by = 0.2), color = "grey", linetype = "dashed")+
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = DisplayRoman), position = position_stack(vjust = 1.1), color = "black", size = 3) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.8), color = "white", size = 2) +
  scale_fill_gradient2(low = "white", high = "lightblue4", limits = c(0, 1)) +
  coord_polar(start = 0) +
  theme_minimal() +
  labs(
    title = "Circular Visualization of Transmission Odds for Each Haplotype normalized for totalGCT",
    x = NULL,
    y = "Odds of Transmission"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )


FigureS9<-a+b
ggsave("FigureS9_OddsOfTransmisson.pdf", width=15, height=7)






######## Investigate not-infectious samples at day 0 ########
###### NEEDS WORK

plotme_1<-csp_merge_clinical_haplotypes_molten[csp_merge_clinical_haplotypes_molten$transmission==1 & csp_merge_clinical_haplotypes_molten$studyvisit==0,]
haplotype_counts_infectious <- table(plotme_1$Haplotype)

plotme_0<-csp_merge_clinical_haplotypes_molten[csp_merge_clinical_haplotypes_molten$transmission==0 & csp_merge_clinical_haplotypes_molten$studyvisit==0,]
haplotype_counts_noninfectious <- table(plotme_0$Haplotype)

# Convert the tables to data frames
df_infectious <- as.data.frame(haplotype_counts_infectious, stringsAsFactors = FALSE)
df_noninfectious <- as.data.frame(haplotype_counts_noninfectious, stringsAsFactors = FALSE)

# Rename columns
colnames(df_infectious) <- c("Haplotype", "Count")
colnames(df_noninfectious) <- c("Haplotype", "Count")

# Add a column to indicate infectious status
df_infectious$Infectious <- "Yes"
df_noninfectious$Infectious <- "No"

# Combine the data frames
combined_df <- rbind(df_infectious, df_noninfectious)

ggplot(combined_df, aes(x = Haplotype, y = Count, fill = Infectious)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Yes" = "red", "No" = "blue")) +
  labs(title = "Haplotype Counts by Infectious Status",
       x = "Haplotype",
       y = "Counts",
       fill = "Infectious Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Improve x-label readability

pie(haplotype_counts, labels = paste(names(haplotype_counts), "(", haplotype_counts, ")", sep = ""), main = "Haplotype Distribution")
pie(haplotype_counts, labels = paste(names(haplotype_counts), "(", haplotype_counts, ")", sep = ""), main = "Haplotype Distribution")


######## Check relation between midgut clonality and oocyst number (mosquitoes only) ########
ggplot(merge_oocyst_haplotypes, aes(x=as.factor(MOI_Combined), y=oocysts)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA,fill="grey")+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(1,200))+
  labs(x= "MOI midgut",y="Oocyst number")+
  ggtitle("Midgut clonality ~ oocyst number")
ggsave("MidgutClonality_oocystnumber.pdf", width=6, height=5)

#Monoclonal versus polyclonal
merge_oocyst_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  ggplot(aes(x=as.factor(MonoPoly), y=oocysts)) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA,fill="grey")+
    theme_classic()+
    geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
    theme(legend.position = "none")+
    scale_y_continuous(limits=c(1,200))+
    labs(x= "MOI midgut",y="Oocyst number")+
    ggtitle("Midgut clonality ~ oocyst number")
ggsave("MidgutMonoPoly_oocystnumber.pdf", width=6, height=5)




######## Check whether multiclonal samples have higher likelihood of causing multiclonal mosquito infections ########

inner_join(csp_trap_finalhaplotypes_human, csp_trap_finalhaplotypes_mosq, by = "Timepoint", suffix = c("_human", "_mosquito")) %>%
  {
    # Calculate Pearson correlation coefficients within the pipe using curly braces
    cor_result <- cor.test(.$MOI_Combined_human, .$MOI_Combined_mosquito, method = "pearson")
    cor_label <- sprintf("Pearson r: %.2f", cor_result$estimate)
    
    ggplot(., aes(x = MOI_Combined_human, y = MOI_Combined_mosquito)) + 
      geom_jitter() +
      labs(
        title = "Comparison of MOI_combined in Humans and Mosquitoes",
        x = "MOI_combined in Humans",
        y = "MOI_combined in Mosquitoes"
      ) +
      annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1, vjust = 1) +
      geom_smooth(method=lm)+
      scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
      scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
      theme_classic()
  }

ggsave("Clonality_human_vs_mosq.pdf", width=6, height=5)







######## Table S3. Coverage across markers ######## 
trap_raw %>%
  mutate(SampleID = gsub("(.*)_.*", "\\1", SampleID),
         Reads = as.numeric(Reads)) %>%
  filter(str_starts(Haplotype, "trap")) %>%
  group_by(SampleID) %>%
  dplyr::summarise(
    TotalReads = sum(Reads, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::summarise(
    MedianReads = median(TotalReads, na.rm = TRUE),
    Q1Reads = quantile(TotalReads, probs = 0.25, na.rm = TRUE),
    Q3Reads = quantile(TotalReads, probs = 0.75, na.rm = TRUE)
  )
  
csp_raw %>%
  mutate(SampleID = gsub("(.*)_.*", "\\1", SampleID),
         Reads = as.numeric(Reads)) %>%
  filter(str_starts(Haplotype, "csp")) %>%
  group_by(SampleID) %>%
  dplyr::summarise(
    TotalReads = sum(Reads, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::summarise(
    MedianReads = median(TotalReads, na.rm = TRUE),
    Q1Reads = quantile(TotalReads, probs = 0.25, na.rm = TRUE),
    Q3Reads = quantile(TotalReads, probs = 0.75, na.rm = TRUE)
  )

combined_coverage %>%
  mutate(Region = factor(Region,levels=c("crt","mdr1","dhfr","dhps","k13"))) %>%
  group_by(Region) %>%
  dplyr::summarise(
    Median = median(Median_Depth, na.rm = TRUE),
    IQR_Lower = quantile(Median_Depth, probs = 0.25, na.rm = TRUE),
    IQR_Upper = quantile(Median_Depth, probs = 0.75, na.rm = TRUE)
  )


######## Table S4. Non-synonymous single nucleotide polymorphisms in genes associated with drug resistance  ######## 
# Create a complete list of combinations of sample_id and protein_change
complete_combinations <- expand.grid(sample_id = unique(combined_variants$sample_id), protein_change = unique(combined_variants$protein_change))

# Merge with original data
merged_data <- left_join(complete_combinations, combined_variants, by = c("sample_id", "protein_change"))

# remove missing positions
merged_data2 <- merged_data %>%
  anti_join(combined_missing_positions, by = c("sample_id", "protein_change"))

# replace NA in freq with zero
merged_data2$freq[is.na(merged_data2$freq)] <- 0

gene_protein_mapping <- merged_data2 %>%
  dplyr::select(protein_change, gene) %>%
  distinct() %>%
  na.omit()

merged_data2 <- merged_data2 %>%
  left_join(gene_protein_mapping, by = "protein_change", suffix = c("", "_ref"))

merged_data2 <- merged_data2 %>%
  mutate(gene = ifelse(is.na(gene), gene_ref, gene)) %>%
  dplyr::select(-gene_ref)

sampleid_species_mapping <- merged_data2 %>%
  dplyr::select(sample_id, species) %>%
  distinct() %>%
  na.omit()

merged_data2 <- merged_data2 %>%
  left_join(sampleid_species_mapping, by = "sample_id", suffix = c("", "_ref"))

merged_data2 <- merged_data2 %>%
  mutate(species = ifelse(is.na(species), species_ref, species)) %>%
  dplyr::select(-species_ref)

merged_data2 <- merged_data2 %>%
  filter(!protein_change %in% c("p.Phe180Phe","p.Gly102Gly")) %>%
  mutate(protein_change = factor(recode(protein_change,
                                        "p.Lys76Thr" = "Lys76Thr",
                                        "p.Tyr184Phe" = "Tyr184Phe",
                                        "p.Asn86Tyr" = "Asn86Tyr",
                                        "p.Asn51Ile" = "Asn51Ile",
                                        "p.Cys59Arg" = "Cys59Arg",
                                        "p.Ser108Asn" = "Ser108Asn",
                                        "p.Ile431Val" = "Ile431Val",
                                        "p.Ser436Ala" = "Ser436Ala",
                                        "p.Ser436Phe" = "Ser436Phe",
                                        "p.SerAla436AlaGly" = "SerAla436AlaGly",
                                        "p.Ala437Gly" = "Ala437Gly",
                                        "p.Lys540Glu" = "Lys540Glu",
                                        "p.Ala581Gly" = "Ala581Gly",
                                        "p.Ala613Ser" = "Ala613Ser",
                                        "p.Asp575Tyr" = "Asp575Tyr",
                                        "p.Ser436Tyr" = "Ser436Tyr",
                                        "p.Arg571Met" = "Arg571Met",
                                        "p.Val494Phe" = "Val494Phe"),
                                 levels = c("Lys76Thr", "Asn86Tyr","Tyr184Phe" ,"Asn51Ile",
                                            "Cys59Arg", "Ser108Asn", "Ile431Val", "Ser436Ala", "Ser436Phe", "SerAla436AlaGly", "Ser436Tyr",
                                            "Ala437Gly", "Lys540Glu", "Arg571Met", "Asp575Tyr","Ala581Gly", "Ala613Ser","Val494Phe")))


summary_data <- merged_data2 %>%
  group_by(protein_change,species) %>%
  dplyr::mutate(count = n()) %>%
  ungroup()%>%
  group_by(protein_change,gene,species,count) %>%
  dplyr::summarise(
    mean_freq = mean(freq, na.rm = TRUE),
    se = sd(freq, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )%>%
  dplyr::mutate(gene= factor(gene,levels=c("CRT","MDR1","DHFR-TS","PPPK-DHPS","K13")))

TableS4 <- summary_data %>%
  dplyr::select(-count) %>%
  pivot_wider(names_from = species,
              values_from = c(mean_freq, se),
              names_glue = "{.value}_{species}",
              values_fn = list(mean_freq = mean, se = mean))%>%
  mutate(mean_freq_human = round(mean_freq_human * 100,2),
         mean_freq_mosquito = round(mean_freq_mosquito * 100,2))

write.csv(TableS4, "TableS4_MAF.csv")

######## Overview of all haplotypes ########
Haplotypes_trap_Molten <- trap_haplotypes_molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Individual = str_extract(Timepoint, "^[^_]+"),
         Reads = as.numeric(Reads),
         Percentage = as.numeric(Percentage),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)),
         Day_Species = paste0(Day,"_",species),
         Day_Species = factor(Day_Species, levels = c("0_human","0_mosquito","2_human","2_mosquito","7_human","7_mosquito","14_human","14_mosquito","21_human","21_mosquito","28_human","28_mosquito","35_human","35_mosquito")),
         Day_Species = as.character(Day_Species),
         Day_Species = ifelse(species == "mosquito",paste0(Day_Species, "_",sub(".*_", "", SampleID)),Day_Species)) %>%
  ungroup()

generate_ordered_levels <- function(values) {
  # Extract unique days from values
  unique_days <- unique(as.numeric(gsub("\\D", "", values)))
  
  # Sort unique days
  sorted_days <- sort(unique_days)
  
  # For each day, order "human" first followed by "mosquito"
  ordered_levels <- unlist(lapply(sorted_days, function(day) {
    human_label <- paste0(day, "_human")
    mosquito_labels <- sort(unique(grep(paste0("^", day, "_mosquito"), values, value = TRUE)))
    c(human_label, mosquito_labels)
  }))
  
  # Ensure no duplicates
  ordered_levels <- unique(ordered_levels)
  
  return(ordered_levels)
}

Haplotypes_trap_Molten$Day_Species <- factor(Haplotypes_trap_Molten$Day_Species, levels = generate_ordered_levels(Haplotypes_trap_Molten$Day_Species))

p <- ggplot(Haplotypes_trap_Molten, aes(x=Day_Species, y=Percentage, fill=Haplotype)) + 
  geom_bar(position="stack", stat="identity") +
  geom_bar_pattern(aes(pattern = species, pattern_density = ifelse(species == "human", 0, 0.1),pattern_fill = Haplotype), 
                   position = "stack", 
                   stat = "identity",
                   pattern="stripe",
                   pattern_spacing = 0.2,
                   pattern_size=0.05,
                   pattern_frequency = 0.05,
                   pattern_angle = 45) +
  facet_wrap(~Individual, scales = "free", ncol=5) +
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("trap_overview.pdf", p, width=40, height=50,limitsize = FALSE)

# Overview of all haplotypes
Haplotypes_csp_Molten <- csp_haplotypes_molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Individual = str_extract(Timepoint, "^[^_]+"),
         Reads = as.numeric(Reads),
         Percentage = as.numeric(Percentage),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)),
         Day_Species = paste0(Day,"_",species),
         Day_Species = factor(Day_Species, levels = c("0_human","0_mosquito","2_human","2_mosquito","7_human","7_mosquito","14_human","14_mosquito","21_human","21_mosquito","28_human","28_mosquito","35_human","35_mosquito")),
         Day_Species = as.character(Day_Species),
         Day_Species = ifelse(species == "mosquito",paste0(Day_Species, "_",sub(".*_", "", SampleID)),Day_Species)) %>%
  ungroup()

generate_ordered_levels <- function(values) {
  # Extract unique days from values
  unique_days <- unique(as.numeric(gsub("\\D", "", values)))
  
  # Sort unique days
  sorted_days <- sort(unique_days)
  
  # For each day, order "human" first followed by "mosquito"
  ordered_levels <- unlist(lapply(sorted_days, function(day) {
    human_label <- paste0(day, "_human")
    mosquito_labels <- sort(unique(grep(paste0("^", day, "_mosquito"), values, value = TRUE)))
    c(human_label, mosquito_labels)
  }))
  
  # Ensure no duplicates
  ordered_levels <- unique(ordered_levels)
  
  return(ordered_levels)
}

Haplotypes_csp_Molten$Day_Species <- factor(Haplotypes_csp_Molten$Day_Species, levels = generate_ordered_levels(Haplotypes_csp_Molten$Day_Species))

p <- ggplot(Haplotypes_csp_Molten, aes(x=Day_Species, y=Percentage, fill=Haplotype)) + 
  geom_bar(position="stack", stat="identity") +
  geom_bar_pattern(aes(pattern = species, pattern_density = ifelse(species == "human", 0, 0.1),pattern_fill = Haplotype), 
                   position = "stack", 
                   stat = "identity",
                   pattern="stripe",
                   pattern_spacing = 0.2,
                   pattern_size=0.05,
                   pattern_frequency = 0.05,
                   pattern_angle = 45) +
  facet_wrap(~Individual, scales = "free", ncol=5) +
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("csp_overview.pdf", p, width=40, height=50,limitsize = FALSE)




######## Calculate numbers of samples ########
table_numbers <- data.frame(
  Description = c("Total samples with csp_hap",
                  "Total samples with trap_hap",
                  "Unique individuals with csp_hap",
                  "Unique individuals with trap_hap",
                  "Timepoints with csp_hap",
                  "Timepoints with trap_hap",
                  "Individuals non-infectious at day 0",
                  "Mosquito samples with csp_hap",
                  "Mosquito samples with trap_hap",
                  "Individuals for which there are matching mosquitoes with csp_hap (any timepoint)",
                  "Individuals for which there are matching mosquitoes with trap_hap (any timepoint)",
                  "Timepoints at which there are matching human-mosquito samples with csp_hap",
                  "Timepoints at which there are matching human-mosquito samples with trap_hap",
                  "Matching samples human-mosquito with csp_hap",
                  "Matching samples human-mosquito with trap_hap",
                  "Matching samples human-mosquito with csp_hap at day14",
                  "Individuals in human drug resistance data",
                  "Individuals in mosquito drug resistance data",
                  "Individuals for which there are matching mosquitoes in drug resistance data",
                  "Mosquitoes in mosquito drug resistance data",
                  "Instances of matching individuals and mosquitoes"),
  Count = c(length(unique(csp_processed$SampleID)), length(unique(trap_processed$SampleID)),
            length(unique(csp_trap_finalhaplotypes_human$Individual[!is.na(csp_trap_finalhaplotypes_human$csp_hap)])), length(unique(csp_trap_finalhaplotypes_human$Individual[!is.na(csp_trap_finalhaplotypes_human$trap_hap)])),
            length(na.omit(csp_trap_finalhaplotypes_human$csp_MOI)), length(na.omit(csp_trap_finalhaplotypes_human$trap_MOI)),
            sum(merge_clinical_haplotypes$percentagemosqinfected == 0 & merge_clinical_haplotypes$Day == 0, na.rm = TRUE),
            length(na.omit(csp_trap_finalhaplotypes_mosq$csp_MOI)),length(na.omit(csp_trap_finalhaplotypes_mosq$trap_MOI)),
            length(unique(csp_trap_finalhaplotypes_mosq$Individual[!is.na(csp_trap_finalhaplotypes_mosq$csp_hap)])),length(unique(csp_trap_finalhaplotypes_mosq$Individual[!is.na(csp_trap_finalhaplotypes_mosq$trap_hap)])),
            length(unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$csp_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$csp_hap)]))),
            length(unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$trap_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$trap_hap)]))),
            nrow(csp_trap_finalhaplotypes_mosq[csp_trap_finalhaplotypes_mosq$Timepoint %in% unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$csp_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$csp_hap)])), ]),
            nrow(csp_trap_finalhaplotypes_mosq[csp_trap_finalhaplotypes_mosq$Timepoint %in% unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$trap_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$trap_hap)])), ]),
            nrow(csp_trap_finalhaplotypes_mosq[csp_trap_finalhaplotypes_mosq$Day == 14 & csp_trap_finalhaplotypes_mosq$Timepoint %in% intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$trap_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$trap_hap)]), ]),
            length(unique(variants_Exp52$uniqueid)),
            length(unique(variants_Exp51$uniqueid)),
            length(intersect(variants_Exp51$uniqueid, variants_Exp52$uniqueid)),
            length(unique(variants_Exp51$mosqid)),
            length(unique(variants_Exp51$mosqid[variants_Exp51$uniqueid%in%variants_Exp52$uniqueid]))
  )
)






