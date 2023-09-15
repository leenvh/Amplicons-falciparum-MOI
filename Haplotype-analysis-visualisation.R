library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

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

# Function to process data
process_data <- function(data, type) {
  # Common preprocessing
  data <- data %>%
    filter(grepl(type, Haplotype)) %>%
    mutate(SampleID = gsub("(.*)_.*", "\\1", SampleID),
           Reads = as.numeric(Reads)) %>%
    group_by(SampleID, SampleName) %>%
    summarise(Haplotype = paste(Haplotype, collapse = " "),
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
    summarise(Haplotype = paste(Matching_Values, collapse = " "),
              Reads = paste(Average_Reads, collapse = " ")) %>%
    mutate(MOI = str_count(Haplotype, " ") + 1)
}

# Function to calculate total reads and percentage
calculate_reads <- function(data) {
  data %>%
    mutate(TotalReads = sapply(strsplit(Reads, " "), function(x) sum(as.numeric(x))),
           Percentage = sapply(Reads, function(x) paste(sprintf("%.2f", as.numeric(strsplit(x, " ")[[1]]) / TotalReads * 100), collapse = " ")))
}

# Function to split and process human haplotypes for further data analysis
process_human <- function(data) {
  data_split <- data %>%
    mutate(Individual = str_extract(SampleID, "^[^_]+"),
           Day = str_extract(SampleID, "_([^_]+)$"),
           Day = str_replace_all(Day, "_", ""), 
           Day = as.numeric(str_sub(Day, 4)),
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)),
           Timepoint = SampleID) %>%
    select(Individual, Day, Timepoint, everything(), -SampleID)
  
  return(data_split)
}

# Function to split and process mosquito haplotypes for further data analysis
process_mosquito <- function(data) {
  data_split <- data %>%
    mutate(Individual = str_extract(SampleID, "^[^_]+"),
           Day = str_extract(SampleID, "_([^_]*)_"),  # Capture value between underscores
           Mosquito = str_extract(SampleID, "_([^_]+)$"),
           Day = str_replace_all(Day, "_", ""), 
           Timepoint = paste0(Individual,"_", Day),
           Mosquito = str_replace_all(Mosquito, "_", ""),
           Day = as.numeric(str_sub(Day, 4)),
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
           
    select(Individual, Day, Mosquito, Timepoint, everything(), -SampleID)
  
  return(data_split)
}

# Load and process the data
TRAP_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/TRAP_finalTab.csv",header=TRUE,sep=',')
CSP_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/CSP_finalTab.csv",header=TRUE,sep=',')
Pfs47_1_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/pfs47_1_finalTab.csv",header=TRUE,sep=',')
Pfs47_2_raw <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/pfs47_2_finalTab.csv",header=TRUE,sep=',')

TRAP_processed <- process_data(TRAP_raw, "trap")
CSP_processed <- process_data(CSP_raw, "csp")
Pfs47_1_processed <- process_data(Pfs47_1_raw, "pfs47_1")
Pfs47_2_processed <- process_data(Pfs47_2_raw, "pfs47_2")

# Calculate numbers of samples
length(unique(TRAP_processed$SampleID))
length(unique(CSP_processed$SampleID))
length(unique(Pfs47_1_processed$SampleID))
length(unique(Pfs47_2_processed$SampleID))

# Process Matching Haplotypes
Matching_Haplotypes_CSP <- process_haplotypes(CSP_processed)
Matching_Haplotypes_TRAP <- process_haplotypes(TRAP_processed)
Matching_Haplotypes_Pfs47_1 <- process_haplotypes(Pfs47_1_processed)
Matching_Haplotypes_Pfs47_2 <- process_haplotypes(Pfs47_2_processed)

# Calculate total reads and percentages
Matching_Haplotypes_CSP <- calculate_reads(Matching_Haplotypes_CSP)
Matching_Haplotypes_TRAP <- calculate_reads(Matching_Haplotypes_TRAP)
Matching_Haplotypes_Pfs47_1 <- calculate_reads(Matching_Haplotypes_Pfs47_1)
Matching_Haplotypes_Pfs47_2 <- calculate_reads(Matching_Haplotypes_Pfs47_2)

# Combine and format final haplotype data
FinalHaplotypes <- merge(Matching_Haplotypes_TRAP, 
                         Matching_Haplotypes_CSP, 
                         by = "SampleID", 
                         all = TRUE)
colnames(FinalHaplotypes) <- c("SampleID", "TRAP_hap", "TRAP_reads", "TRAP_MOI", "TRAP_TotalReads", "TRAP_Perc", "CSP_hap", "CSP_reads", "CSP_MOI", "CSP_TotalReads", "CSP_Perc")

Pfs47_FinalHaplotypes <- merge(Matching_Haplotypes_Pfs47_1, 
                         Matching_Haplotypes_Pfs47_2, 
                         by = "SampleID", 
                         all = TRUE)
colnames(Pfs47_FinalHaplotypes) <- c("SampleID", "Pfs47_1_hap", "Pfs47_1_reads", "Pfs47_1_MOI", "Pfs47_1_TotalReads", "Pfs47_1_Perc", "Pfs47_2_hap", "Pfs47_2_reads", "Pfs47_2_MOI", "Pfs47_2_TotalReads", "Pfs47_2_Perc")

# Split by host type and process SampleID
FinalHaplotypes_Human <- process_human(filter(FinalHaplotypes, !grepl("Mosq", SampleID)))
FinalHaplotypes_Mosq <- process_mosquito(filter(FinalHaplotypes, grepl("Mosq", SampleID)))
FINAL_CSP_TRAP<-bind_rows(FinalHaplotypes_Mosq,FinalHaplotypes_Human)
Pfs47_FinalHaplotypes_Human <- process_human(filter(Pfs47_FinalHaplotypes, !grepl("Mosq", SampleID)))
Pfs47_FinalHaplotypes_Mosq <- process_mosquito(filter(Pfs47_FinalHaplotypes, grepl("Mosq", SampleID)))


# Calculate numbers
# individuals
length(unique(FinalHaplotypes_Human$Individual[!is.na(FinalHaplotypes_Human$CSP_hap)])) #49
length(unique(FinalHaplotypes_Human$Individual[!is.na(FinalHaplotypes_Human$TRAP_hap)])) #49
length(unique(Pfs47_FinalHaplotypes_Human$Individual[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_1_hap)])) #33
length(unique(Pfs47_FinalHaplotypes_Human$Individual[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_2_hap)])) #38
# individuals/timepoints
length(na.omit(FinalHaplotypes_Human$CSP_MOI)) #165
length(na.omit(FinalHaplotypes_Human$TRAP_MOI)) #164
length(na.omit(Pfs47_FinalHaplotypes_Human$Pfs47_1_MOI)) #33
length(na.omit(Pfs47_FinalHaplotypes_Human$Pfs47_2_MOI)) #38
# mosquitoes
length(na.omit(FinalHaplotypes_Mosq$CSP_MOI)) #142
length(na.omit(FinalHaplotypes_Mosq$TRAP_MOI)) #129
length(na.omit(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_MOI)) #23
length(na.omit(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_MOI)) #32
# individuals for which there are matching mosquitoes (any timepoint)
length(unique(FinalHaplotypes_Mosq$Individual[!is.na(FinalHaplotypes_Mosq$CSP_hap)])) #36
length(unique(FinalHaplotypes_Mosq$Individual[!is.na(FinalHaplotypes_Mosq$TRAP_hap)])) #35
length(unique(Pfs47_FinalHaplotypes_Mosq$Individual[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_hap)])) #15
length(unique(Pfs47_FinalHaplotypes_Mosq$Individual[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_hap)])) #21
#	Matching timepoints individuals-mosquitoes
length(unique(intersect(FinalHaplotypes_Human$Timepoint[!is.na(FinalHaplotypes_Human$CSP_hap)], FinalHaplotypes_Mosq$Timepoint[!is.na(FinalHaplotypes_Mosq$CSP_hap)])))
length(unique(intersect(FinalHaplotypes_Human$Timepoint[!is.na(FinalHaplotypes_Human$TRAP_hap)], FinalHaplotypes_Mosq$Timepoint[!is.na(FinalHaplotypes_Mosq$TRAP_hap)])))
length(unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_1_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_hap)])))
length(unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_2_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_hap)])))
#	vislualisation of matching timepoints individuals-mosquitoes
matching_timepoints_CSP<-unique(intersect(FinalHaplotypes_Human$Timepoint[!is.na(FinalHaplotypes_Human$CSP_hap)], FinalHaplotypes_Mosq$Timepoint[!is.na(FinalHaplotypes_Mosq$CSP_hap)]))
Matching_Mosq_Haplotypes_CSP<- FinalHaplotypes_Mosq[FinalHaplotypes_Mosq$Timepoint %in% matching_timepoints_CSP,]

ggplot(Matching_Mosq_Haplotypes_CSP, aes(x=factor(Day, levels=c(0, 2, 7, 14, 21, 28, 35)), fill=Day)) +
  geom_bar() +
  labs(x="Day", y="Number of Occurrences - matching human/mosq") +
  theme_classic()









# Melt data for further analysis
melt_data <- function(data, haplotype_type) {
  data %>%
    separate_rows(Haplotype, Reads, Percentage, sep = " ") %>%
    mutate(Amplicon = ifelse(grepl("trap", Haplotype), "TRAP", "CSP"),
           Amplicon = factor(Amplicon, levels = c("TRAP", "CSP")),
           Host = haplotype_type)
}

# Melt and combine data
FinalHaplotypes_Human_Molten <- melt_data(FinalHaplotypes_Human, "Human")
FinalHaplotypes_Mosq_Molten <- melt_data(FinalHaplotypes_Mosq, "Mosquito")

FinalHaplotypes_ALL <- rbind(FinalHaplotypes_Human_Molten, FinalHaplotypes_Mosq_Molten)
FinalHaplotypes_ALL$Host <- factor(FinalHaplotypes_ALL$Host, levels = c("Mosquito", "Human"))

