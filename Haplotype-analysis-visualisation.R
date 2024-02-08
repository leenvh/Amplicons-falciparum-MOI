---
title: "_Plasmodium falciparum_ transmissibility in asymptomatic gametocyte carriers in Mali "
author: "Leen Vanheer"

---
  
#Load libraries
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

# set seed for reproducibility
set.seed(0)

################################################### Functions ###################################################
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
    mutate(Individual = str_extract(SampleID, "^[^_]+"),
           Day = str_extract(SampleID, "_([^_]+)$"),
           Day = str_replace_all(Day, "_", ""), 
           Day = as.numeric(str_sub(Day, 4)),
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)),
           Timepoint = SampleID,
           MOI_Combined = pmax(trap_MOI,csp_MOI, na.rm = TRUE)) %>%
    select(Individual, Day, Timepoint, everything())
  
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
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35)),
           MOI_Combined = pmax(trap_MOI,csp_MOI, na.rm = TRUE)) %>%
    
    select(Individual, Day, Mosquito, Timepoint, everything())
  
  return(data_split)
}

# Function to melt data for further analysis
melt_data <- function(data) {
  data %>%
    separate_rows(Haplotype, Reads, Percentage, sep = " ")%>%
    mutate(Timepoint = gsub("_Mosq\\d+", "", SampleID),
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
    filter("human" %in% species & "mosquito" %in% species) %>%
    ungroup() %>%
    group_by(Timepoint, Haplotype) %>%
    mutate(Comparison = determine_sample_type(cur_data()),
           Day = str_extract(Timepoint, "_([^_]+)$"),  
           Day = str_replace_all(Day, "_", ""), 
           Day = as.numeric(str_sub(Day, 4)),
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
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
    summarise(count = n_distinct(SampleID), .groups = 'drop')
  
  # Step 3: Calculate the percentage for each species
  total_samples_human <- n_distinct(data$SampleID[data$Species == "human"])
  total_samples_mosquito <- n_distinct(data$SampleID[data$Species == "mosquito"])
  
  haplotype_counts %>%
    mutate(frac_sample = case_when(
      Species == "human" ~ (count / total_samples_human) * 100,
      Species == "mosquito" ~ (count / total_samples_mosquito) * 100
    )) %>%
    select(haplotype = !!sym(haplotype_column), Species, frac_sample)
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



################################################### Data visualisation ############################################################
# QC comparing replicates of the same sample
trap_raw %>%
  filter(grepl("^trap-", Haplotype)) %>%
  group_by(SampleID) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  separate(SampleID, into = c("IndividualID", "Rep"), sep = "_Rep") %>%
  mutate(Rep = paste0("Rep", Rep)) %>%
  pivot_wider(names_from = Rep, values_from = Count, values_fill = list(Count = 0)) %>%
  mutate(Marker = "trap") %>%
  bind_rows(
    csp_raw %>%
      filter(grepl("^csp-", Haplotype)) %>%
      group_by(SampleID) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      separate(SampleID, into = c("IndividualID", "Rep"), sep = "_Rep") %>%
      mutate(Rep = paste0("Rep", Rep)) %>%
      pivot_wider(names_from = Rep, values_from = Count, values_fill = list(Count = 0)) %>%
      mutate(Marker = "csp")
  ) %>%
  { 
    # Calculate Pearson correlation coefficients within the pipe using curly braces
    cor_trap <- cor.test(.$Rep1[.$Marker == "trap"], .$Rep2[.$Marker == "trap"], method = "pearson")$estimate
    cor_csp <- cor.test(.$Rep1[.$Marker == "csp"], .$Rep2[.$Marker == "csp"], method = "pearson")$estimate
    
    ggplot(., aes(x = Rep1, y = Rep2, color = Marker)) + 
      geom_jitter(alpha = 0.5) +
      geom_smooth(method = lm) +
      annotate("text", x = 8, y = 2, label = sprintf("trap: r = %.2f", cor_trap), hjust = 0, color = "#4C67BD") +
      annotate("text", x = 8, y = 1, label = sprintf("csp: r = %.2f", cor_csp), hjust = 0, color = "#D90368") +
      theme_classic() +
      xlab("Rep1 MOI") +
      ylab("Rep2 MOI") +
      scale_x_discrete(limits = 0:12) +
      scale_y_discrete(limits = 0:10) +
      ggtitle("Correlation between replicates MOI") +
      scale_color_manual(values = c("trap" = "#4C67BD", "csp" = "#D90368"))
  }
ggsave("Correlation_Replicates.pdf", width = 6, height = 5)




# Percentage monoclonal versus polyclonal at each timepoint - both markers combined
pal<-c("lightgrey", "grey48")

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
      scale_fill_manual(values = pal)
  }
ggsave("MonoPoly_MarkersCombined.pdf", width=6, height=6)



# Percentage matching clones/human only/mosquito only at each timepoint
pal<-c("#a891cf", "#f8997c","lightgrey")
a<- ggplot(csp_matching_haplotypes_molten[csp_matching_haplotypes_molten$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  xlab("Days after treatment initiation")+
  ylab("Percentage of reads")+
  ggtitle("csp")+
  theme(legend.position="none",
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values=pal)
b<- ggplot(trap_matching_haplotypes_molten[trap_matching_haplotypes_molten$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  #facet_wrap(~Matching_FinalHaplotypes_ALL$Host)
  xlab("Days after treatment initiation")+
  ylab("")+
  scale_fill_manual(values=pal, name =" ")+
  ggtitle("trap")+
  theme(axis.text.x = element_text(size = 10), 
      axis.text.y = element_text(size = 10))
combined_plot <- a + b + plot_layout(ncol=2, widths=c(1.1, 1))
ggsave("Fig1B_Matching_Haplotypes_Percentages.pdf", combined_plot, width=14, height=6)



# Percentage matching clones/human only/mosquito comparing day 0 mosquito clones to all human timepoints
pal<-c("#a891cf","lightgrey")

individuals_with_day0_mosquito<-csp_haplotypes_molten %>%
  filter(species == "mosquito", Day == 0) %>%
  select(Individual) %>%
  distinct()

# Step 2: Filter human samples for those individuals
human_samples_for_matching_individuals <- csp_haplotypes_molten %>%
  filter(species == "human") %>%
  semi_join(individuals_with_day0_mosquito, by = "Individual")

# Step 3: Identify haplotypes in day 0 mosquito samples for matching
mosquito_day0_haplotypes <- csp_haplotypes_molten %>%
  filter(species == "mosquito", Day == 0) %>%
  select(Individual, Haplotype) %>%
  distinct()

# Annotate human samples
annotated_human_samples <- human_samples_for_matching_individuals %>%
  rowwise() %>%
  mutate(Status = case_when(
    Haplotype %in% mosquito_day0_haplotypes$Haplotype[mosquito_day0_haplotypes$Individual == Individual] ~ "matching",
    TRUE ~ "human_only"
  )) %>%
  ungroup() 

a<-ggplot(annotated_human_samples[annotated_human_samples$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Status))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  xlab("Days after treatment initiation")+
  ylab("Percentage of reads")+
  ggtitle("trap")+
  theme(legend.position="none",
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values=pal)


# 3 instances of mosquito-only clones in day 0 mosquitoes (that were not found at any human timepoint)
mosquito_human_comparison <- mosquito_day0_haplotypes %>%
  left_join(annotated_human_samples %>% select(Individual, Haplotype, Status) %>%
              distinct(Individual, Haplotype) %>%
              mutate(human_present = TRUE), 
            by = c("Individual", "Haplotype")) %>%
  mutate(Status = if_else(is.na(human_present), "mosquito_only", "matching"))

#trap
individuals_with_day0_mosquito<-trap_haplotypes_molten %>%
  filter(species == "mosquito", Day == 0) %>%
  select(Individual) %>%
  distinct()

# Step 2: Filter human samples for those individuals
human_samples_for_matching_individuals <- trap_haplotypes_molten %>%
  filter(species == "human") %>%
  semi_join(individuals_with_day0_mosquito, by = "Individual")

# Step 3: Identify haplotypes in day 0 mosquito samples for matching
mosquito_day0_haplotypes <- trap_haplotypes_molten %>%
  filter(species == "mosquito", Day == 0) %>%
  select(Individual, Haplotype) %>%
  distinct()

# Annotate human samples
annotated_human_samples <- human_samples_for_matching_individuals %>%
  rowwise() %>%
  mutate(Status = case_when(
    Haplotype %in% mosquito_day0_haplotypes$Haplotype[mosquito_day0_haplotypes$Individual == Individual] ~ "matching",
    TRUE ~ "human_only"
  )) %>%
  ungroup() 

b<-ggplot(annotated_human_samples[annotated_human_samples$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Status))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  xlab("Days after treatment initiation")+
  ylab("Percentage of reads")+
  ggtitle("csp")+
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values=pal)


# 11 instances of mosquito-only clones in day 0 mosquitoes (that were not found at any human timepoint)
mosquito_human_comparison <- mosquito_day0_haplotypes %>%
  left_join(annotated_human_samples %>% select(Individual, Haplotype, Status) %>%
              distinct(Individual, Haplotype) %>%
              mutate(human_present = TRUE), 
            by = c("Individual", "Haplotype")) %>%
  mutate(Status = if_else(is.na(human_present), "mosquito_only", "matching"))



combined_plot<- a + b + plot_layout(ncol=2, widths=c(1.1, 1))
ggsave("Mosq0_Matching_Haplotypes_Percentages.pdf", combined_plot, width=14, height=6)



#Check correlation between csp MOI and trap MOI
csp_trap_finalhaplotypes %>%
{
  # Calculate Pearson correlation coefficients within the pipe using curly braces
  cor <- cor.test(.$csp_MOI, .$trap_MOI, method = "pearson")$estimate
  cor_label <- sprintf("Pearson r: %.2f", cor)
  
  ggplot(csp_trap_finalhaplotypes, aes(x=csp_MOI, y=trap_MOI)) + 
  geom_jitter()+
  theme_classic()+
  theme(legend.position = "none")+
  annotate("text", x = 6, y = 5, label = cor_label, hjust = 0) +
  xlab("csp MOI")+
  ylab("trap MOI")+
  scale_x_discrete(limits = 0:12) +
  scale_y_discrete(limits = 0:10) +
  geom_smooth(method=lm)+
  ggtitle("Correlation between csp And trap MOI")
}
ggsave("Correlation_csp_trap.pdf", width=6, height=5)




###### Haplotype prevalence in human vs mosquito samples
csp_trap_frac %>%
  group_by(haplotype) %>%
  summarise(frac_sample_human = sum(frac_sample[Species == "human"], na.rm = TRUE),
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



#Check clonality over time in the individuals and mosquitoes
c<-ggplot(csp_trap_finalhaplotypes[csp_trap_finalhaplotypes$Species == "human" & csp_trap_finalhaplotypes$Day%in% c("0","2","7","14","21","28"),], aes(x=Day, y=MOI_Combined)) + 
  geom_boxplot(aes(fill=Day), alpha = 0.5, fill="lightgrey", outlier.shape = NA) +
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Days after treatment initiation", y ="MOI")+
  scale_y_continuous(breaks = seq(0, 10, by = 1)) + 
  ggtitle("Clonality in individuals") 

d <- ggplot(csp_trap_finalhaplotypes[csp_trap_finalhaplotypes$Species == "mosquito",], aes(x=Day, y=MOI_Combined)) + 
  geom_boxplot(alpha = 0.5, fill="lightgrey", outlier.shape = NA) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x= "Days after treatment initiation", y ="MOI") +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) + 
  ggtitle("Clonality in mosquitoes")

combined_clonality_plot <- c + d + plot_layout(ncol=2)
ggsave("Clonality.pdf", combined_clonality_plot, width=10, height=5)


#Check relation between clonality and gametocytemia
ggplot(merge_clinical_haplotypes, aes(x=as.factor(MOI_Combined), y=totalgct_ul)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA,fill="grey")+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Gametocyte density")+
  scale_y_log10()+
  ggtitle("Clonality ~ Gametocytaemia (all timepoints)")
ggsave("Clonaty_GC.pdf", width=6, height=5)

#For Day0 only
ggplot(merge_clinical_haplotypes[merge_clinical_haplotypes$Day=="0",], aes(x=as.factor(MOI_Combined), y=totalgct_ul)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA,fill="grey")+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Gametocyte density")+
  scale_y_log10()+
  ggtitle("Clonality ~ Gametocytaemia")


#Check relation between individual clonality and infection rates
ggplot(merge_clinical_haplotypes, aes(x=as.factor(MOI_Combined), y=percentagemosqinfected)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA,fill="grey")+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Infected Mosquitoes (%)")+
  ggtitle("Individual clonality ~ infection rate")
ggsave("Infection_rate_clonality.pdf", width=6, height=5)

#Monoclonal versus polyclonal
merge_clinical_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  ggplot(aes(x=as.factor(MonoPoly), y=percentagemosqinfected)) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA, fill="grey")+
    theme_classic()+
    geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
    theme(legend.position = "none")+
    labs(x= "MOI",y="Infected Mosquitoes (%)")+
    ggtitle("Individual clonality ~ infection rate")
ggsave("Infection_rate_MonoPolyclonal.pdf", width=6, height=5)

#Check relation between clonality and infection rates (individuals only) - normalised by GC
merge_clinical_haplotypes %>%
  mutate(normalized_infection_rate = (percentagemosqinfected / totalgct_ul) * 100) %>%
  ggplot(aes(x=as.factor(MOI_Combined), y=normalized_infection_rate)) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA, fill="grey")+
    theme_classic()+
    geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
    theme(legend.position = "none")+
    labs(x= "MOI",y="Infected Mosquitoes (%)")+
    ggtitle("Individual clonality ~ infection rate (normalised by GC density)")
ggsave("Infection_rate_clonality_Normalised_by_GC.pdf", width=6, height=5)

#Monoclonal versus polyclonal
merge_clinical_haplotypes %>%
  mutate(normalized_infection_rate = (percentagemosqinfected / totalgct_ul) * 100,
         MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  ggplot(aes(x=as.factor(MonoPoly), y=normalized_infection_rate)) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA, fill="grey")+
    theme_classic()+
    geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
    theme(legend.position = "none")+
    labs(x= "MOI",y="Infected Mosquitoes (%)")+
    ggtitle("Individual clonality ~ infection rate (normalised by GC density)")
ggsave("Infection_rate_MonoPolyclonal_Normalised_by_GC.pdf", width=6, height=5)



#Investigate not-infectious samples at day 0
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


#Check relation between clonality and age (individuals only) at Day 0
merge_clinical_haplotypes %>%
  mutate(
    age_group = case_when(
      age < 5 ~ "<5",
      age >= 5 & age < 10 ~ "5-10",
      age >= 10 & age < 15 ~ "10-15",
      age >= 15 & age < 20 ~ "15-20",
      age >= 20 ~ ">20"
    ),
    age_group = factor(age_group, levels = c("<5", "5-10", "10-15", "15-20", ">20"))
  ) %>%
  group_by(age_group) %>%
  mutate(mean_MOI_by_age = mean(MOI_Combined, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Day == "0") %>%
  ggplot(aes(x = age, y = mean_MOI_by_age)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Age", y = "Mean MOI") +
  scale_y_continuous(limits = c(1, 3), breaks = c(1, 2, 3)) +
  ggtitle("Mean clonality by age")
ggsave("Clonality_by_agegroup.pdf", width=6, height=5)

merge_clinical_haplotypes %>%
  group_by(age) %>%
  mutate(mean_MOI_by_age = mean(MOI_Combined, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Day == "0") %>%
  ggplot(aes(x = age, y = mean_MOI_by_age)) +
  geom_line() +
  geom_point()+
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Age", y = "Mean MOI") +
  scale_y_continuous(limits = c(1, 4), breaks = c(1, 2, 4)) +
  ggtitle("Mean clonality by age")
ggsave("Clonality_by_age.pdf", width=6, height=5)

#Check relation between midgut clonality and oocyst number (mosquitoes only)
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

#Check relation between individual clonality and oocyst number (mosquitoes only)
ggplot(merge_oocyst_individual_haplotypes, aes(x=as.factor(MOI_Combined), y=oocysts)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI individual",y="Oocyst number")+
  ggtitle("Individual clonality ~ oocyst number")
ggsave("Individual_Clonality_Oocyst_number.pdf", width=6, height=5)

#Monoclonal versus polyclonal
merge_oocyst_individual_haplotypes %>%
  mutate(MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  ggplot(aes(x=as.factor(MonoPoly), y=oocysts)) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    theme_classic()+
    geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
    theme(legend.position = "none")+
    labs(x= "MOI individual",y="Oocyst number")+
    ggtitle("Individual clonality ~ oocyst number")
ggsave("Individual_Clonality_MonoPoly_Oocyst_number.pdf", width=6, height=5)

#Check relation between individual clonality and oocyst number (mosquitoes only) - normalised by GC
merge_oocyst_clinical_individual_haplotypes %>%
  mutate(normalized_oocysts = (oocysts / totalgct_ul) * 100) %>%
  ggplot(aes(x=as.factor(MOI_Combined), y=normalized_oocysts)) + 
    geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    theme_classic()+
    geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
    theme(legend.position = "none")+
    labs(x= "MOI individual",y="Oocyst number")+
    ggtitle("Individual clonality ~ oocyst number (normalised by GC)")
ggsave("Individual_Clonality_Oocyst_number_normalised_by_GC.pdf", width=6, height=5)

#Monoclonal versus polyclonal
merge_oocyst_clinical_individual_haplotypes %>%
  mutate(normalized_oocysts = (oocysts / totalgct_ul) * 100) %>%
  mutate(MonoPoly = ifelse(MOI_Combined=="1","monoclonal","polyclonal"))%>%
  ggplot(aes(x=as.factor(MonoPoly), y=normalized_oocysts)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), alpha=0.5)+
  theme(legend.position = "none")+
  labs(x= "MOI individual",y="Oocyst number")+
  ggtitle("Individual clonality ~ oocyst number (normalised by GC)")
ggsave("Individual_Clonality_MonoPoly_Oocyst_number_normalised_by_GC.pdf", width=6, height=5)


#Check relation between clonality and month (individuals only)
merge_clinical_haplotypes %>%
  mutate(visitdate_hb = as.Date(visitdate_hb, format="%d/%m/%Y")) %>%
  group_by(visitdate_hb) %>%
  mutate(mean_MOI_by_time = mean(MOI_Combined, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x=visitdate_hb, y=mean_MOI_by_time)) + 
    geom_line()+
    theme_classic()+
    theme(legend.position = "none")+
    labs(x= "Visitdate",y="Mean MOI")+
    ggtitle("Mean clonality over time (individuals only)")
ggsave("Clonality_by_month.pdf", width=6, height=5)



#Check accross haplotypes which ones are matching or not
csp_grouped_df <- csp_matching_haplotypes_molten %>%
  group_by(Haplotype, Comparison) %>%
  summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
csp_grouped_df <- csp_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("csp-", "", Haplotype))) %>%
  arrange(Haplotype_num)

trap_grouped_df <- trap_matching_haplotypes_molten %>%
  group_by(Haplotype, Comparison) %>%
  summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
trap_grouped_df <- trap_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("trap-", "", Haplotype))) %>%
  arrange(Haplotype_num)

# Plot
e<-ggplot(csp_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#9AC77B', "grey", '#7EC2DE'))+
  labs(title = "csp",
       x = "Haplotype",
       y = "Count") +
  theme_light()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
f<-ggplot(trap_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#9AC77B', "grey", '#7EC2DE'))+
  labs(title = "trap",
       x = "Haplotype",
       y = "Count") +
  theme_light()+
  theme(axis.text.x = element_blank())
combined_plot <- e / f 
ggsave("Haplotypes_transmission.pdf", combined_plot, width=5, height=7)




#Check whether multiclonal samples have higher likelihood of causing multiclonal mosquito infections

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






# Calculate the odds of transmission for each haplotype
Clin_csp <- csp_haplotypes_molten %>%
  left_join(select(clinicaldata, SampleID, totalgct_ul, ring_ul_all), by = "SampleID")

# Normalize read percentages in human data
human_data_csp <- Clin_csp %>% 
  filter(species != 'mosquito') %>%
  mutate(Normalized_csp_Perc = as.numeric(Percentage) / as.numeric(totalgct_ul))

# Prepare mosquito data
mosquito_data_csp <- Clin_csp %>% 
  filter(species == 'mosquito')

# Merge human and mosquito data
merged_data_csp <- merge(human_data_csp, mosquito_data_csp, by = c("Timepoint", "Haplotype"), all = TRUE)

merged_data_csp$Normalized_csp_Perc <- ifelse(is.na(merged_data_csp$Normalized_csp_Perc), 0, merged_data_csp$Normalized_csp_Perc)
merged_data_csp$Percentage.y <- ifelse(is.na(merged_data_csp$Percentage.y), 0, merged_data_csp$Percentage.y)


# Add a column to indicate if the read percentage is higher in mosquito sample
merged_data_csp$Higher_in_Mosquito <- merged_data_csp$Percentage.y > merged_data_csp$Normalized_csp_Perc

# Calculate the odds of transmission for each haplotype
haplotype_odds <- merged_data_csp %>%
  group_by(Haplotype) %>%
  summarize(Odds_of_Transmission = mean(Higher_in_Mosquito, na.rm = TRUE))

haplotype_odds$Haplotype_Num <- as.numeric(sub("csp-", "", haplotype_odds$Haplotype))

# Ordering the dataframe by this numeric part
haplotype_odds <- haplotype_odds[order(haplotype_odds$Haplotype_Num), ]

# Removing the auxiliary numeric column if not needed further
haplotype_odds$Haplotype_Num <- NULL

y_axis_data <- data.frame(Y = seq(0, 1, by = 0.1), Odds_of_Transmission = 0)

# Create the circular plot with the y-axis circles
a<-ggplot(haplotype_odds, aes(x = factor(Haplotype), y = Odds_of_Transmission, fill = Odds_of_Transmission)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = ifelse(Odds_of_Transmission > 0.1, Haplotype, ""), y = Odds_of_Transmission + 0.02), 
            position = position_stack(vjust = 1), 
            color = "black", 
            size = 3) +
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



#########trap

Clin_trap <- trap_haplotypes_molten %>%
  left_join(select(clinicaldata, SampleID, totalgct_ul, ring_ul_all), by = "SampleID")

#write.csv(Clin_csp_trap,"Clin_csp_trap.csv", row.names = FALSE)

# Normalize read percentages in human data
human_data_trap <- Clin_trap %>% 
  filter(species != 'mosquito') %>%
  mutate(Normalized_trap_Perc = as.numeric(Percentage) / as.numeric(totalgct_ul))

# Prepare mosquito data
mosquito_data_trap <- Clin_trap %>% 
  filter(species == 'mosquito')

# Merge human and mosquito data
merged_data_trap <- merge(human_data_trap, mosquito_data_trap, by = c("Timepoint", "Haplotype"), all = TRUE)

merged_data_trap$Normalized_trap_Perc <- ifelse(is.na(merged_data_trap$Normalized_trap_Perc), 0, merged_data_trap$Normalized_trap_Perc)
merged_data_trap$Percentage.y <- ifelse(is.na(merged_data_trap$Percentage.y), 0, merged_data_trap$Percentage.y)


# Add a column to indicate if the read percentage is higher in mosquito sample
merged_data_trap$Higher_in_Mosquito <- merged_data_trap$Percentage.y > merged_data_trap$Normalized_trap_Perc

# Calculate the odds of transmission for each haplotype
haplotype_odds <- merged_data_trap %>%
  group_by(Haplotype) %>%
  summarize(Odds_of_Transmission = mean(Higher_in_Mosquito, na.rm = TRUE))

haplotype_odds$Haplotype_Num <- as.numeric(sub("trap-", "", haplotype_odds$Haplotype))

# Ordering the dataframe by this numeric part
haplotype_odds <- haplotype_odds[order(haplotype_odds$Haplotype_Num), ]

# Removing the auxiliary numeric column if not needed further
haplotype_odds$Haplotype_Num <- NULL

b<-ggplot(haplotype_odds, aes(x = factor(Haplotype), y = Odds_of_Transmission, fill = Odds_of_Transmission)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = Haplotype), 
            position = position_stack(vjust = 1), 
            color = "black", 
            size = 3) +
  scale_fill_gradient2(low = "white",high = "lightblue4",limits=c(0,1)) +
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
    title=element_blank()
  )

combinedplot<-a+b
ggsave("OddsOfTransmisson.pdf", width=15, height=7)


# Overview of all haplotypes
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




#Calculate numbers of samples
table_numbers <- data.frame(
  Description = c("Total samples with csp_hap",
                  "Total samples with trap_hap",
                  "Unique individuals with csp_hap",
                  "Unique individuals with trap_hap",
                  "Timepoints with csp_hap",
                  "Timpeoints with trap_hap",
                  "Individuals non-infectious at day 0",
                  "Mosquito samples with csp_hap",
                  "Mosquito samples with trap_hap",
                  "Individuals for which there are matching mosquitoes with csp_hap (any timepoint)",
                  "Individuals for which there are matching mosquitoes with trap_hap (any timepoint)",
                  "Timepoints at which there are matching human-mosquito samples with csp_hap",
                  "Timepoints at which there are matching human-mosquito samples with trap_hap",
                  "Matching samples human-mosquito with csp_hap",
                  "Matching samples human-mosquito with trap_hap"),
  Count = c(length(unique(csp_processed$SampleID)), length(unique(trap_processed$SampleID)),
            length(unique(csp_trap_finalhaplotypes_human$Individual[!is.na(csp_trap_finalhaplotypes_human$csp_hap)])), length(unique(csp_trap_finalhaplotypes_human$Individual[!is.na(csp_trap_finalhaplotypes_human$trap_hap)])),
            length(na.omit(csp_trap_finalhaplotypes_human$csp_MOI)), length(na.omit(csp_trap_finalhaplotypes_human$trap_MOI)),
            sum(merge_clinical_haplotypes$percentagemosqinfected == 0 & merge_clinical_haplotypes$Day == 0, na.rm = TRUE),
            length(na.omit(csp_trap_finalhaplotypes_mosq$csp_MOI)),length(na.omit(csp_trap_finalhaplotypes_mosq$trap_MOI)),
            length(unique(csp_trap_finalhaplotypes_mosq$Individual[!is.na(csp_trap_finalhaplotypes_mosq$csp_hap)])),length(unique(csp_trap_finalhaplotypes_mosq$Individual[!is.na(csp_trap_finalhaplotypes_mosq$trap_hap)])),
            length(unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$csp_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$csp_hap)]))),
            length(unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$trap_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$trap_hap)]))),
            nrow(csp_trap_finalhaplotypes_mosq[csp_trap_finalhaplotypes_mosq$Timepoint %in% unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$csp_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$csp_hap)])), ]),
            nrow(csp_trap_finalhaplotypes_mosq[csp_trap_finalhaplotypes_mosq$Timepoint %in% unique(intersect(csp_trap_finalhaplotypes_human$Timepoint[!is.na(csp_trap_finalhaplotypes_human$trap_hap)], csp_trap_finalhaplotypes_mosq$Timepoint[!is.na(csp_trap_finalhaplotypes_mosq$trap_hap)])), ])
  )
)




  