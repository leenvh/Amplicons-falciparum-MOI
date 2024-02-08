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

#length(unique(Pfs47_1_processed$SampleID))
#length(unique(Pfs47_2_processed$SampleID))
#length(unique(Pfs47_FinalHaplotypes_Human$Individual[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_1_hap)])) #33
#length(unique(Pfs47_FinalHaplotypes_Human$Individual[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_2_hap)])) #38
#length(na.omit(Pfs47_FinalHaplotypes_Human$Pfs47_1_MOI)) #33
#length(na.omit(Pfs47_FinalHaplotypes_Human$Pfs47_2_MOI)) #38
#length(na.omit(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_MOI)) #23
#length(na.omit(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_MOI)) #32
#length(unique(Pfs47_FinalHaplotypes_Mosq$Individual[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_hap)])) #15
#length(unique(Pfs47_FinalHaplotypes_Mosq$Individual[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_hap)])) #21
#length(unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_1_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_hap)])))
#length(unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_2_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_hap)])))

#matching_timepoints_Pfs47_1<-unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_1_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_hap)]))
#Matching_Mosq_Haplotypes_Pfs47_1<- Pfs47_FinalHaplotypes_Mosq[Pfs47_FinalHaplotypes_Mosq$Timepoint %in% matching_timepoints_Pfs47_1,]
#nrow(Matching_Mosq_Haplotypes_Pfs47_1) #14
#matching_timepoints_Pfs47_2<-unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_2_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_hap)]))
#Matching_Mosq_Haplotypes_Pfs47_2<- Pfs47_FinalHaplotypes_Mosq[Pfs47_FinalHaplotypes_Mosq$Timepoint %in% matching_timepoints_Pfs47_2,]
#nrow(Matching_Mosq_Haplotypes_Pfs47_2) #29



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





#Pfs47
x<-ggplot(FINAL_Pfs47[FINAL_Pfs47$Species == "human",], aes(x=Day, y=Pfs47_1_MOI,fill=Day)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.1), aes(colour=Day), alpha=0.9)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Visit Day", y ="Pfs47_1 haplotypes")+
  scale_y_continuous(limits=c(1,3), breaks = c(0, 1, 2, 3))+
  ggtitle("Individuals")

y<-ggplot(FINAL_Pfs47[FINAL_Pfs47$Species == "mosquito",], aes(x=Day, y=Pfs47_1_MOI,fill=Day)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.1), aes(colour=Day), alpha=0.9)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Visit Day", y ="Pfs47_1 haplotypes")+
  scale_y_continuous(limits=c(1,3), breaks = c(0, 1, 2, 3))+
  ggtitle("Mosquitoes")

z<-ggplot(FINAL_Pfs47[FINAL_Pfs47$Species == "human",], aes(x=Day, y=Pfs47_2_MOI,fill=Day)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.1), aes(colour=Day), alpha=0.9)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Visit Day", y ="Pfs47_2 haplotypes")+
  scale_y_continuous(limits=c(1,3), breaks = c(0, 1, 2, 3))+
  ggtitle("Individuals")

w<-ggplot(FINAL_Pfs47[FINAL_Pfs47$Species == "mosquito",], aes(x=Day, y=Pfs47_2_MOI,fill=Day)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.1), aes(colour=Day), alpha=0.9)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Visit Day", y ="Pfs47_2 haplotypes")+
  scale_y_continuous(limits=c(1,3), breaks = c(0, 1, 2, 3))+
  ggtitle("Mosquitoes")

combined_clonality_plot <- x + y + z + w + plot_layout(ncol=4)
combined_clonality_plot
ggsave("Pfs47_Clonality.pdf", combined_clonality_plot, width=10, height=5)



#Check accross Pfs47 haplotypes which ones are matching or not
Pfs47_1_grouped_df <- Matching_Haplotypes_Pfs47_1_Molten_matching %>%
  group_by(Haplotype, Comparison) %>%
  summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
Pfs47_1_grouped_df <- Pfs47_1_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("pfs47_1-", "", Haplotype))) %>%
  arrange(Haplotype_num)

Pfs47_2_grouped_df <- Matching_Haplotypes_Pfs47_2_Molten_matching %>%
  group_by(Haplotype, Comparison) %>%
  summarise(Count = n()) %>%
  ungroup()

# Extract the number from "pfs-" to order the haplotypes
Pfs47_2_grouped_df <- Pfs47_2_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("pfs47_2-", "", Haplotype))) %>%
  arrange(Haplotype_num)

# Plot
g<-ggplot(Pfs47_1_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Pfs47_1",
       x = "Haplotype",
       y = "Count") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
h<-ggplot(Pfs47_2_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Pfs47_2",
       x = "Haplotype",
       y = "Count") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
combined_plot <- g + h + plot_annotation(title ="Are any Pfs47 haplotypes more likely to transmit than others?") 
ggsave("Pfs47_haplotypes_transmission.pdf", combined_plot, width=11, height=5)


