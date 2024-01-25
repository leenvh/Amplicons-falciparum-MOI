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

#QC to check that replicates are similar
trap_counts <- TRAP_raw %>%
  filter(grepl("^trap-", Haplotype)) %>%
  group_by(SampleID) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  separate(SampleID, into = c("IndividualID", "Rep"), sep = "_Rep") %>%
  mutate(Rep = paste0("Rep", Rep)) %>%
  pivot_wider(names_from = Rep, values_from = Count, values_fill = list(Count = 0)) %>%
  mutate(Marker = "TRAP")

# Process CSP data (similar to TRAP, adjust the grepl pattern and mutate(Source = "CSP"))
csp_counts <- CSP_raw %>%
  filter(grepl("^csp-", Haplotype)) %>%
  group_by(SampleID) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  separate(SampleID, into = c("IndividualID", "Rep"), sep = "_Rep") %>%
  mutate(Rep = paste0("Rep", Rep)) %>%
  pivot_wider(names_from = Rep, values_from = Count, values_fill = list(Count = 0)) %>%
  mutate(Marker = "CSP")

# Combine TRAP and CSP data
combined_counts <- bind_rows(trap_counts, csp_counts)

# Plotting the combined data
ggplot(combined_counts, aes(x=Rep1, y=Rep2, color=Marker)) + 
  geom_jitter(alpha=0.5) +
  theme_classic() +
  xlab("Rep1 MOI") +
  ylab("Rep2 MOI") +
  xlim(0, 10) +
  ylim(0, 10)+
  geom_smooth(method=lm) +
  ggtitle("Correlation between replicates MOI") +
  scale_color_manual(values = c("TRAP" = "#4C67BD", "CSP" = "#D90368"))
ggsave("Correlation_Replicates.pdf", width=6, height=5)

cor.test(combined_counts$Rep1[combined_counts$Source=="TRAP"], combined_counts$Rep2[combined_counts$Source=="TRAP"], method = c("pearson")) #0.8672038 
cor.test(combined_counts$Rep1[combined_counts$Source=="CSP"], combined_counts$Rep2[combined_counts$Source=="CSP"], method = c("pearson")) #0.8422626 


#heatmap
ggplot(combined_counts, aes(x = Rep1, y = Rep2)) +
  geom_tile(stat = "bin_2d", aes(fill = ..count..), bins = 10) +  # Adjust 'bins' as needed
  scale_fill_gradientn(colours = c("orange", "#D90368", "blue","green"),  # Add more colors as needed
                       values = scales::rescale(c(0, 0.3, 0.6, 1)),  # Adjust these values as needed
                       limits = c(0, 80)) +
  theme_classic() +
  xlab("Rep1 MOI") +
  ylab("Rep2 MOI") +
  facet_wrap(~Marker)+
  geom_smooth(method=lm) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10))+
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10))+
  ggtitle("Heatmap of Correlation between Replicates MOI")


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
           Timepoint = SampleID) %>%
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
           Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
           
    select(Individual, Day, Mosquito, Timepoint, everything())
  
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
Haplotypes_CSP <- process_haplotypes(CSP_processed)
Haplotypes_TRAP <- process_haplotypes(TRAP_processed)
Haplotypes_Pfs47_1 <- process_haplotypes(Pfs47_1_processed)
Haplotypes_Pfs47_2 <- process_haplotypes(Pfs47_2_processed)

# Calculate total reads and percentages
Haplotypes_CSP <- calculate_reads(Haplotypes_CSP)
Haplotypes_TRAP <- calculate_reads(Haplotypes_TRAP)
Haplotypes_Pfs47_1 <- calculate_reads(Haplotypes_Pfs47_1)
Haplotypes_Pfs47_2 <- calculate_reads(Haplotypes_Pfs47_2)

# Assign species
Haplotypes_CSP$species <- ifelse(grepl("Mosq", Haplotypes_CSP$SampleID), "mosquito", "human")
Haplotypes_TRAP$species <- ifelse(grepl("Mosq", Haplotypes_TRAP$SampleID), "mosquito", "human")
Haplotypes_Pfs47_1$species <- ifelse(grepl("Mosq", Haplotypes_Pfs47_1$SampleID), "mosquito", "human")
Haplotypes_Pfs47_2$species <- ifelse(grepl("Mosq", Haplotypes_Pfs47_2$SampleID), "mosquito", "human")

# Combine and format final haplotype data
FinalHaplotypes <- merge(Haplotypes_TRAP, 
                         Haplotypes_CSP, 
                         by = "SampleID", 
                         all = TRUE)
FinalHaplotypes$species.x <- NULL
FinalHaplotypes$species.y <- NULL
colnames(FinalHaplotypes) <- c("SampleID", "TRAP_hap", "TRAP_reads", "TRAP_MOI", "TRAP_TotalReads", "TRAP_Perc", "CSP_hap", "CSP_reads", "CSP_MOI", "CSP_TotalReads", "CSP_Perc")

Pfs47_FinalHaplotypes <- merge(Haplotypes_Pfs47_1, 
                         Haplotypes_Pfs47_2, 
                         by = "SampleID", 
                         all = TRUE)
Pfs47_FinalHaplotypes$species.x <- NULL
Pfs47_FinalHaplotypes$species.y <- NULL
colnames(Pfs47_FinalHaplotypes) <- c("SampleID", "Pfs47_1_hap", "Pfs47_1_reads", "Pfs47_1_MOI", "Pfs47_1_TotalReads", "Pfs47_1_Perc", "Pfs47_2_hap", "Pfs47_2_reads", "Pfs47_2_MOI", "Pfs47_2_TotalReads", "Pfs47_2_Perc")



# Split by host type and process SampleID
FinalHaplotypes_Human <- process_human(filter(FinalHaplotypes, !grepl("Mosq", SampleID)))
FinalHaplotypes_Mosq <- process_mosquito(filter(FinalHaplotypes, grepl("Mosq", SampleID)))
FINAL_CSP_TRAP<-bind_rows(FinalHaplotypes_Mosq,FinalHaplotypes_Human)
FINAL_CSP_TRAP$Species <- ifelse(is.na(FINAL_CSP_TRAP$Mosquito), "human", "mosquito")
write.csv(FINAL_CSP_TRAP, "FINAL_CSP_TRAP.csv", row.names = FALSE)
Pfs47_FinalHaplotypes_Human <- process_human(filter(Pfs47_FinalHaplotypes, !grepl("Mosq", SampleID)))
Pfs47_FinalHaplotypes_Mosq <- process_mosquito(filter(Pfs47_FinalHaplotypes, grepl("Mosq", SampleID)))
FINAL_Pfs47<-bind_rows(Pfs47_FinalHaplotypes_Mosq,Pfs47_FinalHaplotypes_Human)
FINAL_Pfs47$Species <- ifelse(is.na(FINAL_Pfs47$Mosquito), "human", "mosquito")

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
# number of those human timepoints that are infectious/non-infectious
ClinicalData<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Paper_AmpSeq/ClinicalData.csv",header=TRUE,sep=',')
ClinicalData$SampleID<-paste(ClinicalData$studycode,"_Day",ClinicalData$studyvisit,sep="")
FinalHaplotypes_Human_clinical <- FinalHaplotypes_Human %>%
  left_join(ClinicalData, by = "SampleID")
FinalHaplotypes_Human_clinical$Infectious<- "Yes"
FinalHaplotypes_Human_clinical$Infectious[FinalHaplotypes_Human_clinical$mosq_pos == 0]<- "No"
sum(!is.na(FinalHaplotypes_Human_clinical$TRAP_hap[FinalHaplotypes_Human_clinical$Infectious == "Yes"])) #109
sum(!is.na(FinalHaplotypes_Human_clinical$TRAP_hap[FinalHaplotypes_Human_clinical$Infectious == "No"])) #55
sum(!is.na(FinalHaplotypes_Human_clinical$CSP_hap[FinalHaplotypes_Human_clinical$Infectious == "Yes"])) #109
sum(!is.na(FinalHaplotypes_Human_clinical$CSP_hap[FinalHaplotypes_Human_clinical$Infectious == "No"])) #56
sum(FinalHaplotypes_Human_clinical$Infectious == "No" & FinalHaplotypes_Human_clinical$Day == 0, na.rm = TRUE) #14 Number of non-infectious individuals at Day0
sum(!is.na(FinalHaplotypes_Human_clinical$CSP_hap[FinalHaplotypes_Human_clinical$Infectious == "No" & FinalHaplotypes_Human_clinical$Day == 0])) #14
sum(!is.na(FinalHaplotypes_Human_clinical$TRAP_hap[FinalHaplotypes_Human_clinical$Infectious == "No" & FinalHaplotypes_Human_clinical$Day == 0])) #14
FinalHaplotypes_Human_clinical$Individual[FinalHaplotypes_Human_clinical$Infectious == "No" & FinalHaplotypes_Human_clinical$Day == 0] #14
ClinicalData$studycode[ClinicalData$mosq_pos == 0 & ClinicalData$studyvisit == 0 & ClinicalData$arm_num%in% c(1,3)] #17
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
#	Matching samples individuals-mosquitoes
matching_timepoints_CSP<-unique(intersect(FinalHaplotypes_Human$Timepoint[!is.na(FinalHaplotypes_Human$CSP_hap)], FinalHaplotypes_Mosq$Timepoint[!is.na(FinalHaplotypes_Mosq$CSP_hap)]))
Matching_Mosq_Haplotypes_CSP<- FinalHaplotypes_Mosq[FinalHaplotypes_Mosq$Timepoint %in% matching_timepoints_CSP,]
nrow(Matching_Mosq_Haplotypes_CSP) #130
matching_timepoints_TRAP<-unique(intersect(FinalHaplotypes_Human$Timepoint[!is.na(FinalHaplotypes_Human$TRAP_hap)], FinalHaplotypes_Mosq$Timepoint[!is.na(FinalHaplotypes_Mosq$TRAP_hap)]))
Matching_Mosq_Haplotypes_TRAP<- FinalHaplotypes_Mosq[FinalHaplotypes_Mosq$Timepoint %in% matching_timepoints_TRAP,]
nrow(Matching_Mosq_Haplotypes_TRAP) #121
matching_timepoints_Pfs47_1<-unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_1_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_1_hap)]))
Matching_Mosq_Haplotypes_Pfs47_1<- Pfs47_FinalHaplotypes_Mosq[Pfs47_FinalHaplotypes_Mosq$Timepoint %in% matching_timepoints_Pfs47_1,]
nrow(Matching_Mosq_Haplotypes_Pfs47_1) #14
matching_timepoints_Pfs47_2<-unique(intersect(Pfs47_FinalHaplotypes_Human$Timepoint[!is.na(Pfs47_FinalHaplotypes_Human$Pfs47_2_hap)], Pfs47_FinalHaplotypes_Mosq$Timepoint[!is.na(Pfs47_FinalHaplotypes_Mosq$Pfs47_2_hap)]))
Matching_Mosq_Haplotypes_Pfs47_2<- Pfs47_FinalHaplotypes_Mosq[Pfs47_FinalHaplotypes_Mosq$Timepoint %in% matching_timepoints_Pfs47_2,]
nrow(Matching_Mosq_Haplotypes_Pfs47_2) #29

#	Visualisation of timepoints at which there are both individuals and mosquitoes
p<-ggplot(Matching_Mosq_Haplotypes_CSP, aes(x=factor(Day, levels=c(0, 2, 7, 14, 21, 28, 35)), fill=Day)) +
  geom_bar() +
  labs(x="Day", y="Number of Occurrences - matching human/mosq") +
  theme_classic()+
  scale_fill_viridis_d() +
  theme(legend.position="none")+
  ggtitle("CSP")
q<-ggplot(Matching_Mosq_Haplotypes_TRAP, aes(x=factor(Day, levels=c(0, 2, 7, 14, 21, 28, 35)), fill=Day)) +
  geom_bar() +
  labs(x="Day", y="Number of Occurrences - matching human/mosq") +
  theme_classic()+
  scale_fill_viridis_d() +
  theme(legend.position="none")+
  ggtitle("TRAP")
grid.arrange(p, q, ncol=2)
combined_plot <- arrangeGrob(p, q, ncol=2)
ggsave ("Barchart_matching_hum_mosq.pdf", combined_plot,width = 14, height = 6)





# Percentage monoclonal versus polyclonal at each timepoint and each marker
MonoPoly<-FINAL_CSP_TRAP
MonoPoly$TRAP_MonoPoly[MonoPoly$TRAP_MOI==1]<-"Mono"
MonoPoly$TRAP_MonoPoly[MonoPoly$TRAP_MOI>1]<-"Poly"
MonoPoly$CSP_MonoPoly[MonoPoly$CSP_MOI==1]<-"Mono"
MonoPoly$CSP_MonoPoly[MonoPoly$CSP_MOI>1]<-"Poly"
MonoPoly$Species<-factor(MonoPoly$Species,levels = c("mosquito","human"))
MonoPoly$Mosquito<-NULL
CSP_MonoPoly <- MonoPoly[!is.na(MonoPoly$CSP_MonoPoly), ]
TRAP_MonoPoly <- MonoPoly[!is.na(MonoPoly$TRAP_MonoPoly), ]
pal<-c("lightgrey", "grey48")

# Calculate counts and positions
CSP_counts <- CSP_MonoPoly %>%
  group_by(Day, Species, CSP_MonoPoly) %>%
  dplyr::summarize(Count = n(), .groups = 'drop') %>%
  group_by(Day, Species) %>%
  mutate(LabelPos = cumsum(Count) - (0.5 * Count))


TRAP_counts <- TRAP_MonoPoly %>%
  group_by(Day, Species, TRAP_MonoPoly) %>%
  dplyr::summarize(Count = n(), .groups = 'drop') %>%
  group_by(Day, Species) %>%
  mutate(LabelPos = cumsum(Count) - (0.5 * Count))

a<-ggplot(CSP_MonoPoly, aes(x=as.factor(Day), fill=factor(CSP_MonoPoly))) + 
  geom_bar(position="fill", stat="count") +
  #geom_text(data = CSP_counts, aes(x = as.factor(Day), y = LabelPos/sum(Count), label = Count), position = position_fill(), vjust = 0.5) +
  xlab("Days after treatment initiation") +
  ylab("Percentage of samples") +
  ggtitle("CSP") +
  facet_wrap(~Species, ncol=1) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10)) +
  scale_fill_manual(values=pal)

b<-ggplot(TRAP_MonoPoly, aes(x=as.factor(Day),fill=factor(TRAP_MonoPoly))) + 
  geom_bar(position="fill", stat="count")+
  geom_text(data = TRAP_counts, aes(x = as.factor(Day), y = LabelPos/sum(Count), label = Count), position = position_fill(), vjust = 0.5) +
  xlab("Days after treatment initiation")+
  ylab("Percentage of samples")+
  ggtitle("TRAP")+
  facet_wrap(~Species,ncol=1)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values=pal)
combined_plot <- a + b + plot_layout(ncol=2, widths=c(1.1, 1))
ggsave("Fig1A_MonoPoly.pdf", combined_plot, width=14, height=6)


# Percentage monoclonal versus polyclonal at each timepoint - both markers combined
MonoPoly<-FINAL_CSP_TRAP
MonoPoly$TRAP_MonoPoly[MonoPoly$TRAP_MOI==1]<-"Mono"
MonoPoly$TRAP_MonoPoly[MonoPoly$TRAP_MOI>1]<-"Poly"
MonoPoly$CSP_MonoPoly[MonoPoly$CSP_MOI==1]<-"Mono"
MonoPoly$CSP_MonoPoly[MonoPoly$CSP_MOI>1]<-"Poly"
MonoPoly$Species<-factor(MonoPoly$Species,levels = c("mosquito","human"))
MonoPoly$Mosquito<-NULL
MonoPoly$MonoPoly<-"Mono"
MonoPoly$MonoPoly[MonoPoly$CSP_MonoPoly == "Poly" |MonoPoly$TRAP_MonoPoly == "Poly" ] <-"Poly"

pal<-c("lightgrey", "grey48")

# Calculate counts and positions
counts <- MonoPoly %>%
  group_by(Day, Species, MonoPoly) %>%
  dplyr::summarize(Count = n(), .groups = 'drop') %>%
  group_by(Day, Species) %>%
  mutate(LabelPos = cumsum(Count) - (0.5 * Count))

c<-ggplot(MonoPoly[MonoPoly$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(MonoPoly))) + 
  geom_bar(position="fill", stat="count")+
  #geom_text(data = counts, aes(x = as.factor(Day), y = LabelPos/sum(Count), label = Count), position = position_fill(), vjust = 1.5) +
  xlab("Days after treatment initiation")+
  ylab("Percentage of samples")+
  facet_wrap(~Species,ncol=1)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values=pal)
ggsave("Fig1A_MonoPoly_MarkersCombined.pdf", c, width=6, height=6)



# Melt data for further analysis
melt_data <- function(data) {
  data %>%
    separate_rows(Haplotype, Reads, Percentage, sep = " ")
}

# Melt and combine data
Haplotypes_CSP_Molten <- melt_data(Haplotypes_CSP)
Haplotypes_CSP_Molten$Timepoint <- gsub("_Mosq\\d+", "", Haplotypes_CSP_Molten$SampleID)
Haplotypes_TRAP_Molten <- melt_data(Haplotypes_TRAP)
Haplotypes_TRAP_Molten$Timepoint <- gsub("_Mosq\\d+", "", Haplotypes_TRAP_Molten$SampleID)
write.csv(Haplotypes_CSP_Molten,"Haplotypes_CSP_Molten.csv",row.names = FALSE)
Haplotypes_Pfs47_1_Molten <- melt_data(Haplotypes_Pfs47_1)
Haplotypes_Pfs47_1_Molten$Timepoint <- gsub("_Mosq\\d+", "", Haplotypes_Pfs47_1_Molten$SampleID)
Haplotypes_Pfs47_2_Molten <- melt_data(Haplotypes_Pfs47_2)
Haplotypes_Pfs47_2_Molten$Timepoint <- gsub("_Mosq\\d+", "", Haplotypes_Pfs47_2_Molten$SampleID)

# Only keep rows of timepoints at which there are both human and mosquito samples
Matching_Haplotypes_CSP_Molten <- Haplotypes_CSP_Molten %>%
  group_by(Timepoint) %>%
  filter("human" %in% species & "mosquito" %in% species) %>%
  ungroup()
Matching_Haplotypes_TRAP_Molten <- Haplotypes_TRAP_Molten %>%
  group_by(Timepoint) %>%
  filter("human" %in% species & "mosquito" %in% species) %>%
  ungroup()
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

Matching_Haplotypes_CSP_Molten <- Matching_Haplotypes_CSP_Molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
  ungroup()
Matching_Haplotypes_TRAP_Molten <- Matching_Haplotypes_TRAP_Molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
  ungroup()
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

# Plot Figure 1B
Matching_Haplotypes_CSP_Molten$Comparison<-factor(Matching_Haplotypes_CSP_Molten$Comparison,levels=c("human_only","mosquito_only","matching"))
pal<-c("#a891cf", "#f8997c","lightgrey")
a<- ggplot(Matching_Haplotypes_CSP_Molten[Matching_Haplotypes_CSP_Molten$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  xlab("Days after treatment initiation")+
  ylab("Percentage of reads")+
  ggtitle("CSP")+
  theme(legend.position="none",
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10))+
  scale_fill_manual(values=pal)
Matching_Haplotypes_TRAP_Molten$Comparison<-factor(Matching_Haplotypes_TRAP_Molten$Comparison,levels=c("human_only","mosquito_only","matching"))
b<- ggplot(Matching_Haplotypes_TRAP_Molten[Matching_Haplotypes_TRAP_Molten$Day%in% c("0","2","7","14"),], aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  #facet_wrap(~Matching_FinalHaplotypes_ALL$Host)
  xlab("Days after treatment initiation")+
  ylab("")+
  scale_fill_manual(values=pal, name =" ")+
  ggtitle("TRAP")+
  theme(axis.text.x = element_text(size = 10), 
      axis.text.y = element_text(size = 10))
combined_plot <- a + b + plot_layout(ncol=2, widths=c(1.1, 1))
ggsave("Fig1B_Matching_Haplotypes_Percentages.pdf", combined_plot, width=14, height=6)


#Check if matching haplotypes in multiclonal individuals more often have high or low coverage at Day0? CSP
Matching_Haplotypes_CSP_Molten_Day0<-Matching_Haplotypes_CSP_Molten[Matching_Haplotypes_CSP_Molten$Day==0,]
Matching_Haplotypes_CSP_Molten_Day0<-Matching_Haplotypes_CSP_Molten_Day0[Matching_Haplotypes_CSP_Molten_Day0$MOI>2,]
Matching_Haplotypes_CSP_Molten_Day0<-Matching_Haplotypes_CSP_Molten_Day0[Matching_Haplotypes_CSP_Molten_Day0$species=="human",]
human_only_df <- Matching_Haplotypes_CSP_Molten_Day0 %>%
  filter(Comparison == "human_only") %>%
  mutate(Percentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison) %>%
  summarise(TotalPercentage = sum(Percentage, na.rm = TRUE))

# Filtering the non-human_only rows
non_human_only_df <- Matching_Haplotypes_CSP_Molten_Day0 %>%
  filter(Comparison != "human_only")%>%
  mutate(TotalPercentage = as.numeric(Percentage))

# Binding the two dataframes together
filtered_df <- bind_rows(human_only_df, non_human_only_df)

ggplot(filtered_df, aes(x=Comparison,y=as.numeric(TotalPercentage))) + 
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 3) +
  theme_classic()+
  xlab("Type")+
  ylab("Percentage of total reads")+
  ggtitle("CSP")
ggsave("Percentage_Reads_Human_only_vs_matching_CSP.pdf", width=6, height=5)

#Check if matching haplotypes in multiclonal individuals more often have high or low coverage at Day0 per GC category. CSP
Matching_Haplotypes_CSP_Molten_Day0_GC_Cat<-merge(Matching_Haplotypes_CSP_Molten_Day0,ClinicalData,by="SampleID")
Matching_Haplotypes_CSP_Molten_Day0_GC_Cat$GC_cat[Matching_Haplotypes_CSP_Molten_Day0_GC_Cat$totalgct_ul>100]<-2
Matching_Haplotypes_CSP_Molten_Day0_GC_Cat$GC_cat[Matching_Haplotypes_CSP_Molten_Day0_GC_Cat$totalgct_ul<100]<-1
human_only_df <- Matching_Haplotypes_CSP_Molten_Day0_GC_Cat %>%
  filter(Comparison == "human_only") %>%
  mutate(Percentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison) %>%
  summarise(
    TotalPercentage = sum(Percentage, na.rm = TRUE),
    GC_cat = mean(GC_cat, na.rm = TRUE)
  )

# Filtering the non-human_only rows
non_human_only_df <- Matching_Haplotypes_CSP_Molten_Day0_GC_Cat %>%
  filter(Comparison != "human_only")%>%
  mutate(TotalPercentage = as.numeric(Percentage))

# Binding the two dataframes together
filtered_df <- bind_rows(human_only_df, non_human_only_df)

ggplot(filtered_df, aes(x=Comparison,y=as.numeric(TotalPercentage))) + 
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 3) +
  facet_wrap(~GC_cat)+
  theme_classic()+
  xlab("Type")+
  ylab("Percentage of total reads")+
  ggtitle("CSP")
ggsave("Percentage_Reads_Human_only_vs_matching_CSP_GC_cat.pdf", width=6, height=5)


#Check if matching haplotypes in multiclonal individuals more often have high or low coverage at Day0? TRAP
Matching_Haplotypes_TRAP_Molten_Day0<-Matching_Haplotypes_TRAP_Molten_matching[Matching_Haplotypes_TRAP_Molten$Day==0,]
Matching_Haplotypes_TRAP_Molten_Day0<-Matching_Haplotypes_TRAP_Molten_matching_Day0[Matching_Haplotypes_TRAP_Molten_Day0$MOI>2,]
Matching_Haplotypes_TRAP_Molten_Day0<-Matching_Haplotypes_TRAP_Molten_matching_Day0[Matching_Haplotypes_TRAP_Molten_Day0$species=="human",]
human_only_df <- Matching_Haplotypes_TRAP_Molten_Day0 %>%
  filter(Comparison == "human_only") %>%
  mutate(Percentage = as.numeric(Percentage)) %>%
  group_by(SampleID, Comparison) %>%
  summarise(TotalPercentage = sum(Percentage, na.rm = TRUE))

# Filtering the non-human_only rows
non_human_only_df <- Matching_Haplotypes_TRAP_Molten_Day0 %>%
  filter(Comparison != "human_only")%>%
  mutate(TotalPercentage = as.numeric(Percentage))

# Binding the two dataframes together
filtered_df <- bind_rows(human_only_df, non_human_only_df)

ggplot(filtered_df, aes(x=Comparison,y=as.numeric(TotalPercentage))) + 
  geom_boxplot()+
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 3) +
  theme_classic()+
  xlab("Type")+
  ylab("Percentage of total reads")+
  ggtitle("TRAP")
ggsave("Percentage_Reads_Human_only_vs_matching_TRAP.pdf", width=6, height=5)


#Check correlation between CSP MOI and TRAP MOI
ggplot(FINAL_CSP_TRAP, aes(x=CSP_MOI, y=TRAP_MOI)) + 
  geom_jitter()+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("CSP MOI")+
  ylab("TRAP MOI")+
  xlim(0,10)+
  geom_smooth(method=lm)+
  ggtitle("Correlation between CSP And TRAP MOI")
ggsave("Correlation_CSP_TRAP.pdf", width=6, height=5)


cor.test(FINAL_CSP_TRAP$CSP_MOI, FINAL_CSP_TRAP$TRAP_MOI, method = c("pearson")) #0.745256

#Heatmap
ggplot(FINAL_CSP_TRAP, aes(x = CSP_MOI, y = TRAP_MOI)) +
  geom_tile(stat = "bin_2d", aes(fill = ..count..), bins = 10) +  # Adjust 'bins' as needed
  scale_fill_gradientn(colours = c("yellow", "orange", "red", "blue","purple","green"),  # Add more colors as needed
                       values = scales::rescale(c(0, 0.3, 0.6, 1)),  # Adjust these values as needed
                       limits = c(0, 120)) +
  theme_classic() +
  xlab("CSP MOI") +
  ylab("TRAP MOI") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10))+
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10))+
  ggtitle("Correlation between CSP And TRAP MOI")

###### Haplotype prevalence in human vs mosquito samples

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

# Calculate fractions for TRAP_hap and CSP_hap
trap_results <- calculate_haplotype_fractions(FINAL_CSP_TRAP, "TRAP_hap")
csp_results <- calculate_haplotype_fractions(FINAL_CSP_TRAP, "CSP_hap")

# Combine the results
combined_results <- bind_rows(trap_results, csp_results)


final_results <- combined_results %>%
  group_by(haplotype) %>%
  summarise(frac_sample_human = sum(frac_sample[Species == "human"], na.rm = TRUE),
            frac_sample_mosquito = sum(frac_sample[Species == "mosquito"], na.rm = TRUE),
            .groups = 'drop')
final_results$marker <- sub("-.*", "", final_results$haplotype)
final_results$difference <- final_results$frac_sample_human-final_results$frac_sample_mosquito

final_results<- final_results%>%  
mutate(hap_lab = ifelse(haplotype == 'csp-22', "RTS,S/AS01\nvaccine haplotype", NA),
       hap_lab2 = ifelse(haplotype == 'trap-8', "3D7 trap", NA))
  
ggplot(final_results,aes(x = frac_sample_human, y = frac_sample_mosquito,col=difference)) +
  geom_abline(intercept = 0, col = 'grey80') +
  geom_point(aes(shape = marker), size = 3) +
  scale_shape_manual(values = c(17, 19)) +
  #geom_text_repel(aes(label = hap_lab), color = 'black') + 
  #geom_text_repel(aes(label = hap_lab2), color = 'black') +
  geom_text_repel(aes(label = haplotype), color = 'black') +
  labs(x = 'Haplotype prevalence in human samples', y = 'Haplotype prevalence in mosquito samples',
       shape = 'Marker', col = 'Difference in prevalence') +
  coord_cartesian(xlim = c(0,40), ylim = c(0,40)) + 
  scale_colour_gradient2(low = "red",mid = "white",high = "darkblue",midpoint = 0)+
  theme_light()
ggsave("Haplotype_prevalence_humans_mosquitoes.pdf", width=7, height=5)



#Check clonality over time in the individuals and mosquitoes
FINAL_CSP_TRAP<-FINAL_CSP_TRAP %>%
  mutate(MOI_Combined = pmax(FINAL_CSP_TRAP$TRAP_MOI,FINAL_CSP_TRAP$CSP_MOI, na.rm = TRUE))

c<-ggplot(FINAL_CSP_TRAP[FINAL_CSP_TRAP$Species == "human" & FINAL_CSP_TRAP$Day%in% c("0","2","7","14","21","28"),], aes(x=Day, y=MOI_Combined)) + 
  geom_boxplot(aes(fill=Day), alpha = 0.5, fill="lightgrey", outlier.shape = NA) +
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Days after treatment initiation", y ="MOI")+
  scale_y_continuous(breaks = seq(0, 10, by = 1)) + 
  ggtitle("Clonality in individuals") 

d <- ggplot(FINAL_CSP_TRAP[FINAL_CSP_TRAP$Species == "mosquito",], aes(x=Day, y=MOI_Combined)) + 
  geom_boxplot(alpha = 0.5, fill="lightgrey", outlier.shape = NA) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x= "Days after treatment initiation", y ="MOI") +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) + 
  ggtitle("Clonality in mosquitoes")

combined_clonality_plot <- c + d + plot_layout(ncol=2)
ggsave("Clonality.pdf", combined_clonality_plot, width=10, height=5)


#Presence of clones in population
#All haplotypes at any timepoint
Total_haplotypes_CSP <- unique(na.omit(str_split(FinalHaplotypes$CSP_hap, " ") %>% unlist())) #57
ordered_Total_haplotypes_CSP <- Total_haplotypes_CSP[order(as.integer(str_extract(Total_haplotypes_CSP, "\\d+")))]
print(ordered_Total_haplotypes_CSP)
Total_haplotypes_TRAP <- unique(na.omit(str_split(FinalHaplotypes$TRAP_hap, " ") %>% unlist())) #53
ordered_Total_haplotypes_TRAP <- Total_haplotypes_TRAP[order(as.integer(str_extract(Total_haplotypes_TRAP, "\\d+")))]
print(ordered_Total_haplotypes_TRAP)
#human
#TRAP
FinalHaplotypes_Human_Day0<-FinalHaplotypes_Human[FinalHaplotypes_Human$Day=="0",]
all_haplotypes <- na.omit(str_split(FinalHaplotypes_Human_Day0$TRAP_hap, " ") %>% unlist())
haplotypes <- unique(all_haplotypes)

# Calculate percentages
percentage_list <- sapply(haplotypes, function(hap) {
  detected <- str_detect(FinalHaplotypes_Human_Day0$TRAP_hap, regex(hap, ignore_case = TRUE))
  occurrences <- sum(detected, na.rm = TRUE)
  total_rows <- sum(!is.na(detected))
  (occurrences / total_rows) * 100
})

percentage_df_TRAP_human <- data.frame(
  Haplotype = names(percentage_list),
  Percentage = as.numeric(percentage_list)
)
sorted_df_TRAP_human <- percentage_df_TRAP_human %>%
  mutate(SortOrder = as.integer(str_extract(Haplotype, "\\d+"))) %>%
  arrange(SortOrder) %>%
  select(-SortOrder)
sorted_df_TRAP_human$Haplotype <- factor(sorted_df_TRAP_human$Haplotype, levels = sorted_df_TRAP_human$Haplotype)

a<-ggplot(sorted_df_TRAP_human, aes(x = Haplotype, y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distribution of TRAP Haplotypes in human baseline population",
    x = "Haplotype",
    y = "Percentage"
  ) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#CSP
all_haplotypes <- na.omit(str_split(FinalHaplotypes_Human_Day0$CSP_hap, " ") %>% unlist())
haplotypes <- unique(all_haplotypes)

# Calculate percentages
percentage_list <- sapply(haplotypes, function(hap) {
  detected <- str_detect(FinalHaplotypes_Human_Day0$CSP_hap, regex(hap, ignore_case = TRUE))
  occurrences <- sum(detected, na.rm = TRUE)
  total_rows <- sum(!is.na(detected))
  (occurrences / total_rows) * 100
})

percentage_df_CSP_human <- data.frame(
  Haplotype = names(percentage_list),
  Percentage = as.numeric(percentage_list)
)
sorted_df_CSP_human <- percentage_df_CSP_human %>%
  mutate(SortOrder = as.integer(str_extract(Haplotype, "\\d+"))) %>%
  arrange(SortOrder) %>%
  select(-SortOrder)
sorted_df_CSP_human$Haplotype <- factor(sorted_df_CSP_human$Haplotype, levels = sorted_df_CSP_human$Haplotype)

b<-ggplot(sorted_df_CSP_human, aes(x = Haplotype, y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distribution of CSP Haplotypes in human baseline population",
    x = "Haplotype",
    y = "Percentage"
  ) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
 

#mosquito
#TRAP
FinalHaplotypes_Mosq_Day0<-FinalHaplotypes_Mosq[FinalHaplotypes_Mosq$Day=="0",]
all_haplotypes <- na.omit(str_split(FinalHaplotypes_Mosq_Day0$TRAP_hap, " ") %>% unlist())
haplotypes <- unique(all_haplotypes)

# Calculate percentages
percentage_list <- sapply(haplotypes, function(hap) {
  detected <- str_detect(FinalHaplotypes_Mosq_Day0$TRAP_hap, regex(hap, ignore_case = TRUE))
  occurrences <- sum(detected, na.rm = TRUE)
  total_rows <- sum(!is.na(detected))
  (occurrences / total_rows) * 100
})

percentage_df_TRAP_mosq <- data.frame(
  Haplotype = names(percentage_list),
  Percentage = as.numeric(percentage_list)
)
sorted_df_TRAP_mosq <- percentage_df_TRAP_mosq %>%
  mutate(SortOrder = as.integer(str_extract(Haplotype, "\\d+"))) %>%
  arrange(SortOrder) %>%
  select(-SortOrder)
sorted_df_TRAP_mosq$Haplotype <- factor(sorted_df_TRAP_mosq$Haplotype, levels = sorted_df_TRAP_mosq$Haplotype)

c<-ggplot(sorted_df_TRAP_mosq, aes(x = Haplotype, y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distribution of TRAP Haplotypes in mosquito baseline population",
    x = "Haplotype",
    y = "Percentage"
  ) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#CSP
all_haplotypes <- na.omit(str_split(FinalHaplotypes_Mosq_Day0$CSP_hap, " ") %>% unlist())
haplotypes <- unique(all_haplotypes)

# Calculate percentages
percentage_list <- sapply(haplotypes, function(hap) {
  detected <- str_detect(FinalHaplotypes_Mosq_Day0$CSP_hap, regex(hap, ignore_case = TRUE))
  occurrences <- sum(detected, na.rm = TRUE)
  total_rows <- sum(!is.na(detected))
  (occurrences / total_rows) * 100
})

percentage_df_CSP_mosq <- data.frame(
  Haplotype = names(percentage_list),
  Percentage = as.numeric(percentage_list)
)
sorted_df_CSP_mosq <- percentage_df_CSP_mosq %>%
  mutate(SortOrder = as.integer(str_extract(Haplotype, "\\d+"))) %>%
  arrange(SortOrder) %>%
  select(-SortOrder)
sorted_df_CSP_mosq$Haplotype <- factor(sorted_df_CSP_mosq$Haplotype, levels = sorted_df_CSP_mosq$Haplotype)

d<-ggplot(sorted_df_CSP_mosq, aes(x = Haplotype, y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distribution of CSP Haplotypes in mosquito baseline population",
    x = "Haplotype",
    y = "Percentage"
  ) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

combined_distribution_plots_haplotyptes_Day0 <- b + d + a +c + plot_layout(ncol=2)
combined_distribution_plots_haplotyptes_Day0
ggsave("combined_distribution_plots_haplotyptes_Day0.pdf", combined_distribution_plots_haplotyptes_Day0, width=15, height=10)


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


#Check relation between clonality and infection rates (individuals only)
ClinicalData<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Paper_AmpSeq/ClinicalData.csv",header=TRUE,sep=',')
ClinicalData$SampleID<-paste(ClinicalData$studycode,"_Day",ClinicalData$studyvisit,sep="")

merge<-merge(FINAL_CSP_TRAP,ClinicalData,by="SampleID")
merge$percentagemosqinfected<-(merge$mosq_pos/merge$mosq_total)*100

ggplot(merge, aes(x=as.factor(MOI_Combined), y=percentagemosqinfected,fill=factor(MOI_Combined))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(MOI_Combined)), alpha=0.9)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Infected Mosquitoes (%)")+
  ggtitle("Infection rate ~ clonality")
ggsave("Infection_rate_clonality.pdf", width=6, height=5)

#Check relation between clonality and infection rates (individuals only) - normalised by GC
merge$normalized_infection_rate <- (merge$percentagemosqinfected / merge$totalgct_ul) * 100

ggplot(merge, aes(x=as.factor(MOI_Combined), y=normalized_infection_rate,fill=factor(MOI_Combined))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(MOI_Combined)), alpha=0.9)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Infected Mosquitoes (%)")+
  ggtitle("Infection rate ~ clonality")
ggsave("Infection_rate_clonality_Normalised_by_GC.pdf", width=6, height=5)


#Check relation between clonality and gametocytemia
ggplot(merge, aes(x=as.factor(MOI_Combined), y=totalgct_ul,fill=factor(MOI_Combined))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(MOI_Combined)), alpha=0.9)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Gametocyte density")+
  scale_y_log10()+
  ggtitle("Clonality ~ Gametocytaemia - all timepoints")
ggsave("Clonaty_GC.pdf", width=6, height=5)

#For Day0 only
ggplot(merge[merge$Day=="0",], aes(x=as.factor(MOI_Combined), y=totalgct_ul,fill=factor(MOI_Combined))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(MOI_Combined)), alpha=0.9)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Gametocyte density")+
  scale_y_log10()+
  ggtitle("Clonality ~ Gametocytaemia")

#Check relation between clonality and age (individuals only)
merge <- merge %>%
  mutate(age_group = case_when(
    age < 5 ~ "<5",
    age >= 5 & age < 10 ~ "5-10",
    age >= 10 & age < 15 ~ "10-15",
    age >= 15 & age < 20 ~ "15-20",
    age >= 20 ~ ">20"
  ),
  age_group = factor(age_group, levels = c("<5", "5-10", "10-15", "15-20", ">20")))

merge <- merge %>%
  group_by(age) %>%
  mutate(mean_MOI_by_age = mean(MOI_Combined, na.rm = TRUE)) %>%
  ungroup()

ggplot(merge, aes(x=age, y=mean_MOI_by_age)) + 
  geom_line()+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Age",y="Mean MOI")+
  ggtitle("Mean clonality by age")
ggsave("Clonality_by_age.pdf", width=6, height=5)


#Check relation between midgut clonality and oocyst number (mosquitoes only)
OocystData<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Paper_AmpSeq/Oocyst_data.csv",header=TRUE,sep=',')
OocystData$SampleID<-paste(OocystData$studycode,"_Day",OocystData$studyvisit,sep="")


transformed_df <- data.frame()

for (i in 1:nrow(OocystData)) {
  row <- OocystData[i,]
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
      transformed_df <- rbind(transformed_df, new_row)
    }
  }
}


merge<-merge(FINAL_CSP_TRAP,transformed_df,by="SampleID")
ggplot(merge, aes(x=as.factor(MOI_Combined), y=oocysts,fill=factor(MOI_Combined))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(MOI_Combined)), alpha=0.9)+
  theme(legend.position = "none")+
  labs(x= "MOI",y="Oocyst number")+
  ggtitle("Oocyst number ~ clonality")
ggsave("Oocyst_clonality.pdf", width=6, height=5)

#Check relation between individual clonality and oocyst number (mosquitoes only)
transformed_df$Individual_ID <- gsub("(.*)_(.*)", "\\1", transformed_df$SampleID)
FINAL_CSP_TRAP_human<-FINAL_CSP_TRAP[FINAL_CSP_TRAP$Species == "human",]
merge <- merge(x = FINAL_CSP_TRAP_human, y = transformed_df, by.x = 'SampleID', by.y = 'Individual_ID', all.y = TRUE)
merge <- merge[!is.na(merge$MOI_Combined), ]
  
ggplot(merge, aes(x=as.factor(MOI_Combined), y=oocysts,fill=factor(MOI_Combined))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(MOI_Combined)), alpha=0.9)+
  theme(legend.position = "none")+
  labs(x= "MOI individual",y="Oocyst number")+
  ggtitle("Individual clonality ~ Oocyst number")
ggsave("Individual_Clonality_Oocyst_number.pdf", width=6, height=5)

#Check relation between individual clonality and oocyst number (mosquitoes only) - normalised by GC
merge<-merge(merge,ClinicalData,by="SampleID")
merge$normalized_oocysts <- (merge$oocysts / merge$totalgct_ul) * 100

ggplot(merge, aes(x=as.factor(MOI_Combined), y=normalized_oocysts,fill=factor(MOI_Combined))) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
  theme_classic()+
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(MOI_Combined)), alpha=0.9)+
  theme(legend.position = "none")+
  labs(x= "MOI individual",y="Oocyst number")+
  ggtitle("Individual clonality ~ Oocyst number")
ggsave("Individual_Clonality_Oocyst_number_normalised_by_GC.pdf", width=6, height=5)

#Check relation between clonality and month (individuals only)
merge$visitdate_hb <- as.Date(merge$visitdate_hb, format="%d/%m/%Y")  # Adjust the format argument as needed
merge$visitdate_hb <- as.numeric(merge$visitdate_hb)

merge <- merge %>%
  group_by(visitdate_hb) %>%
  mutate(mean_MOI_by_time = mean(MOI_Combined, na.rm = TRUE)) %>%
  ungroup()

ggplot(merge, aes(x=visitdate_hb, y=mean_MOI_by_time)) + 
  geom_line()+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Visitdate",y="Mean MOI")+
  ggtitle("Mean clonality over time (individuals only)")
ggsave("Clonality_by_age.pdf", width=6, height=5)




#Check accross haplotypes which ones are matching or not
CSP_grouped_df <- Matching_Haplotypes_CSP_Molten %>%
  group_by(Haplotype, Comparison) %>%
  summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
CSP_grouped_df <- CSP_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("csp-", "", Haplotype))) %>%
  arrange(Haplotype_num)

TRAP_grouped_df <- Matching_Haplotypes_TRAP_Molten %>%
  group_by(Haplotype, Comparison) %>%
  summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
TRAP_grouped_df <- TRAP_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("trap-", "", Haplotype))) %>%
  arrange(Haplotype_num)

# Plot
e<-ggplot(CSP_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#9AC77B', "grey", '#7EC2DE'))+
  labs(title = "CSP",
       x = "Haplotype",
       y = "Count") +
  theme_light()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
f<-ggplot(TRAP_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#9AC77B', "grey", '#7EC2DE'))+
  labs(title = "TRAP",
       x = "Haplotype",
       y = "Count") +
  theme_light()+
  theme(axis.text.x = element_blank())
combined_plot <- e / f 
ggsave("Haplotypes_transmission.pdf", combined_plot, width=5, height=7)

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

# Overview of all haplotypes
Matching_Haplotypes_TRAP_Molten <- Matching_Haplotypes_TRAP_Molten %>%
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

Matching_Haplotypes_TRAP_Molten$Day_Species <- factor(Matching_Haplotypes_TRAP_Molten$Day_Species, levels = generate_ordered_levels(Matching_Haplotypes_TRAP_Molten$Day_Species))

p <- ggplot(Matching_Haplotypes_TRAP_Molten, aes(x=Day_Species, y=Percentage, fill=Haplotype)) + 
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

ggsave("TRAP_overview.pdf", p, width=40, height=50,limitsize = FALSE)
write.csv(Matching_Haplotypes_TRAP_Molten,"Matching_Haplotypes_TRAP_Molten.csv")

# Overview of all haplotypes
Matching_Haplotypes_CSP_Molten <- Matching_Haplotypes_CSP_Molten %>%
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

Matching_Haplotypes_CSP_Molten$Day_Species <- factor(Matching_Haplotypes_CSP_Molten$Day_Species, levels = generate_ordered_levels(Matching_Haplotypes_CSP_Molten$Day_Species))

p <- ggplot(Matching_Haplotypes_CSP_Molten, aes(x=Day_Species, y=Percentage, fill=Haplotype)) + 
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

ggsave("CSP_overview.pdf", p, width=40, height=50,limitsize = FALSE)
write.csv(Matching_Haplotypes_CSP_Molten,"Matching_Haplotypes_CSP_Molten.csv")

#Make nice plot for CSP PQ-04-043
Matching_Haplotypes_CSP_Molten_043<-Matching_Haplotypes_CSP_Molten[Matching_Haplotypes_CSP_Molten$Individual=="PQ-04-043",]

Matching_Haplotypes_CSP_Molten_043$Day <- as.numeric(as.character(Matching_Haplotypes_CSP_Molten_043$Day))
Matching_Haplotypes_CSP_Molten_043$Percentage <- as.numeric(as.character(Matching_Haplotypes_CSP_Molten_043$Percentage))

Matching_Haplotypes_CSP_Molten_043 <- Matching_Haplotypes_CSP_Molten_043 %>%
  complete(Day, nesting(Haplotype, species), fill = list(Reads = 0))

Matching_Haplotypes_CSP_Molten_043$Percentage <- ifelse(
  is.na(Matching_Haplotypes_CSP_Molten_043$Percentage),
  0,
  Matching_Haplotypes_CSP_Molten_043$Percentage
)

Matching_Haplotypes_CSP_Molten_043$SampleID[Matching_Haplotypes_CSP_Molten_043$species == "human"] <- paste0("PQ-04-043_Day",Matching_Haplotypes_CSP_Molten_043$Day[Matching_Haplotypes_CSP_Molten_043$species == "human"])
merge <- merge(Matching_Haplotypes_CSP_Molten_043, ClinicalData, by = "SampleID")
merge$totalPC<-merge$ring_ul_all + merge$totalgct_ul
merge$Perc_totalPC<- (merge$totalPC * merge$Percentage)/100

cols <- unikn::usecol(c("lavender", "rosybrown1", "indianred3", "palegreen3", "plum4", "cyan3"), n = 9)
#cols<- c("#FBB4AE","#FFF2AE","#C2E699" ,"#FDCDAC", "#CCEBC5" ,"#FDDAEC" ,"#FEEBE2", "#CBD5E8","#F1E2CC","#EDF8E9")
unique_haplotypes <- unique(merge$Haplotype)
names(cols) <- unique_haplotypes 

ordered_haplotypes <- merge %>%
  filter(Day == 0) %>%
  arrange(Percentage) %>%
  pull(Haplotype)
merge$Haplotype <- factor(merge$Haplotype, levels = unique(ordered_haplotypes))

area_plot<-ggplot(merge, aes(x = Day, y = Percentage, fill = Haplotype)) +
  facet_wrap(~species, scales = "free_y") + 
  geom_area(position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(0,2,7,21))+
  theme_classic() +
  theme(legend.position = "bottom")

#ggplot(merge, aes(x = Day, y = Perc_totalPC, fill = Haplotype)) +
#facet_wrap(~species, scales = "free_y") + 
#geom_area(position = "stack") +
# scale_fill_viridis_d() +
# labs(x = "Day", y = "Reads") +
# theme_classic() +
# theme(legend.position = "bottom") +
# scale_y_sqrt()

Matching_Haplotypes_CSP_Molten_043 <- Matching_Haplotypes_CSP_Molten_043 %>%
  mutate(mosquito_id = gsub(".*_", "", SampleID))  

bar_plot <- ggplot(Matching_Haplotypes_CSP_Molten_043[Matching_Haplotypes_CSP_Molten_043$species == "mosquito",], 
                   aes(x = as.factor(Day), y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~mosquito_id, scales = "free_x", strip.position = "bottom")


combined_plot <- grid.arrange(bar_plot, area_plot, ncol=1)
ggsave("PQ-04-043.pdf", combined_plot, width=10, height=10,limitsize = FALSE)

ggplot(merge, aes(x=Day, y=ring_ul_all))+ 
  geom_line()+
  geom_point()+
  labs(x="Timepoint",y="rings/uL")+
  scale_x_continuous(breaks=c(0,2,7,21))+
  theme_classic()
ggsave("PQ04043_PC.pdf", width=18, height=1.4)

ggplot(merge, aes(x=Day, y=totalgct_ul))+ 
  geom_line()+
  geom_point()+
  labs(x="Timepoint",y="Gametocytes/uL")+
  scale_x_continuous(breaks=c(0,2,7,21))+
  theme_classic()
ggsave("PQ04043_GC.pdf", width=18, height=1.4)


#Make nice plot for CSP PQ-04-026
Matching_Haplotypes_CSP_Molten_026<-Matching_Haplotypes_CSP_Molten[Matching_Haplotypes_CSP_Molten$Individual=="PQ-04-026",]

Matching_Haplotypes_CSP_Molten_026$Day <- as.numeric(as.character(Matching_Haplotypes_CSP_Molten_026$Day))
Matching_Haplotypes_CSP_Molten_026$Percentage <- as.numeric(as.character(Matching_Haplotypes_CSP_Molten_026$Percentage))

Matching_Haplotypes_CSP_Molten_026 <- Matching_Haplotypes_CSP_Molten_026 %>%
  complete(Day, nesting(Haplotype, species), fill = list(Reads = 0))

Matching_Haplotypes_CSP_Molten_026$Percentage <- ifelse(
  is.na(Matching_Haplotypes_CSP_Molten_026$Percentage),
  0,
  Matching_Haplotypes_CSP_Molten_026$Percentage
)

Matching_Haplotypes_CSP_Molten_026$SampleID[Matching_Haplotypes_CSP_Molten_026$species == "human"] <- paste0("PQ-04-026_Day",Matching_Haplotypes_CSP_Molten_026$Day[Matching_Haplotypes_CSP_Molten_026$species == "human"])
merge <- merge(Matching_Haplotypes_CSP_Molten_026, ClinicalData, by = "SampleID")
merge$totalPC<-merge$ring_ul_all + merge$totalgct_ul
merge$Perc_totalPC<- (merge$totalPC * merge$Percentage)/100

cols <- unikn::usecol(c("lavender", "khaki2", "lightskyblue2", "plum4"), n = 3)

unique_haplotypes <- unique(merge$Haplotype)
names(cols) <- unique_haplotypes 

ordered_haplotypes <- merge %>%
  filter(Day == 0) %>%
  arrange(Percentage) %>%
  pull(Haplotype)
merge$Haplotype <- factor(merge$Haplotype, levels = unique(ordered_haplotypes))

area_plot<-ggplot(merge, aes(x = Day, y = Percentage, fill = Haplotype)) +
  facet_wrap(~species, scales = "free_y") + 
  geom_area(position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(0,2,7,21))+
  theme_classic() +
  theme(legend.position = "bottom")

#ggplot(merge, aes(x = Day, y = Perc_totalPC, fill = Haplotype)) +
#facet_wrap(~species, scales = "free_y") + 
#geom_area(position = "stack") +
# scale_fill_viridis_d() +
# labs(x = "Day", y = "Reads") +
# theme_classic() +
# theme(legend.position = "bottom") +
# scale_y_sqrt()

Matching_Haplotypes_CSP_Molten_026 <- Matching_Haplotypes_CSP_Molten_026 %>%
  mutate(mosquito_id = gsub(".*_", "", SampleID))  

bar_plot <- ggplot(Matching_Haplotypes_CSP_Molten_026[Matching_Haplotypes_CSP_Molten_026$species == "mosquito",], 
                   aes(x = as.factor(Day), y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~mosquito_id, scales = "free_x", strip.position = "bottom")


combined_plot <- grid.arrange(bar_plot, area_plot, ncol=1)
ggsave("PQ-04-026.pdf", combined_plot, width=10, height=10,limitsize = FALSE)

ggplot(merge, aes(x=Day, y=ring_ul_all))+ 
  geom_line()+
  geom_point()+
  labs(x="Timepoint",y="rings/uL")+
  scale_x_continuous(breaks=c(0,2,7,14,21,28))+
  theme_classic()
ggsave("PQ04026_PC.pdf", width=18, height=1.4)

ggplot(merge, aes(x=Day, y=totalgct_ul))+ 
  geom_line()+
  geom_point()+
  labs(x="Timepoint",y="Gametocytes/uL")+
  scale_x_continuous(breaks=c(0,2,7,14,21,28))+
  theme_classic()
ggsave("PQ04026_GC.pdf", width=18, height=1.4)


#Make nice plot for TRAP PQ-04-043
Matching_Haplotypes_TRAP_Molten_043<-Matching_Haplotypes_TRAP_Molten[Matching_Haplotypes_TRAP_Molten$Individual=="PQ-04-043",]

Matching_Haplotypes_TRAP_Molten_043$Day <- as.numeric(as.character(Matching_Haplotypes_TRAP_Molten_043$Day))
Matching_Haplotypes_TRAP_Molten_043$Percentage <- as.numeric(as.character(Matching_Haplotypes_TRAP_Molten_043$Percentage))

Matching_Haplotypes_TRAP_Molten_043 <- Matching_Haplotypes_TRAP_Molten_043 %>%
  complete(Day, nesting(Haplotype, species), fill = list(Reads = 0))

Matching_Haplotypes_TRAP_Molten_043$Percentage <- ifelse(
  is.na(Matching_Haplotypes_TRAP_Molten_043$Percentage),
  0,
  Matching_Haplotypes_TRAP_Molten_043$Percentage
)

Matching_Haplotypes_TRAP_Molten_043$SampleID[Matching_Haplotypes_TRAP_Molten_043$species == "human"] <- paste0("PQ-04-043_Day",Matching_Haplotypes_TRAP_Molten_043$Day[Matching_Haplotypes_TRAP_Molten_043$species == "human"])
merge <- merge(Matching_Haplotypes_TRAP_Molten_043, ClinicalData, by = "SampleID")
merge$totalPC<-merge$ring_ul_all + merge$totalgct_ul
merge$Perc_totalPC<- (merge$totalPC * merge$Percentage)/100

cols <- unikn::usecol(c("goldenrod2", "indianred3", "rosybrown1", "palegreen3", "plum4", "dodgerblue"), n = 11)

unique_haplotypes <- unique(merge$Haplotype)
names(cols) <- unique_haplotypes 

ordered_haplotypes <- merge %>%
  filter(Day == 0) %>%
  arrange(Percentage) %>%
  pull(Haplotype)
merge$Haplotype <- factor(merge$Haplotype, levels = unique(ordered_haplotypes))

area_plot<-ggplot(merge, aes(x = Day, y = Percentage, fill = Haplotype)) +
  facet_wrap(~species, scales = "free_y") + 
  geom_area(position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(0,2,7,21))+
  theme_classic() +
  theme(legend.position = "bottom")

#ggplot(merge, aes(x = Day, y = Perc_totalPC, fill = Haplotype)) +
#facet_wrap(~species, scales = "free_y") + 
#geom_area(position = "stack") +
# scale_fill_viridis_d() +
# labs(x = "Day", y = "Reads") +
# theme_classic() +
# theme(legend.position = "bottom") +
# scale_y_sqrt()

Matching_Haplotypes_TRAP_Molten_043 <- Matching_Haplotypes_TRAP_Molten_043 %>%
  mutate(mosquito_id = gsub(".*_", "", SampleID))  

bar_plot <- ggplot(Matching_Haplotypes_TRAP_Molten_043[Matching_Haplotypes_TRAP_Molten_043$species == "mosquito",], 
                   aes(x = as.factor(Day), y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~mosquito_id, scales = "free_x", strip.position = "bottom")


combined_plot <- grid.arrange(bar_plot, area_plot, ncol=1)
ggsave("PQ-04-043_TRAP.pdf", combined_plot, width=10, height=10,limitsize = FALSE)


#Make nice plot for CSP PQ-04-047
Matching_Haplotypes_CSP_Molten_047<-Matching_Haplotypes_CSP_Molten[Matching_Haplotypes_CSP_Molten$Individual=="PQ-04-047",]

Matching_Haplotypes_CSP_Molten_047$Day <- as.numeric(as.character(Matching_Haplotypes_CSP_Molten_047$Day))
Matching_Haplotypes_CSP_Molten_047$Percentage <- as.numeric(as.character(Matching_Haplotypes_CSP_Molten_047$Percentage))

Matching_Haplotypes_CSP_Molten_047 <- Matching_Haplotypes_CSP_Molten_047 %>%
  complete(Day, nesting(Haplotype, species), fill = list(Reads = 0))

Matching_Haplotypes_CSP_Molten_047$Percentage <- ifelse(
  is.na(Matching_Haplotypes_CSP_Molten_047$Percentage),
  0,
  Matching_Haplotypes_CSP_Molten_047$Percentage
)

Matching_Haplotypes_CSP_Molten_047$SampleID[Matching_Haplotypes_CSP_Molten_047$species == "human"] <- paste0("PQ-04-047_Day",Matching_Haplotypes_CSP_Molten_047$Day[Matching_Haplotypes_CSP_Molten_047$species == "human"])
merge <- merge(Matching_Haplotypes_CSP_Molten_047, ClinicalData, by = "SampleID")
merge$totalPC<-merge$ring_ul_all + merge$totalgct_ul
merge$Perc_totalPC<- (merge$totalPC * merge$Percentage)/100

cols <- unikn::usecol(c("goldenrod2", "indianred3", "rosybrown1", "palegreen3", "plum4", "dodgerblue"), n = 11)

unique_haplotypes <- unique(merge$Haplotype)
names(cols) <- unique_haplotypes 

ordered_haplotypes <- merge %>%
  filter(Day == 0) %>%
  arrange(Percentage) %>%
  pull(Haplotype)
merge$Haplotype <- factor(merge$Haplotype, levels = unique(ordered_haplotypes))

area_plot<-ggplot(merge, aes(x = Day, y = Percentage, fill = Haplotype)) +
  facet_wrap(~species, scales = "free_y") + 
  geom_area(position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(0,2,7,21))+
  theme_classic() +
  theme(legend.position = "bottom")

#ggplot(merge, aes(x = Day, y = Perc_totalPC, fill = Haplotype)) +
#facet_wrap(~species, scales = "free_y") + 
#geom_area(position = "stack") +
# scale_fill_viridis_d() +
# labs(x = "Day", y = "Reads") +
# theme_classic() +
# theme(legend.position = "bottom") +
# scale_y_sqrt()

Matching_Haplotypes_CSP_Molten_047 <- Matching_Haplotypes_CSP_Molten_047 %>%
  mutate(mosquito_id = gsub(".*_", "", SampleID))  

bar_plot <- ggplot(Matching_Haplotypes_CSP_Molten_047[Matching_Haplotypes_CSP_Molten_047$species == "mosquito",], 
                   aes(x = as.factor(Day), y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Day", y = "Reads") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~mosquito_id, scales = "free_x", strip.position = "bottom")


combined_plot <- grid.arrange(bar_plot, area_plot, ncol=1)
ggsave("PQ-04-047.pdf", combined_plot, width=10, height=10,limitsize = FALSE)

#Check whether multiclonal samples have higher likelihood of causing multiclonal mosquito infections
human_df <- FINAL_CSP_TRAP %>% filter(Species == "human")
mosquito_df <- FINAL_CSP_TRAP %>% filter(Species == "mosquito")
joined_df <- inner_join(human_df, mosquito_df, by = "Timepoint", suffix = c("_human", "_mosquito"))

ggplot(joined_df, aes(x = MOI_Combined_human, y = MOI_Combined_mosquito)) +
  geom_jitter() +
  labs(
    title = "Comparison of MOI_combined in Humans and Mosquitoes",
    x = "MOI_combined in Humans",
    y = "MOI_combined in Mosquitoes"
  ) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  theme_classic()
ggsave("Clonality_human_vs_mosq.pdf", width=6, height=5)

cor.test(joined_df$MOI_Combined_human, joined_df$MOI_Combined_mosquito, method = c("pearson")) #0.5278912 


#Make overview table
OverviewTable<-ClinicalData
OverviewTable$GC<-as.numeric(OverviewTable$femalegct_ul)+as.numeric(OverviewTable$malegct_ul)
OverviewTable<-OverviewTable[-which(OverviewTable$studyvisit==1),]
#hist(OverviewTable$GC, main="Histogram of GC", xlab="GC", col="lightblue", border="black")
#plot(density(na.omit(OverviewTable$GC)), main="Density Plot of GC", xlab="GC")
OverviewTable$GC_cat[OverviewTable$GC>100]<-2
OverviewTable$GC_cat[OverviewTable$GC<100]<-1
#hist(OverviewTable$ring_ul_all, main="Histogram of ring_ul_all", xlab="ring_ul_all", col="lightblue", border="black")
#plot(density(na.omit(OverviewTable$ring_ul_all)), main="Density Plot of ring_ul_all", xlab="ring_ul_all")
OverviewTable$Asex_cat[OverviewTable$ring_ul_all>100]<-2
OverviewTable$Asex_cat[OverviewTable$ring_ul_all<100]<-1
OverviewTable <- OverviewTable %>%
  left_join(select(Matching_Haplotypes_CSP_Molten, SampleID, MOI), by = "SampleID")
OverviewTable <- OverviewTable %>%
  left_join(select(Matching_Haplotypes_TRAP_Molten, SampleID, MOI), by = "SampleID")
colnames(OverviewTable)[c(17:18)]<-c("MOI_CSP","MOI_TRAP")

# First summarize Matching_Haplotypes_CSP_Molten_matching
summary_table_CSP <- Matching_Haplotypes_CSP_Molten %>% 
  group_by(SampleID) %>%
  summarise(
    count_matching = sum(Comparison == "matching"),
    count_human_only = sum(Comparison == "human_only"),
    count_mosquito_only = sum(Comparison == "mosquito_only")
  )

# Now join this summary_table with OverviewTable
OverviewTable <- OverviewTable %>%
  left_join(summary_table_CSP, by = "SampleID")

# First summarize Matching_Haplotypes_TRAP_Molten_matching
summary_table_TRAP <- Matching_Haplotypes_TRAP_Molten_matching %>% 
  group_by(SampleID) %>%
  summarise(
    count_matching = sum(Comparison == "matching"),
    count_human_only = sum(Comparison == "human_only"),
    count_mosquito_only = sum(Comparison == "mosquito_only")
  )

# Now join this summary_table with OverviewTable
OverviewTable <- OverviewTable %>%
  left_join(summary_table_TRAP, by = "SampleID")

colnames(OverviewTable)[c(16:21)]<-c("count_matching_CSP","count_human_only_CSP","count_mosquito_only_CSP","count_matching_TRAP","count_human_only_TRAP","count_mosquito_only_TRAP")

OverviewTable_filtered <- OverviewTable %>%
  group_by(studycode) %>%
  filter(!(all(is.na(MOI_CSP)) & all(is.na(MOI_TRAP)))) %>%
  ungroup()

#Individuals non-infectious at D0
Individuals_non_inf_Day0<-unique(OverviewTable_filtered$studycode[OverviewTable_filtered$mosq_pos==0 & OverviewTable_filtered$studyvisit=="0"])
Matching_Haplotypes_CSP_Molten$Infectious_D0<- "Yes"
Matching_Haplotypes_CSP_Molten$Infectious_D0[Matching_Haplotypes_CSP_Molten$Individual %in% Individuals_non_inf_Day0]<- "No"

#Haplotypes in individuals non-infectious at D0
ggplot(Matching_Haplotypes_CSP_Molten[Matching_Haplotypes_CSP_Molten$species=="human",], aes(x = Haplotype, y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Haplotype in infectious vs non-infectious individuals D0",
    x = "Haplotype",
    y = "Percentage"
  ) +
  facet_wrap(~Infectious_D0)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("Haplotype in infectious vs non-infectious individuals D0.pdf", width=6, height=5)

#Individuals non-infectious at certain timepoint
Matching_Haplotypes_CSP_Molten_infectious <- Matching_Haplotypes_CSP_Molten %>%
  left_join(ClinicalData, by = "SampleID")
Matching_Haplotypes_CSP_Molten_infectious$Infectious<- 1
Matching_Haplotypes_CSP_Molten_infectious$Infectious[Matching_Haplotypes_CSP_Molten_infectious$mosq_pos == 0]<- 0

ggplot(Matching_Haplotypes_CSP_Molten_infectious[Matching_Haplotypes_CSP_Molten_infectious$species=="human",], aes(x = Haplotype, y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Haplotype in infectious vs non-infectious individuals at all timepoints",
    x = "Haplotype",
    y = "Percentage"
  ) +
  facet_wrap(~Infectious)+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("Haplotype in infectious vs non-infectious individuals all timepoints.pdf", width=6, height=5)

#model <- glm(Infectious ~ Haplotype + totalgct_ul, data = Matching_Haplotypes_CSP_Molten_infectious, family = binomial())
#summary(model)

model <- logistf(Infectious ~ Haplotype + totalgct_ul, 
                 data = Matching_Haplotypes_CSP_Molten_infectious)
summary(model)


coefficients <- coef(summary(model))
df <- data.frame(
  term = rownames(coefficients),
  estimate = coefficients[, "Estimate"],
  lower = coefficients[, "Estimate"] - 1.96 * coefficients[, "Std. Error"],
  upper = coefficients[, "Estimate"] + 1.96 * coefficients[, "Std. Error"]
)

# Remove the intercept for the plot
df <- df[df$term != "(Intercept)", ]

# Create forest plot using ggplot2
ggplot(df, aes(x=estimate, y=term)) + 
  geom_point() +
  geom_segment(aes(x=lower, xend=upper, y=term, yend=term)) +
  geom_vline(xintercept=0, linetype="dashed") +
  labs(title="Forest Plot", x="Effect Size", y="")+
  theme_classic()


#To analyze if the likelihood that a certain haplotype is "matching" depends on the haplotype environment (other haplotypes present in the same SampleID) - logistic regression
Matching_Haplotypes_CSP_Molten_infectious$IsMatching <- ifelse(Matching_Haplotypes_CSP_Molten_infectious$Comparison == "matching", 1, 0)
wide_df <- dcast(Matching_Haplotypes_CSP_Molten_infectious, SampleID ~ Haplotype, value.var = "IsMatching", fun.aggregate = max)
wide_df[wide_df == -Inf] <- NA

# Fit a logistic regression model
model <- glm(IsMatching ~ Hap1 + Hap2 + Hap3, data = wide_df, family = "binomial")

# Print the summary of the model
summary(model)



######  Combine ClinicalData with FINAL_CSP_TRAP
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

Haplotypes_CSP_Molten_ <- Haplotypes_CSP_Molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
  ungroup()

Clin_CSP <- Haplotypes_CSP_Molten_ %>%
  left_join(select(ClinicalData, SampleID, totalgct_ul, ring_ul_all), by = "SampleID")

#write.csv(Clin_CSP_TRAP,"Clin_CSP_TRAP.csv", row.names = FALSE)

# Normalize read percentages in human data
human_data_CSP <- Clin_CSP %>% 
  filter(species != 'mosquito') %>%
  mutate(Normalized_CSP_Perc = as.numeric(Percentage) / as.numeric(totalgct_ul))

# Prepare mosquito data
mosquito_data_CSP <- Clin_CSP %>% 
  filter(species == 'mosquito')

# Merge human and mosquito data
merged_data_CSP <- merge(human_data_CSP, mosquito_data_CSP, by = c("Timepoint", "Haplotype"), all = TRUE)

merged_data_CSP$Normalized_CSP_Perc <- ifelse(is.na(merged_data_CSP$Normalized_CSP_Perc), 0, merged_data_CSP$Normalized_CSP_Perc)
merged_data_CSP$Percentage.y <- ifelse(is.na(merged_data_CSP$Percentage.y), 0, merged_data_CSP$Percentage.y)


# Add a column to indicate if the read percentage is higher in mosquito sample
merged_data_CSP$Higher_in_Mosquito <- merged_data_CSP$Percentage.y > merged_data_CSP$Normalized_CSP_Perc

# Calculate the odds of transmission for each haplotype
haplotype_odds <- merged_data_CSP %>%
  group_by(Haplotype) %>%
  summarize(Odds_of_Transmission = mean(Higher_in_Mosquito, na.rm = TRUE))

haplotype_odds$Haplotype_Num <- as.numeric(sub("csp-", "", haplotype_odds$Haplotype))

# Ordering the dataframe by this numeric part
haplotype_odds <- haplotype_odds[order(haplotype_odds$Haplotype_Num), ]

# Removing the auxiliary numeric column if not needed further
haplotype_odds$Haplotype_Num <- NULL

y_axis_data <- data.frame(Y = seq(0, 1, by = 0.1), Odds_of_Transmission = 0)

# Create the circular plot with the y-axis circles
ggplot(haplotype_odds, aes(x = factor(Haplotype), y = Odds_of_Transmission, fill = Odds_of_Transmission)) +
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

#what the odds are of transmission for each haplotype in a certain clonal environment (eg less likely if other haplotypes are present as well at that certain individual and timepoint)

#Keep only individuals and timepoints that have csp-1
Timepoint_csp1<-merged_data_CSP$Timepoint[merged_data_CSP$Haplotype=="csp-1"]
merged_data_CSP<-merged_data_CSP[merged_data_CSP$Timepoint%in%Timepoint_csp1,]
merged_data_CSP$Transmission<-ifelse(merged_data_CSP$Comparison.x=="matching",1,0)
merged_data_CSP$Transmission <- ifelse(is.na(merged_data_CSP$Transmission), 0, merged_data_CSP$Transmission)
#merged_data_CSP$Transmission <- ifelse(merged_data_CSP$Haplotype=="csp-1", merged_data_CSP$Transmission, NA)

# Remove haplotypes that have less than 5 occurences
# Count the occurrences of each level in the Haplotype column
haplotype_counts <- table(merged_data_CSP$Haplotype)

# Identify the names of haplotypes that occur at least 5 times
haplotypes_to_keep <- names(haplotype_counts[haplotype_counts >= 5])

# Subset the dataframe to keep only the rows with those haplotypes
merged_data_CSP <- merged_data_CSP[merged_data_CSP$Haplotype %in% haplotypes_to_keep, ]

# Run the logistic regression model with transmission =1 or 0
model <- glm(Transmission ~ MOI.x + factor(Haplotype) + totalgct_ul.x, 
             family = binomial, data = merged_data_CSP)

# Summary of the model
summary(model)
#write.csv(merged_data_CSP,"merged_data_CSP.csv", row.names = FALSE)

anova(model)

# Getting model summary
model_summary <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)

# Creating the data frame for plotting
plot_data <- model_summary %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(term = gsub("factor\\(Haplotype\\)", "", term)) # Clean up term names

# Filter out the extreme values by excluding any estimates above a certain threshold
threshold <- 10^10
plot_data_filtered <- subset(plot_data, estimate < threshold)

# Convert estimates and confidence intervals to log scale for the plot
plot_data_filtered$log_estimate <- log(plot_data_filtered$estimate)
plot_data_filtered$log_conf.low <- log(plot_data_filtered$conf.low)
plot_data_filtered$log_conf.high <- log(plot_data_filtered$conf.high)
plot_data_filtered<-plot_data_filtered[-which(plot_data_filtered$term=="csp-10"),]

# Creating the forest plot
ggplot(plot_data_filtered, aes(y = term, x = log_estimate)) +
  geom_point() +
  geom_errorbarh(aes(xmin = log_conf.low, xmax = log_conf.high), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Log Odds Ratio") +
  ylab("") +
  scale_x_continuous(limits=c(-10,20))+
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


#########TRAP
Haplotypes_TRAP_Molten <- Haplotypes_TRAP_Molten %>%
  group_by(Timepoint, Haplotype) %>%
  mutate(Comparison = determine_sample_type(cur_data()),
         Day = str_extract(Timepoint, "_([^_]+)$"),  
         Day = str_replace_all(Day, "_", ""), 
         Day = as.numeric(str_sub(Day, 4)),
         Day = factor(Day, levels = c(0, 2, 7, 14, 21, 28, 35))) %>%
  ungroup()

Clin_TRAP <- Haplotypes_TRAP_Molten %>%
  left_join(select(ClinicalData, SampleID, totalgct_ul, ring_ul_all), by = "SampleID")

#write.csv(Clin_CSP_TRAP,"Clin_CSP_TRAP.csv", row.names = FALSE)

# Normalize read percentages in human data
human_data_TRAP <- Clin_TRAP %>% 
  filter(species != 'mosquito') %>%
  mutate(Normalized_TRAP_Perc = as.numeric(Percentage) / as.numeric(totalgct_ul))

# Prepare mosquito data
mosquito_data_TRAP <- Clin_TRAP %>% 
  filter(species == 'mosquito')

# Merge human and mosquito data
merged_data_TRAP <- merge(human_data_TRAP, mosquito_data_TRAP, by = c("Timepoint", "Haplotype"), all = TRUE)

merged_data_TRAP$Normalized_TRAP_Perc <- ifelse(is.na(merged_data_TRAP$Normalized_TRAP_Perc), 0, merged_data_TRAP$Normalized_TRAP_Perc)
merged_data_TRAP$Percentage.y <- ifelse(is.na(merged_data_TRAP$Percentage.y), 0, merged_data_TRAP$Percentage.y)


# Add a column to indicate if the read percentage is higher in mosquito sample
merged_data_TRAP$Higher_in_Mosquito <- merged_data_TRAP$Percentage.y > merged_data_TRAP$Normalized_TRAP_Perc

# Calculate the odds of transmission for each haplotype
haplotype_odds <- merged_data_TRAP %>%
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
  geom_point(data = y_axis_data, aes(x = 0, y = Y), shape = 1, size = 4, color = "gray") + 
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
