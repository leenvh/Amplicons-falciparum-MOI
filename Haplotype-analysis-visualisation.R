library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(ggpubr)
library(viridis)

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

#	Visualisation of matching timepoints individuals-mosquitoes
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


# Melt data for further analysis
Matching_Haplotypes_CSP$species <- ifelse(grepl("Mosq", Matching_Haplotypes_CSP$SampleID), "mosquito", "human")
Matching_Haplotypes_TRAP$species <- ifelse(grepl("Mosq", Matching_Haplotypes_TRAP$SampleID), "mosquito", "human")
Matching_Haplotypes_Pfs47_1$species <- ifelse(grepl("Mosq", Matching_Haplotypes_Pfs47_1$SampleID), "mosquito", "human")
Matching_Haplotypes_Pfs47_2$species <- ifelse(grepl("Mosq", Matching_Haplotypes_Pfs47_2$SampleID), "mosquito", "human")



melt_data <- function(data) {
  data %>%
    separate_rows(Haplotype, Reads, Percentage, sep = " ")
}

# Melt and combine data
Matching_Haplotypes_CSP_Molten <- melt_data(Matching_Haplotypes_CSP)
Matching_Haplotypes_CSP_Molten$Timepoint <- gsub("_Mosq\\d+", "", Matching_Haplotypes_CSP_Molten$SampleID)
Matching_Haplotypes_TRAP_Molten <- melt_data(Matching_Haplotypes_TRAP)
Matching_Haplotypes_TRAP_Molten$Timepoint <- gsub("_Mosq\\d+", "", Matching_Haplotypes_TRAP_Molten$SampleID)
write.csv(Matching_Haplotypes_CSP_Molten,"Matching_Haplotypes_CSP_Molten.csv",row.names = FALSE)
Matching_Haplotypes_Pfs47_1_Molten <- melt_data(Matching_Haplotypes_Pfs47_1)
Matching_Haplotypes_Pfs47_1_Molten$Timepoint <- gsub("_Mosq\\d+", "", Matching_Haplotypes_Pfs47_1_Molten$SampleID)
Matching_Haplotypes_Pfs47_2_Molten <- melt_data(Matching_Haplotypes_Pfs47_2)
Matching_Haplotypes_Pfs47_2_Molten$Timepoint <- gsub("_Mosq\\d+", "", Matching_Haplotypes_Pfs47_2_Molten$SampleID)

# Only keep rows of timepoints at which there are both human and mosquito samples
Matching_Haplotypes_CSP_Molten <- Matching_Haplotypes_CSP_Molten %>%
  group_by(Timepoint) %>%
  filter("human" %in% species & "mosquito" %in% species) %>%
  ungroup()
Matching_Haplotypes_TRAP_Molten <- Matching_Haplotypes_TRAP_Molten %>%
  group_by(Timepoint) %>%
  filter("human" %in% species & "mosquito" %in% species) %>%
  ungroup()
Matching_Haplotypes_Pfs47_1_Molten <- Matching_Haplotypes_Pfs47_1_Molten %>%
  group_by(Timepoint) %>%
  filter("human" %in% species & "mosquito" %in% species) %>%
  ungroup()
Matching_Haplotypes_Pfs47_2_Molten <- Matching_Haplotypes_Pfs47_2_Molten %>%
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

# Plot
Matching_Haplotypes_CSP_Molten$Comparison<-factor(Matching_Haplotypes_CSP_Molten$Comparison,levels=c("human_only","mosquito_only","matching"))
pal<-c("#2c84a5", "#6f4e7c","#d4d084")
a<- ggplot(Matching_Haplotypes_CSP_Molten, aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  xlab("Day")+
  ylab("")+
  theme(legend.position="none")+
  ggtitle("CSP")+
  scale_fill_manual(values=pal)
Matching_Haplotypes_TRAP_Molten$Comparison<-factor(Matching_Haplotypes_TRAP_Molten$Comparison,levels=c("human_only","mosquito_only","matching"))
b<- ggplot(Matching_Haplotypes_TRAP_Molten, aes(x=as.factor(Day),fill=factor(Comparison))) + 
  geom_bar(position="fill", stat="count")+
  theme_classic()+
  #facet_wrap(~Matching_FinalHaplotypes_ALL$Host)
  xlab("Day")+
  ylab("")+
  scale_fill_manual(values=pal, name =" ")+
  ggtitle("TRAP")
combined_plot <- a + b + plot_layout(ncol=2, widths=c(1.1, 1))
ggsave("Matching_Haplotypes_Percentages.pdf", combined_plot, width=14, height=6)



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



#Check clonality over time in the individuals and mosquitoes
FINAL_CSP_TRAP<-FINAL_CSP_TRAP %>%
  mutate(MOI_Combined = pmax(FINAL_CSP_TRAP$TRAP_MOI,FINAL_CSP_TRAP$CSP_MOI, na.rm = TRUE))

c<-ggplot(FINAL_CSP_TRAP[FINAL_CSP_TRAP$Species == "human",], aes(x=Day, y=MOI_Combined,fill=Day)) + 
  geom_boxplot(aes(fill=Day), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(position=position_jitter(width=0.4, height=0.2), aes(colour=Day), alpha=0.9)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Visit Day", y ="MOI")+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  ggtitle("Clonality in individuals")

d<-ggplot(FINAL_CSP_TRAP[FINAL_CSP_TRAP$Species == "mosquito",], aes(x=Day, y=MOI_Combined,fill=Day)) + 
  geom_boxplot(aes(fill=Day), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(position=position_jitter(width=0.4, height=0.2), aes(colour=Day), alpha=0.9)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x= "Visit Day", y ="MOI")+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  ggtitle("Clonality in mosquitoes")

combined_clonality_plot <- c + d + plot_layout(ncol=2)
ggsave("Clonality.pdf", combined_clonality_plot, width=10, height=5)


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
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(title = "CSP",
       x = "Haplotype",
       y = "Count") +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
f<-ggplot(TRAP_grouped_df, aes(x = reorder(Haplotype, Haplotype_num), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  labs(title = "TRAP",
       x = "Haplotype",
       y = "Count") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
combined_plot <- e + f + plot_annotation(title ="Are some haplotypes more likely to transmit than others?") 
ggsave("Haplotypes_transmission.pdf", combined_plot, width=11, height=5)

#Check accross Pfs47 haplotypes which ones are matching or not
Pfs47_1_grouped_df <- Matching_Haplotypes_Pfs47_1_Molten %>%
  group_by(Haplotype, Comparison) %>%
  summarise(Count = n()) %>%
  ungroup()

# Extract the number from "csp-" to order the haplotypes
Pfs47_1_grouped_df <- Pfs47_1_grouped_df %>%
  mutate(Haplotype_num = as.numeric(gsub("pfs47_1-", "", Haplotype))) %>%
  arrange(Haplotype_num)

Pfs47_2_grouped_df <- Matching_Haplotypes_Pfs47_2_Molten %>%
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
