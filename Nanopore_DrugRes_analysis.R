require(dplyr)
require(plyr)
require(scales)
require(data.table)
require(ggplot2)
library("ggpattern")


setwd("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION")
Variants_Exp51<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/variants_Exp51.csv")
Variants_Exp52<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/variants_Exp52.csv")

#Coverage
ggplot(Variants_Exp51, aes(x = as.factor(gene), y = depth)) +
  geom_boxplot()+
  scale_y_sqrt(breaks = c(1, 50, 100, 500, 1000, 30000, 60000,90000))+
  theme_classic()

ggplot(Variants_Exp52, aes(x = as.factor(gene), y = depth)) +
  geom_boxplot()+
  scale_y_sqrt(breaks = c(1, 50, 100, 500, 1000, 30000, 60000,90000))+
  theme_classic()

Variants_Exp51$mosqid<-sub("^([^_]+_[^_]+_)(.*)", "\\2", Variants_Exp51$sample_id)
Variants_Exp51$uniqueid<-sub("^([^_]+).*", "\\1",  Variants_Exp51$mosqid)
Variants_Exp52$uniqueid<-sub("^([^_]+_[^_]+_)(.*)", "\\2", Variants_Exp52$sample_id)


#Overall frequency of SNPs
Variants_Exp52 <- Variants_Exp52 %>%
  group_by(genome_pos, gene, protein_change) %>%
  mutate(species="human",
         mosqid=NA)
Variants_Exp52 <- Variants_Exp52 %>%
  filter(!grepl("fs", protein_change))
Variants_Exp52<-Variants_Exp52[,c(1:10,13,11:12)]

Variants_Exp51 <- Variants_Exp51 %>%
  group_by(genome_pos, gene, protein_change) %>%
  mutate(species="mosquito")
Variants_Exp51 <- Variants_Exp51 %>%
  filter(!grepl("fs", protein_change))

both<-rbind(Variants_Exp51,Variants_Exp52)

# Create a complete list of combinations of sample_id and protein_change
complete_combinations <- expand.grid(sample_id = unique(both$sample_id), protein_change = unique(both$protein_change))

# Merge with original data
merged_data <- left_join(complete_combinations, both, by = c("sample_id", "protein_change"))

# Step 3: Replace NA in freq with zero
merged_data$freq[is.na(merged_data$freq)] <- 0

gene_protein_mapping <- merged_data %>%
  select(protein_change, gene) %>%
  distinct() %>%
  na.omit()

merged_data <- merged_data %>%
  left_join(gene_protein_mapping, by = "protein_change", suffix = c("", "_ref"))

merged_data <- merged_data %>%
  mutate(gene = ifelse(is.na(gene), gene_ref, gene)) %>%
  select(-gene_ref)

sampleid_species_mapping <- merged_data %>%
  select(sample_id, species) %>%
  distinct() %>%
  na.omit()

merged_data <- merged_data %>%
  left_join(sampleid_species_mapping, by = "sample_id", suffix = c("", "_ref"))

merged_data <- merged_data %>%
  mutate(species = ifelse(is.na(species), species_ref, species)) %>%
  select(-species_ref)

merged_data$protein_change <- factor(merged_data$protein_change,levels = c("p.Lys76Thr","p.Tyr184Phe","p.Gly102Gly","p.Asn86Tyr","p.Asn51Ile","p.Cys59Arg","p.Ser108Asn","p.Ile431Val","p.Ser436Ala","p.Ser436Phe","p.Ala437Gly","p.Lys540Glu","p.Ala581Gly","p.Ala613Ser"))
merged_data <- merged_data %>%
  mutate(protein_change = recode(protein_change,
                                 "p.Lys76Thr" = "CRT-Lys76Thr",
                                 "p.Tyr184Phe" = "MDR1-Tyr184Phe",
                                 "p.Gly102Gly" = "MDR1-Gly102Gly",
                                 "p.Asn86Tyr" = "MDR1-Asn86Tyr",
                                 "p.Asn51Ile" = "DHFR-Asn51Ile",
                                 "p.Cys59Arg" = "DHFR-Cys59Arg",
                                 "p.Ser108Asn" = "DHFR-Ser108Asn",
                                 "p.Ile431Val" = "DHPS-Ile431Val",
                                 "p.Ser436Ala" = "DHPS-Ser436Ala",
                                 "p.Ser436Phe" = "DHPS-Ser436Phe",
                                 "p.Ala437Gly" = "DHPS-Ala437Gly",
                                 "p.Lys540Glu" = "DHPS-Lys540Glu",
                                 "p.Ala581Gly" = "DHPS-Ala581Gly",
                                 "p.Ala613Ser" = "DHPS-Ala613Ser"))

summary_data <- merged_data %>%
  group_by(protein_change,gene,species) %>%
  summarise(
    mean_freq = mean(freq, na.rm = TRUE),
    se = sd(freq, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )


mycols <- c("#8dd3c6", "#f98073", "#fccde5", "#80b0d3")
mycols2 <- c("#4c9184", "#a8392d", "#bf6b96", "#3373a1")
ggplot(summary_data, aes(x = as.factor(protein_change), y = mean_freq, fill = gene, pattern = species,color =gene)) +
  geom_bar_pattern(aes(pattern_fill=gene,pattern_color=gene),position = "dodge", stat = "identity", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  geom_errorbar(aes(ymin = mean_freq - se, ymax = mean_freq + se,color = gene), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Protein change", y = "Overall frequency in entire sample set", fill = "Gene") +
  scale_fill_manual(values = mycols) +  
  scale_pattern_manual(values = c("none", "stripe"),labels=c("Human", "Mosquito")) +
  scale_color_manual(values = mycols2) +
  scale_pattern_fill_manual(values = mycols2) +
  scale_pattern_color_manual(values = mycols2)+
  theme_light()+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),limits=c(0,1))
ggsave("Overall_frequencies_bar.pdf",width=20, height=8)



#Fisher exact test
count_df1 <- table(Variants_Exp51$protein_change)
count_df2 <- table(Variants_Exp52$protein_change)

# Get unique protein_changes from both dataframes
unique_protein_changes <- unique(c(Variants_Exp51$protein_change, Variants_Exp52$protein_change))

# Initialize a list to store Fisher's exact test results
results_list <- list()

# Loop through each unique protein_change and perform Fisher's exact test
for (change in unique_protein_changes) {
  # Count uniqueids for the specific protein_change in each dataframe
  count_change_df1 <- sum(Variants_Exp51$protein_change == change)
  count_change_df2 <- sum(Variants_Exp52$protein_change == change)
  
  # Total uniqueids in each dataframe
  total_uniqueids_df1 <- length(unique(Variants_Exp51$mosqid))
  total_uniqueids_df2 <- length(unique(Variants_Exp52$uniqueid))
  
  # Create a 2x2 contingency table if counts are non-negative
  if (count_change_df1 >= 0 && count_change_df2 >= 0) {
    table_change <- matrix(c(count_change_df1, total_uniqueids_df1 - count_change_df1,
                             count_change_df2, total_uniqueids_df2 - count_change_df2),
                           nrow = 2, byrow = TRUE)
  
    
    # Perform Fisher's exact test
    result <- tryCatch(fisher.test(table_change), error = function(e) NA)
    
    # Store test result for each protein_change in the list if test is successful
      results_list[[change]] <- result
  }
}

# Display the results for Fisher's exact test
results_list




#Count number of samples
# individuals in human data
length(unique(Variants_Exp52$uniqueid)) #48
# individuals in mosq data
length(unique(Variants_Exp51$uniqueid)) #32
# individuals for which there are matching mosquitoes
length(intersect(Variants_Exp51$uniqueid, Variants_Exp52$uniqueid)) #30
# mosquitoes in mosq data
length(unique(Variants_Exp51$mosqid)) #73
#	Matching cases individuals-mosquitoes
length(unique(Variants_Exp51$mosqid[Variants_Exp51$uniqueid%in%Variants_Exp52$uniqueid])) #68


#Difference in frequency (human-mosquito) for matching samples
human_uniqueid<-unique(Variants_Exp52$uniqueid)
mosq_uniqueid<-unique(Variants_Exp51$uniqueid)
Variants_Hum_matching<-Variants_Exp52[Variants_Exp52$uniqueid %in% mosq_uniqueid,]
Variants_Mosq_matching<-Variants_Exp51[Variants_Exp51$uniqueid %in% human_uniqueid,]
Variants_Hum_matching<-Variants_Hum_matching[-which(Variants_Hum_matching$type=="frameshift_variant"),]
Variants_Mosq_matching<-Variants_Mosq_matching[-which(Variants_Mosq_matching$type=="frameshift_variant"),]
Variants_Hum_matching$species<-"human"
Variants_Mosq_matching$species<-"mosquito"


#write.csv(Variants_Hum_matching,"Variants_Hum_matching.csv",row.names=FALSE)
#write.csv(Variants_Mosq_matching,"Variants_Mosq_matching.csv",row.names=FALSE)
#write.csv(merged_data,"merged_data.csv",row.names=FALSE)



hum_df <- Variants_Hum_matching
mosq_df <- Variants_Mosq_matching

# Unique mosqid's
unique_mosqids <- unique(mosq_df$mosqid)
number_of_unique_mosqids <- length(unique_mosqids)


# Initialize an empty dataframe to store the results
final_hum_df <- data.frame()

# Duplicate rows in Hum dataframe based on the number of mosqids and assign a mosqid to each set of rows
for (mosqid in unique_mosqids) {
  # Create a temporary dataframe with duplicated rows and the mosqid column
  temp_df <- hum_df
  temp_df$mosqid <- mosqid
  
  # Append the temporary dataframe to the final dataframe
  final_hum_df <- rbind(final_hum_df, temp_df)
}


# Perform left joins to combine 'human' and 'mosquito' data on uniqueid and genome_pos
merged_data <- merge(Variants_Mosq_matching, final_hum_df, by = c('mosqid', 'uniqueid','genome_pos','protein_change','gene'), suffixes = c('_mosquito', '_human'), all = TRUE)

# Replace NA values with 0 in freq columns after merge
merged_data[is.na(merged_data)] <- 0

# Calculate the frequency differences
merged_data$freq_diff <- as.numeric(merged_data$freq_human - merged_data$freq_mosquito)

# Select necessary columns for the output dataframe
output_df <- merged_data[c('uniqueid', 'genome_pos', 'freq_diff','protein_change','gene','mosqid')]


# Calculate mean and standard error for each protein_change
protein_summary <- output_df %>%
  group_by(protein_change,gene) %>%
  summarise(
    mean_freq_diff = mean(freq_diff, na.rm = TRUE),
    se = sd(freq_diff, na.rm = TRUE) / sqrt(n()),
    count = n()
  )

# Create the bar plot with error bars
mycols <- c("#8dd3c6", "#f98073", "#fccde5", "#80b0d3")
ggplot(protein_summary, aes(x = protein_change, y = mean_freq_diff,fill=gene)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_freq_diff - se, ymax = mean_freq_diff + se), width = .2, position = position_dodge(.9)) +
  geom_text(aes(label = count, y = 0.05),
            vjust = 0, position = position_dodge(0.9)) +
  scale_fill_manual(values=mycols)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Protein Change", y = "Mean Freq Diff", title = "Mean Freq Diff by Protein Change with Observation Count")+
  theme_light()
ggsave("Pairwise_comparison.pdf",width=15, height=10)


# Create a mapping of uniqueid to mosqid
uniqueid_to_mosqid <- unique(mosq_df[, c('uniqueid', 'mosqid')])

# Initialize an empty dataframe to store the results
final_hum_df <- data.frame()

# Duplicate rows in Hum dataframe based on the mapping
for (row in 1:nrow(uniqueid_to_mosqid)) {
  uniqueid <- uniqueid_to_mosqid$uniqueid[row]
  mosqid <- uniqueid_to_mosqid$mosqid[row]
  
  # Create a temporary dataframe with duplicated rows for the current uniqueid and the mosqid column
  temp_df <- hum_df[hum_df$uniqueid == uniqueid, ]
  temp_df$mosqid <- mosqid
  
  # Append the temporary dataframe to the final dataframe
  final_hum_df <- rbind(final_hum_df, temp_df)
}

# Perform the merge with the corrected duplication
merged_data <- merge(mosq_df, final_hum_df, by = c('mosqid', 'uniqueid', 'genome_pos', 'protein_change', 'gene'), suffixes = c('_mosquito', '_human'), all = TRUE)

# Replace NA values with 0 in freq columns after merge
merged_data[is.na(merged_data)] <- 0

# Calculate the frequency differences
merged_data$freq_diff <- as.numeric(merged_data$freq_human - merged_data$freq_mosquito)

# Select necessary columns for the output dataframe
merged_data$protein_change <- factor(merged_data$protein_change,levels = c("p.Lys76Thr","p.Tyr184Phe","p.Gly102Gly","p.Asn86Tyr","p.Asn51Ile","p.Cys59Arg","p.Ser108Asn","p.Ile431Val","p.Ser436Ala","p.Ala437Gly","p.Lys540Glu","p.Ala581Gly","p.Ala613Ser"))
output_df <- merged_data[c('uniqueid', 'genome_pos', 'freq_diff', 'protein_change', 'gene', 'mosqid')]
                                        
plotme<-ggplot(output_df, aes(x = as.factor(protein_change), y = freq_diff*100,fill=as.factor(gene))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(x = "Genome Position", y = "Difference in frequency (human-mosquito)", fill = "Gene") +
  scale_y_continuous(limits=c(-100,100))+
  facet_wrap(~uniqueid+mosqid,scales="free_x")+
  theme_classic()+
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Difference_Frequency_matching.pdf",plotme , width=14, height=30)
write.csv(merged_data, "merged_data.csv")


 
