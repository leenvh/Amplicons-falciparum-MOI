require(dplyr)
require(plyr)
require(scales)
require(data.table)
require(ggplot2)

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

length(unique(Variants_Exp52$sample_id)) #48
#Calculate percentages
Variants_Exp52_percent <- Variants_Exp52 %>%
  group_by(genome_pos, gene, protein_change) %>%
  dplyr::summarise(count = n()) %>%
  mutate(percentage = ((count / 48) * 100),
    species="human")
Variants_Exp52_percent <- Variants_Exp52_percent %>%
filter(!grepl("fs", protein_change))

length(unique(Variants_Exp51$sample_id)) #73
#Calculate percentages
Variants_Exp51_percent <- Variants_Exp51 %>%
  group_by(genome_pos, gene, protein_change) %>%
  dplyr::summarise(count = n()) %>%
  mutate(percentage = ((count / 73) * 100),
  species="mosquito")
Variants_Exp51_percent <- Variants_Exp51_percent %>%
filter(!grepl("fs", protein_change))

both<-rbind(Variants_Exp51_percent,Variants_Exp52_percent)

ggplot(both, aes(x = as.factor(protein_change), y = percentage, fill = as.factor(species))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(x = "Protein change", y = "Percentage of Total", fill = "Gene") +
  facet_wrap(~gene,scales="free_x")+
  theme_classic()
ggsave("Overall_frequencies.pdf",width=15, height=10)


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

ggplot(output_df, aes(x = as.factor(protein_change), y = freq_diff*100,fill=as.factor(gene))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(x = "Genome Position", y = "Difference in frequency (human-mosquito)", fill = "Gene") +
  scale_y_continuous(limits=c(-100,100))+
  facet_grid(mosqid~uniqueid,scales="free_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




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


 
