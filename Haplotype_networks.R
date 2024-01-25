#Haplotype networks

library("ape")
library("pegas")
library("RColorBrewer")
library("dplyr")
library("Biostrings")
library(lubridate)
setwd("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork")

#Construct FASTA files to plot
TRAP <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/Matching_Haplotypes_TRAP_Molten.csv",header=TRUE,sep=',')
CSP <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Haplotypr_Analysis/Matching_Haplotypes_CSP_Molten.csv",header=TRUE,sep=',')
TRAP_haplotypes <- readDNAStringSet("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork/trap_HaplotypeSeq_merge.fasta", format="fasta")
CSP_haplotypes <- readDNAStringSet("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork/csp_HaplotypeSeq_merge.fasta", format="fasta")
ClinicalData <- read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/Paper_AmpSeq/ClinicalData.csv",header=TRUE,sep=',')

ClinicalData$Timepoint <- paste0(ClinicalData$studycode,"_Day",ClinicalData$studyvisit)
TRAP <- TRAP %>%
  left_join(ClinicalData %>% select(Timepoint, visitdate_hb), by = "Timepoint") %>%
  mutate(visitdate = visitdate_hb) %>%
  select(-visitdate_hb) 
TRAP$visitdate <- as.Date(TRAP$visitdate)
TRAP$visitdate <- month(TRAP$visitdate, label = TRUE, abbr = FALSE)

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

# Function to extract sequence by haplotype and write to file
write_haplotype_sequences <- function(sample_id, haplotype, day, species, sequence, file_path) {
  # Construct the FASTA header
  header <- paste0(">", sample_id, "_", haplotype, "_", day, "_",species)
  # Combine header and sequence
  fasta_entry <- paste(header, sequence, sep="\n")
  # Write to file
  cat(fasta_entry, file=file_path, append=TRUE, sep="\n")
}

# Create a unique key for each SampleID and Haplotype
TRAP <- TRAP %>%
  mutate(HaploKey = paste(SampleID, Haplotype, sep="_"))
CSP <- CSP %>%
  mutate(HaploKey = paste(SampleID, Haplotype, sep="_"))

# Iterate over the dataframe and write sequences to the FASTA file
for (i in 1:nrow(TRAP)) {
  individual <- TRAP$Individual[i]
  haplotype <- TRAP$Haplotype[i]
  species <- TRAP$species[i]
  day <- TRAP$Day[i]
  
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
for (i in 1:nrow(CSP)) {
  individual <- CSP$Individual[i]
  haplotype <- CSP$Haplotype[i]
  species <- CSP$species[i]
  day <- CSP$Day[i]
  
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

seqs <- read.FASTA("TRAP_plotme.fasta")
seqs <- read.FASTA("CSP_plotme.fasta")
haps <- haplotype(seqs)
dist <- dist.dna(haps, "N")
net <- rmst(dist, B=100)
(sz <- summary(haps))
nt.labs <- names(seqs)
#sz <- sz[nt.labs]

# Extract species names from the label names
#species <- sapply(strsplit(nt.labs, "_"), function(x) x[length(x)])
# Create a color vector based on the species
#species_colors <- ifelse(species == "human", "#db7376", "#465c7a")


matrix <- as.matrix(seqs)
rownames(matrix) <- nt.labs
#rownames(matrix)<- sub("_","", rownames(matrix))
#rownames(matrix) <- sub(".*_(.*)$", "\\1", rownames(matrix))
regions <- haploFreq(matrix,split="_",what=3)
write.csv(regions, file = "CSV_regions.csv")
#plot(net, size=sz, scale.ratio =2, cex=0.75, pie = regions, bg=brewer.pal(n=12, name="Set3"), threshold = 0, col.link="black")
#legend("topright", colnames(regions),col=brewer.pal(ncol(regions), name="Set3"), pch=20, cex=0.7)

#col=usecol(pal_unikn_pref)
mycols <- c("#fabfbd", "#fce6cc", "#e0eed3", "#e1d8e9", "#d6edf8","#eaccd8","#fefbde")
mycols <- brewer.pal(7, "BrBG")
plot(net, size=sz, scale.ratio =2, cex=0.75, pie = regions, bg = mycols,threshold = 0, col.link="black")
xy_trap <- list(x=c(71.9171256, 15.3003829, 115.2096901, 113.1091161, 38.4959351,
  53.8231244, 56.9306340, 118.7295029, 94.1127099, 83.2170195,
  -8.4520770, 106.9047401, -15.0606097, 87.3291249, 0.0000000,
  17.7186059, 42.4950721, 0.3525274, 27.8691270, 77.1458774,
  112.1127129, 45.5733324, 175.8792380, 146.4886241, 194.7341999,
  -1.3540001, 82.1096508, 108.2070211, 24.9897178, -7.3361081,
  220.7173830, 8.9436075, 30.6230154, 91.5318280, 164.9220837,
  64.3787494, -3.7313406, 141.5066080, 91.6864469, 13.7775298,
  86.7044307, 32.1715445, 12.5346775, 80.1759477, 92.6828501,
  -24.5655948, 128.6965272, -25.4610692, -32.8016691, 113.7274122,
  134.5317854, 101.3335515, 64.8853865), y=c(
  187.44148841, 151.24629904, 159.05254276, 78.65494766, -17.80759037,
  75.16753639, 34.99541398, 132.85919261, 169.89149155, 80.64775411,
  61.71609289, 200.86517541, 5.88088137, 182.56314339, 0.00000000,
  0.06946311, 1.61799225, 10.86896656, 22.70293509, 8.40339392,
  30.82759300, 223.27504104, 148.40317322, 125.48589911, 178.79347150,
  -26.07582636, 132.45619242, 190.39235447, 69.18911705, 23.38975179,
  141.42835067, 18.13563635, 90.40032930, 206.02383803, 107.55064111,
  85.64638862, -11.62798904, 107.05243950, 16.87794788, -11.62798904,
  99.57941533, 11.94151981, 46.27184294, 156.98708209, 63.70889933,
  26.92973095, 159.56796398, -18.49436650, 2.92830277, 169.37531517,
  57.23227839, 118.00531094, 176.82172192))
replot(xy=xy_trap)
legend("left", colnames(regions),col=mycols, pch=20, cex=0.7)

current_plot <- recordPlot()
pdf("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/HaplotypeNetwork/CSP_network.pdf")
replayPlot(current_plot)
dev.off()

#population genetics
nuc.div(seqs)
tajima.test(seqs)
hap.div(seqs)






