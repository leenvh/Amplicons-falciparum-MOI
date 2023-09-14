# Load libraries
library("HaplotypR")
library("ShortRead")

# Define output directory
outputDir <- "/mnt/storage9/leen/MOI_ampliconseq/14_09_23_Pfs47_0.005"
# Create output directory
if(!dir.exists(outputDir))
  dir.create(outputDir, recursive=T)

# Specify number of pools to analyse
Poolnum <- 14

# Specify mismatch options
minMMrate <- 0.5
minOccGen <- 2

# Specify haplotype options
minCov <- 3
detectionLimit <- 1/200
minOccHap <- 2
minCovSample <- 25

# source edited funs
source("/mnt/storage9/leen/MOI_ampliconseq/scripts/deplexBySample_ext.R")
source("/mnt/storage9/leen/MOI_ampliconseq/scripts/checkBarcode_ext.R")
source("/mnt/storage9/leen/MOI_ampliconseq/scripts/createHaplotypeTable_ext.R")

##############Run demultiplexing by sample 
# Define a function to perform demultiplexing
demultiplex <- function(pattern, outputDir, sampleNum) {
  # Set input file path
  primerFile <- paste0("/mnt/storage9/leen/MOI_ampliconseq/primerfiles/markerFile_Pfs47.txt")
  sampleFile <- paste0("/mnt/storage9/leen/MOI_ampliconseq/samplefiles/Pfs47/", "sampleFile_", pattern, ".txt")
  fnBarcodeF <- paste0("/mnt/storage9/leen/MOI_ampliconseq/barcodefiles/barcode_Fwd.fasta")
  fnBarcodeR <- paste0("/mnt/storage9/leen/MOI_ampliconseq/barcodefiles/barcode_Rev.fasta")
  
  # Set reads
  reads <- list.files("/mnt/storage9/leen/MOI_ampliconseq/fastq/", pattern = paste0(pattern, "_R"), full.names = TRUE)
  
  # Create output subdirectory
  outDeplexSample <- file.path(outputDir, paste0("dePlexSample", sampleNum))
  dir.create(outDeplexSample, showWarnings = FALSE)
  
  # Demultiplex by samples - no SNPs or deletions allowed
  dePlexSample <- deplexSampleExtended(reads[1], reads[2], fnBarcodeF, fnBarcodeR, outDeplexSample, max.mismatch = 0, with.indels = FALSE)
  
  sampleTab <- read.delim(sampleFile, stringsAsFactors = FALSE)
  dePlexSample <- renameDemultiplexedFiles(sampleTab, dePlexSample)
  
  write.table(dePlexSample, file.path(outputDir, paste0("demultiplex", sampleNum, "SampleSummary.txt")), sep = "\t", row.names = FALSE)
  dePlexSample <- na.omit(dePlexSample)
  # Save the dePlexSample dataframe to the R environment with a variable name based on sampleNum
  assign(paste0("dePlexSample", sampleNum), dePlexSample, envir = .GlobalEnv)
  # Save the primerFile to the R environment
  assign("primerFile", primerFile, envir = .GlobalEnv)
}


# Set the folder path to fastq files
folder_path <- "/mnt/storage9/leen/MOI_ampliconseq/fastq/"

# List all the files in the folder
file_list <- list.files(folder_path, full.names = TRUE)
file_list <- file_list[grep("AUG", file_list)]

# Define a function to extract the sample name from the file name
extract_sample_name <- function(file_name) {
  # Extract the sample name from the file name (e.g., Pool1Rep1Oct22_Combined_R1)
  # Modify this logic as needed based on your file naming pattern
  sample_name <- gsub("_R[12].*\\.fastq\\.gz", "", file_name)
  return(sample_name)
}

generate_and_execute_demultiplex_lines <- function(file_list, outputDir) {
  unique_names <- character(0)
  
  # Iterate over the file list and generate demultiplex lines
  for (i in 1:length(file_list)) {
    file_name <- basename(file_list[i])
    sample_name <- extract_sample_name(file_name)
    unique_name <- sample_name
    
    # Check if the unique_name has already been encountered
    if (unique_name %in% unique_names) {
      # If it's a duplicate, skip generating the line
      next  # Skip to the next iteration
    }
    
    # Add the unique_name to the set of encountered names
    unique_names <- union(unique_names, unique_name)
    
    demultiplex_line <- paste0("demultiplex(\"", unique_name, "\", outputDir, ", i, ")")
    cat(demultiplex_line, "\n")
    # Execute the generated line as R code
    eval(parse(text = demultiplex_line))
  }
}

# Call the function with file_list and outputDir
generate_and_execute_demultiplex_lines(file_list, outputDir)


###################merge dePlexSample dataframes
# Create an empty list to store the dataframes
dePlexSampleList <- list()

# Loop to load dataframes and add them to the list
for (i in seq(1, Poolnum*2, by = 2)) {
  # Assuming your dataframes are named dePlexSample1, dePlexSample2, etc.
  df_name <- paste0("dePlexSample", i)
  
  # Check if the dataframe exists in your environment
  if (exists(df_name)) {
    dePlexSampleList[[i]] <- get(df_name)
  }
}

# Merge dataframes using do.call and rbind
if (length(dePlexSampleList) > 0) {
  dePlexSample <- do.call(rbind, dePlexSampleList)
} else {
  # Handle the case where no dataframes were found
  cat("No dataframes found.")
}

write.table(dePlexSample, file.path(outputDir, "demultiplexSampleSummary.txt"), sep="\t", row.names=F)

################Run demultiplex by marker and truncate primer sequence
# create output subdirectory
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)

# process each marker
markerTab <- read.delim(primerFile, stringsAsFactors=F)
dePlexMarker <- demultiplexByMarker(dePlexSample, markerTab, outDeplexMarker)

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
#remove lines without sampleID
dePlexMarker <- na.omit(dePlexMarker)

#############Second method work for non-overlapping sequence read pairs 
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

# Trim options
numNtF <- 221
numNtR <- 221
postfix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)

# Adjust reference to trim options and save as fasta file
refSeq <- as.character(markerTab$ReferenceSequence)
refSeq <- DNAStringSet(paste(substr(refSeq, 1,numNtF), substr(refSeq, nchar(refSeq)+1-numNtR, nchar(refSeq)), sep=""))
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

# Fuse paired read
procReads <- bindAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), outProcFiles, 
                         read1Length=numNtF, read2Length=numNtR)
procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F)


###############Calculate mismatch rate and call SNPs
# process each marker
snpLst <- lapply(markerTab$MarkerID, function(marker){
  # Calculate mismatch rate
  seqErrLst <- calculateMismatchFrequencies(as.character(procReads[procReads$MarkerID == marker, "ReadFile"]),
                                            refSeq[marker],
                                            method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"),
                                            minCoverage=100L)
  names(seqErrLst) <- procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)

  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  snps <- data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N", stringsAsFactors=F)
  write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt",
                                                      minMMrate*100, minOccGen, marker, postfix)),
              row.names=F, col.names=T, sep="\t", quote=F)

  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png",
                                   minMMrate*100, minOccGen, marker, postfix)),
      width=1500 , height=600)
  matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
  abline(v=snps[,"Pos"], lty=2, col="grey")
  abline(h=minMMrate, lty=1, col="red")
  dev.off()

  return(snps)
})
names(snpLst) <- markerTab$MarkerID

##################call haplotypes
# remove samples without reads
procReads <- procReads[procReads$numRead>0,]

# call final haplotypes
finalTab <- createFinalHaplotypTableExtended(
  outputDir = outputDir, sampleTable = procReads, markerTable = markerTab, referenceSeq = refSeq,
  snpList = snpLst, postfix = postfix, minHaplotypCoverage = minCov, minReplicate = minOccHap,
  detectability = detectionLimit, minSampleCoverage = minCovSample)

pfs47_1_finalTab<-finalTab[["pfs47_1"]]
pfs47_2_finalTab<-finalTab[["pfs47_2"]]

write.table(pfs47_1_finalTab, file.path(outputDir, "pfs47_1_finalTab.txt"), sep="\t", row.names=F, quote=F)
write.table(pfs47_2_finalTab, file.path(outputDir, "pfs47_2_finalTab.txt"), sep="\t", row.names=F, quote=F)
