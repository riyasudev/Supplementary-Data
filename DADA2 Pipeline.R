if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
library(dada2); packageVersion("dada2")
getwd()

#I created a folder with FastQ files for all samples (_1.fastq and _2.fastq)

path <- "C:/Users/riyas/Documents/Biomed Year 3/In Silico Research Project/SRA Reads/FinalReadsASD"
list.files(path)

#read in the names of the fastq files, and perform some string manipulation to get matched lists of the forward and reverse fastq files. #careful this will read all files in the folder you specified… (so remove old ones or ones you arent interested in)
#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
# Set file paths
fnFs <- sort(list.files(path, pattern = "_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_002.fastq", full.names = TRUE))

# Check they were loaded correctly: tells you how many file of each type have been loaded
length(fnFs)  
length(fnRs)  

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#check/plot quality of F reads - below will show first 2 samples: #ok to skip, unless want to check
plotQualityProfile(fnFs[1:4])
#check/plot quality of F reads:
plotQualityProfile(fnRs[1:4])

#Filter and trim:
#The tutorial is using 2x250 V4 sequence data
#Assign the filenames for the filtered fastq.gz files that will be created:
#AND place filtered files in filtered/ subdirectoryc
# 1. Create 'filtered' directory if it doesn't exist
dir.create(file.path(path, "filtered"))

# 2. Define the paths for filtered files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#check for duplicate file names – if True is produced, check input file names! (see above about hyphens!)
any(duplicated(c(fnFs, fnRs)))
any(duplicated(c(filtFs, filtRs)))

list.files(path)  # Check if the 'filtered' folder is listed

#use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed
#in a read, which is a better filter than simply averaging quality scores.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,250),trimLeft=c(20,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)

#On Windows set multithread=FALSE

#above default for truncLen was 270,250 for truncLen
#NOTE: If you want to speed up downstream computation, consider tightening maxEE. If too
#few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse
#reads (eg. maxEE=c(2,5))/default is 2:2, and reducing the truncLen to remove low quality tails.
#Remember though, when choosing truncLen for paired-end reads you must maintain overlap
#after truncation in order to merge them later.

#Learn the Error Rates:
#calculate the error rates for F reads (#can take some time):
errF <- learnErrors(filtFs, multithread=TRUE)
#visualise the error rates:
plotErrors(errF, nominalQ=TRUE)
#calculate the error rates for R reads (#can take some time):
errR <- learnErrors(filtRs, multithread=TRUE)
#visualise the error rates:
plotErrors(errR, nominalQ=TRUE)


#Sample Inference (making ASVs):
#We are now ready to apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


#Inspecting the returned dada-class object:
dadaFs[[1]]

#+Extensions: By default, the dada function processes each sample independently.
#However, pooling information across samples can increase sensitivity to sequence variants that may be
#present at very low frequencies in multiple samples. The dada2 package offers two types of
#pooling. dada(..., pool=TRUE) performs standard pooled processing, in which all samples are
#pooled together for sample inference. dada(..., pool="pseudo") performs pseudo-pooling, in
#which samples are processed independently after sharing information between samples,
#approximating pooled sample inference in linear time.


#Merge paired reads:
#We now merge the forward and reverse reads together to obtain the full denoised sequences.
#Merging is performed by aligning the denoised forward reads with the reverse-complement of
#the corresponding denoised reverse reads, and then constructing the merged “contig”
#sequences. By default, merged sequences are only output if the forward and reverse reads
#overlap by at least 12 bases, and are identical to each other in the overlap region (but these
#conditions can be changed via function arguments).


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

head(mergers[[1]])


#prints unique, merged sequences:
uniquesToFasta(getUniques(mergers), fout="C:/FastQ files _raw_clean/uniqueSeqs.fasta",
ids=paste0("Seq", seq(length(getUniques(mergers)))))


#Construct sequence table
#We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
#below will make a csv file for excel of the table (doesnt give ASV numbers, but sequences)
write.csv(seqtab, file = "my_data_dada2_FastQ files _raw_clean.csv")


#Remove chimeras
#The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the
#accuracy of sequence variants after denoising makes identifying chimeric ASVs simpler than when
#dealing with fuzzy OTUs. Chimeric sequences are identified if they can be exactly reconstructed by
#combining a left-segment and a right-segment from two more abundant “parent” sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE,
                                    verbose=TRUE)
dim(seqtab.nochim)


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "read-count-tracking.csv")


#write ASV table
write.csv(seqtab.nochim, file = "dada2_FastQ files _raw_clean.csv")


#make an ASV multi-fasta and give seq headers more manageable names (>ASV_1, >ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_dada2_FastQ files _raw_clean.fasta")

######## Assign Taxonomy – add this to the end of your DADA2 script#################


taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/riyas/Documents/Biomed Year 3/In Silico Research Project/silva_nr99_v138.2_toGenus_trainset.fa", multithread=TRUE)


#Just have a look if it seems correct

taxa.print <- taxa # Removing sequence rownames for display only

rownames(taxa.print) <- NULL

head(taxa.print)


#Save the table in your directory

write.csv(taxa, file = "taxonomy_table.csv")


####################################################
# Read in the taxonomy table you saved
taxa <- read.csv("taxonomy_table.csv", row.names = 1)

# Convert it back to a matrix (assignTaxonomy creates a matrix, not a dataframe)
taxa <- as.matrix(taxa)

# Check it loaded correctly
head(taxa)

# Now run the backfill function
backfill_taxonomy <- function(tax_table) {
  tax_filled <- tax_table
  
  for (i in 1:nrow(tax_filled)) {
    last_known <- NULL
    
    for (j in 1:ncol(tax_filled)) {
      if (is.na(tax_filled[i, j])) {
        if (!is.null(last_known)) {
          rank_name <- colnames(tax_filled)[j]
          tax_filled[i, j] <- paste0(last_known, "_", rank_name)
        }
      } else {
        last_known <- tax_filled[i, j]
      }
    }
  }
  
  return(tax_filled)
}

# Apply backfill
taxa_backfilled <- backfill_taxonomy(taxa)

# Check results
head(taxa_backfilled)
cat("NAs before:", sum(is.na(taxa)), "\n")
cat("NAs after:", sum(is.na(taxa_backfilled)), "\n")

# Save
write.csv(taxa_backfilled, file = "taxonomy_table_backfilled.csv", row.names = TRUE)

# Backfill NA values in taxonomy table
# This propagates the last known taxonomic assignment down to lower ranks

# Function to backfill NAs
backfill_taxonomy <- function(tax_table) {
  # Make a copy to avoid modifying original
  tax_filled <- tax_table
  
  # For each row (ASV)
  for (i in 1:nrow(tax_filled)) {
    # Start from Kingdom and move through ranks
    last_known <- NULL
    
    for (j in 1:ncol(tax_filled)) {
      # If current level is NA
      if (is.na(tax_filled[i, j])) {
        # Fill with last known + rank name
        if (!is.null(last_known)) {
          rank_name <- colnames(tax_filled)[j]
          tax_filled[i, j] <- paste0(last_known, "_", rank_name)
        }
      } else {
        # Update last known taxonomy
        last_known <- tax_filled[i, j]
      }
    }
  }
  
  return(tax_filled)
}

# Apply backfill to your taxonomy table
taxa_backfilled <- backfill_taxonomy(taxa)

# Check the results
head(taxa_backfilled)

# Compare before and after
cat("NAs before backfill:", sum(is.na(taxa)), "\n")
cat("NAs after backfill:", sum(is.na(taxa_backfilled)), "\n")

# Save the backfilled taxonomy table
write.csv(taxa_backfilled, file = "taxonomy_table_backfilled.csv", row.names = TRUE)











