library(dada2)
library(tidyverse)
library(phyloseq)
set.seed(1312)

####################### ITS ASV TABLE #######################

#Load in fastq files and plot quality profile as a sanity check
fnFs <- sort(list.files(path = "path/to/your/files", pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path = "path/to/your/files", pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])

filtFs <- file.path(path = "path/to/your/files", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path = "path/to/your/files", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter reads and truncate low quality reads (in this case, not setting a truncation length due to the variable length of the ITS region)
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, matchIDs = FALSE, multithread = FALSE)

#Fit error model, then plot as a sanity check
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Denoise reads, pooling for greater recovery of rare reads
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = TRUE)

#Merge three different ways: Merge with no mismatches, merge with some mismatches allowed, and concatenation
mergedData <- mergePairs(dadaFs,filtFs, dadaRs, filtRs,maxMismatch=0, trimOverhang=TRUE, returnRejects=TRUE, verbose=TRUE)
mismatch_mergers <- mergePairs(dadaFs,filtFs, dadaRs, filtRs, maxMismatch=2, trimOverhang=TRUE, returnRejects=TRUE, verbose=TRUE)
concats <- mergePairs(dadaFs,filtFs, dadaRs, filtRs, justConcatenate=TRUE, verbose=TRUE)

#This function keeps the zero mismatch merging where available, then adds either mismatched merging (where available) or concatenated reads to rows that did not merge under the no mismatch parameters
for(i in names(mergedData)) {
  # store row index to drop certain ASVs later
  rowsToDelete = vector()
  
  mergedDf = mergedData[[i]]
  cat(i, "Out of total", sum(mergedDf$abundance), "paired-reads (in", nrow(mergedDf), "unique pairings), retained ")
  mismatchDf = mismatch_mergers[[i]]
  rownames(mismatchDf) = paste(mismatchDf$forward, mismatchDf$reverse, sep="_")
  concatDf = concats[[i]]
  rownames(concatDf) = paste(concatDf$forward, concatDf$reverse, sep="_")
  
  for (row in 1:nrow(mergedDf)) {
    if (mergedDf[row,]$accept) { next }
    
    uniquePairID = paste(mergedDf[row,]$forward, mergedDf[row,]$reverse, sep="_")
    
    if (mismatchDf[uniquePairID,]$nmatch <= 12) {
      mergedDf[row,] = concatDf[uniquePairID,]
    }
    else{
      misMatchIndels = mismatchDf[row,]$nmismatch + mismatchDf[row,]$nindel
      cutOff = 0
      if ( mismatchDf[uniquePairID,]$nmatch > 12 && mismatchDf[uniquePairID,]$nmatch <= 50 ) {
        cutOff = 1
      }
      else if ( mismatchDf[uniquePairID,]$nmatch > 50 && mismatchDf[uniquePairID,]$nmatch <= 100 ) {
        cutOff = 2
      }
      else if (mismatchDf[uniquePairID,]$nmatch > 100) {
        cutOff = 3
      }
      # check if mismatches are below cut off
      if (misMatchIndels <= cutOff) {
        mergedDf[row,] = mismatchDf[uniquePairID,]
      } 

      else {
        trimLength = mismatchDf[uniquePairID,]$nmatch
        concatDf[uniquePairID,]$sequence = gsub(paste0("N{10}[A-z]{", trimLength, "}"), "NNNNNNNNNN", concatDf[uniquePairID,]$sequence)
        mergedDf[row,] = concatDf[uniquePairID,]
      }
    }
  }
  
  if (length(rowsToDelete) == 0){
    mergedData[[i]] = mergedDf
  } else {
    mergedData[[i]] = mergedDf[-rowsToDelete,]
  }
  cat(sum(mergedData[[i]]$abundance), "paired-reads (in", nrow(mergedData[[i]]), "unique pairings)\n")
}

seqtab <- makeSequenceTable(mergedData)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "Gut ITS Sequence Table No Chimera.rds")

#Tracks number of reads through each step as a sanity check
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergedData, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy using Unite 9.0 database
set.seed(1312)
taxa <- assignTaxonomy(seqtab.nochim, "path/to/database.fasta", tryRC = TRUE, multithread=TRUE, verbose = TRUE)

seqtab.nochim <- as.data.frame(seqtab.nochim)
colnames(seqtab.nochim) <- paste("ASV", 1:4974, sep="")

taxa_its <- as.data.frame(taxa)
repseqs_its <- as.matrix(row.names(taxa_its))
row.names(repseqs_its) <- paste("ASV", 1:4974, sep="")
row.names(taxa_its) <- colnames(seqtab.nochim)

#Save ASV table
#Save taxonomy table
#Save representative sequences

#Remove Chytridiomycota or taxa undescribed past the Phylum level (sequences found to be host contamination by BLAST), also remove sequences with NA at every taxonomic level
fung_taxa <- taxa_its[taxa_its$Phylum != "p__Chytridiomycota" & taxa_its$Phylum != "p__Fungi_phy_Incertae_sedis",]
fung_taxa <- fung_taxa[!apply(is.na(fung_taxa) | fung_taxa == "", 1, all),]

#Remove everything from the OTU table that isn't in the cleaned up taxa table
its <- seqtab.nochim[,colnames(seqtab.nochim) %in% row.names(fung_taxa)]
its <- t(as_tibble(its, rownames = NA))
its.taxa <- as.matrix(fung_taxa)

#Save filtered ASV table
#Save filtered taxonomy table


####################### 16S ASV TABLE #######################

#Load in files and make a directory for filtered files
fnFs2 <- sort(list.files(path = "path/to/your/files", pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs2 <- sort(list.files(path = "path/to/your/files", pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)
plotQualityProfile(fnFs2[1:2])

filtFs2 <- file.path(path = "path/to/your/files", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs2 <- file.path(path = "path/to/your/files", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs2) <- sample.names
names(filtRs2) <- sample.names

#Filter and truncate based on length (ensure you still have the necessary length for merging, and check quality profile to find a good truncation length)
out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2,truncLen=c(220,220),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=FALSE)

#Fit error model
errF2 <- learnErrors(filtFs2, multithread=TRUE)
errR2 <- learnErrors(filtRs2, multithread=TRUE)
plotErrors(errF2, nominalQ=TRUE)

#Dereplicate, pseudo-pooling to make computation easier and since rare ASVs aren't as necessary as ITS
dadaFs2 <- dada(filtFs2, err=errF2, multithread=TRUE, pool = "pseudo")
saveRDS(dadaFs2, "dadaFs2_16s_gut.rds")
dadaRs2 <- dada(filtRs2, err=errR2, multithread=TRUE, pool = "pseudo")
saveRDS(dadaRs2, "dadaRs2_16s_gut.rds")

#Merge using same strategy as in ITS
mergedData2 <- mergePairs(dadaFs2,filtFs2, dadaRs2, filtRs2,maxMismatch=0, trimOverhang=TRUE, returnRejects=TRUE, verbose=TRUE)
mismatch_mergers2 <- mergePairs(dadaFs2,filtFs2, dadaRs2, filtRs2, maxMismatch=2, trimOverhang=TRUE, returnRejects=TRUE, verbose=TRUE)
concats2 <- mergePairs(dadaFs2,filtFs2, dadaRs2, filtRs2, justConcatenate=TRUE, verbose=TRUE)

# replace mismatched or concatenated ASVs in the main mergedData2
for(i in names(mergedData2)) {
  # store row index to drop certain ASVs later
  rowsToDelete = vector()
  
  mergedDf = mergedData2[[i]]
  cat(i, "Out of total", sum(mergedDf$abundance), "paired-reads (in", nrow(mergedDf), "unique pairings), retained ")
  mismatchDf = mismatch_mergers2[[i]]
  rownames(mismatchDf) = paste(mismatchDf$forward, mismatchDf$reverse, sep="_")
  concatDf = concats2[[i]]
  rownames(concatDf) = paste(concatDf$forward, concatDf$reverse, sep="_")
  
  for (row in 1:nrow(mergedDf)) {
    # skipping rows that are good to go from default analysis
    if (mergedDf[row,]$accept) { next }
    
    uniquePairID = paste(mergedDf[row,]$forward, mergedDf[row,]$reverse, sep="_")
    
    if (mismatchDf[uniquePairID,]$nmatch <= 12) {
      mergedDf[row,] = concatDf[uniquePairID,]
    }
    else{
      misMatchIndels = mismatchDf[row,]$nmismatch + mismatchDf[row,]$nindel
      cutOff = 0
      if ( mismatchDf[uniquePairID,]$nmatch > 12 && mismatchDf[uniquePairID,]$nmatch <= 50 ) {
        cutOff = 1
      }
      else if ( mismatchDf[uniquePairID,]$nmatch > 50 && mismatchDf[uniquePairID,]$nmatch <= 100 ) {
        cutOff = 2
      }
      else if (mismatchDf[uniquePairID,]$nmatch > 100) {
        cutOff = 3
      }
      if (misMatchIndels <= cutOff) {
        mergedDf[row,] = mismatchDf[uniquePairID,]
      } 
      else {
        trimLength = mismatchDf[uniquePairID,]$nmatch
        concatDf[uniquePairID,]$sequence = gsub(paste0("N{10}[A-z]{", trimLength, "}"), "NNNNNNNNNN", concatDf[uniquePairID,]$sequence)
        mergedDf[row,] = concatDf[uniquePairID,]
      }
    }
  }
  
  if (length(rowsToDelete) == 0){
    mergedData2[[i]] = mergedDf
  } else {
    mergedData2[[i]] = mergedDf[-rowsToDelete,]
  }
  cat(sum(mergedData2[[i]]$abundance), "paired-reads (in", nrow(mergedData2[[i]]), "unique pairings)\n")
}

seqtab2 <- makeSequenceTable(mergedData2)
table(nchar(getSequences(seqtab2)))
#Save sequence table
saveRDS(seqtab2, "16s_gut_seqtab_all.rds")

#Remove chimeras
seqtab.nochim2 <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim2)

#Track read counts through pipeline
getN <- function(x) sum(getUniques(x))
track2 <- cbind(out2, sapply(dadaFs2, getN), sapply(dadaRs2, getN), sapply(mergedData2, getN), rowSums(seqtab.nochim2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs2, getN) with getN(dadaFs2)
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track2) <- sample.names
head(track2)


#Save sequence table as csv (and as R dataset to back up your progress before taxonomy assignment, it is memory intensive and may crash R)
saveRDS(seqtab.nochim2, "16s_gut_seqtab_nochim_all.rds")
set.seed(1312)

#Assign taxonomy using SILVA v138.1 Database
taxa2 <- assignTaxonomy(seqtab.nochim2, "path/to/database.fa.gz", tryRC = TRUE, multithread=TRUE, verbose = TRUE)

seqtab.nochim2 <- as.data.frame(seqtab.nochim2)
colnames(seqtab.nochim2) <- paste("ASV", 1:18047, sep="")

taxa_16s <- as.data.frame(taxa2)
repseqs_16s <- as.matrix(row.names(taxa_16s))
row.names(repseqs_16s) <- paste("ASV", 1:18047, sep="")
row.names(taxa_16s) <- colnames(seqtab.nochim2)

#Save ASV table
#Save taxonomy table
#Save repseqs

#Remove Mitochondria, Cyanobacteria (usually also mitochondrial contamination), and Eukaryotes (usually host contamination), as well as rows with NA at all taxonomic levels
bact_taxa <- taxa_16s[taxa_16s$Phylum != "Cyanobacteria" & taxa_16s$Kingdom != "Eukaryota" & taxa_16s$Family != "Mitochondria",]
bact_taxa <- bact_taxa[!apply(is.na(bact_taxa) | bact_taxa == "", 1, all),]

#Remove everything from the ASV table that isn't in the cleaned up taxa table
bact <- seqtab.nochim2[,colnames(seqtab.nochim2) %in% row.names(bact_taxa)]
bact <- t(as_tibble(bact), rownames = NA)
bact.taxa <- as.matrix(bact_taxa)

#Save filtered ASV table
#Save filtered taxonomy table
