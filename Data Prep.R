library(tidyverse)
library(phyloseq)

####ITS Data Prep####

#Load ASV table
fung <- fung %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")
fung.taxa <- as.matrix(fung.taxa)

#Load metadata and the rest of this is all the same
fung.sort <- fung[,order(colnames(fung))]
md.fung <- roughmetadata_fung[which(unlist(roughmetadata_fung$Name) %in% colnames(fung.sort)),]
#Summarizes read counts of samples, useful sanity check
summary(colSums(fung.sort[,-1]))

md.fung <- md.fung %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <250, 21 removed
fung.sort2 <- fung.sort[,-c(1,which(colSums(fung.sort[,-1])<250))]
md.fung.sort2 <- md.fung[which(row.names(md.fung) %in% colnames(fung.sort2)),]

#Sums up number of samples a given ASV appears in (function written by Dr. Tipton)
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

#Apply function and remove ASVs only present in 1 or 0 samples
fung.nsampls <- apply(fung.sort2, 1, gt0)
fung.clean <- fung.sort2[which(fung.nsampls>1),]


#Setup to create the phyloseq object
OTU = otu_table(fung.clean, taxa_are_rows = TRUE)
TAX = tax_table(fung.taxa)
samples = sample_data(md.fung.sort2)

#Convert ITS data to phyloseq object 
exp2_fung <- phyloseq(OTU, TAX, samples)

#Aggregate to Genus level
exp2_fung_gen <- tax_glom(exp2_fung, taxrank = "Genus")

# 10% prevalence filter to remove spurious ASVs from ASV and genus-aggregated phyloseq objects
exp2_fung_prev <- filter_taxa(exp2_fung, function(x) sum(x >= 1) > (0.10*length(x)), TRUE)
exp2_fung_gen_prev <- filter_taxa(exp2_fung_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)

#Subset to BPD or non-BPD only
exp2_fung_gen_prev_bpd <- subset_samples(exp2_fung_gen_prev, BPD=="BPD")
exp2_fung_gen_prev_pprd <- subset_samples(exp2_fung_gen_prev, BPD=="PPRD")


####MULTIKINGDOM DATA PREP####

#Load ASV table
#Set sample names as the row names
combined <- combined %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Set sample names as row names and convert to a matrix to make converting to phyloseq easier
combined.taxa <- combined.taxa %>%
  tibble::column_to_rownames("Name")
combined.taxa <- as.matrix(combined.taxa)

#Load metadata and the rest of this is all the same
combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined[which(unlist(roughmetadata_combined$Name) %in% colnames(combined.sort)),]
summary(colSums(combined.sort[,-1]))


md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <250, 2 removed
combined.sort2 <- combined.sort[,-c(1,which(colSums(combined.sort[,-1])<250))]
md.combined.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.sort2)),]
combined.nsampls <- apply(combined.sort2, 1, gt0)
combined.clean <- combined.sort2[which(combined.nsampls>1),]


#Setup to create the phyloseq object
OTU = otu_table(combined.clean, taxa_are_rows = TRUE)
TAX = tax_table(combined.taxa)
samples = sample_data(md.combined.sort2)
exp2_combined <- phyloseq(OTU, TAX, samples)

#Aggregate to Genus level
exp2_combined_gen <- tax_glom(exp2_combined, taxrank = "Genus")

# 10% prevalence filter to remove spurious ASVs
exp2_combined_prev <- filter_taxa(exp2_combined , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)
exp2_combined_gen_prev <- filter_taxa(exp2_combined_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)

#Subset to BPD or non-BPD only
exp2_combined_gen_prev_bpd <- subset_samples(exp2_combined_gen_prev, BPD=="BPD")
exp2_combined_gen_prev_pprd <- subset_samples(exp2_combined_gen_prev, BPD=="PPRD")
