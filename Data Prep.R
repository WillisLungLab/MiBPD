library(tidyverse)
library(phyloseq)
library(plyr)

fung <- read_csv("Human_ITS_ASVs_Decontaminated.csv")
fung.taxa <- read_csv("Human_ITS_Taxa2.csv")
roughmetadata_fung <- read_csv("Human_Metadata2.csv")

#This removes all samples with no BPD severity listed. These samples have been removed in the public version of the metadata, uncomment it if using the complete metadata
#roughmetadata_fung <- roughmetadata_fung[!(is.na(roughmetadata_fung$BPD_Severity)),]

bact <- read_csv("Human_16S_ASVs_Decontaminated.csv")
bact.taxa <- read_csv("Human_16s_Taxa2.csv")
roughmetadata_bact <- read_csv("Human_Metadata2.csv.csv")

#These samples have been removed in the public version of the metadata, uncomment it if using the complete metadata
#roughmetadata_bact <- roughmetadata_bact[!(is.na(roughmetadata_bact$BPD_Severity)),]

####ITS Data Prep####

fung <- fung[rowSums(is.na(fung)) == 0, ]
fung <- fung %>%
  tibble::column_to_rownames("Name")

#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")
fung.taxa <- as.matrix(fung.taxa)

fung.sort <- fung[,order(colnames(fung))]
md.fung <- roughmetadata_fung[which(unlist(roughmetadata_fung$Name) %in% colnames(fung.sort)),]
#Summarizes read counts of samples, useful sanity check
summary(colSums(fung.sort))

md.fung <- md.fung %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <100, 5 removed
fung.sort2 <- fung.sort[,c(which(colSums(fung.sort)>=100))]
md.fung.sort2 <- md.fung[which(row.names(md.fung) %in% colnames(fung.sort2)),]

OTU_rough = otu_table(fung.sort2, taxa_are_rows = TRUE)
TAX_rough = tax_table(fung.taxa)
samples_rough = sample_data(md.fung.sort2)

#Convert ITS data to phyloseq object 
exp2_fung_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

#Sums up number of samples a given ASV appears in (function written by Dr. Tipton)
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

#Apply function and remove ASVs only present in 1 or 0 samples
fung.nsampls2 <- apply(fung.sort2,1,gt0)
fung.clean <- fung.sort2[which(fung.nsampls2>1),]


#Setup to create the phyloseq object
OTU = otu_table(fung.clean, taxa_are_rows = TRUE)
TAX = tax_table(fung.taxa)
samples = sample_data(md.fung.sort2)

#Convert ITS data to phyloseq object 
exp2_fung <- phyloseq(OTU, TAX, samples)

#Aggregate to Genus level
exp2_fung_gen <- tax_glom(exp2_fung, taxrank = "Genus")

# 5% prevalence filter to remove spurious ASVs from ASV and genus-aggregated phyloseq objects
exp2_fung_prev <- filter_taxa(exp2_fung, function(x) sum(x >= 1) > (0.05*length(x)), TRUE)
exp2_fung_gen_prev <- filter_taxa(exp2_fung_gen , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)

#Subset to BPD or non-BPD only
exp2_fung_prev_bpd <- subset_samples(exp2_fung_prev, BPD=="BPD")
exp2_fung_prev_pprd <- subset_samples(exp2_fung_prev, BPD=="No_BPD")



####16S DATA PREP####

#Set sample names as the row names
bact <- bact %>%
  tibble::column_to_rownames("Name")

#Set sample names as row names and convert to a matrix to make converting to phyloseq easier
bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("Name")
bact.taxa <- as.matrix(bact.taxa)

bact.sort <- bact[,order(colnames(bact))]
md.bact <- roughmetadata_bact[which(unlist(roughmetadata_bact$Name) %in% colnames(bact.sort)),]
summary(colSums(bact.sort))


md.bact <- md.bact %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <1000, 0 removed
bact.sort2 <- bact.sort[,c(which(colSums(bact.sort)>=1000))]
md.bact.sort2 <- md.bact[which(row.names(md.bact) %in% colnames(bact.sort2)),]

OTU_rough = otu_table(bact.sort2, taxa_are_rows = TRUE)
TAX_rough = tax_table(bact.taxa)
samples_rough = sample_data(md.bact.sort2)

#Convert ITS data to phyloseq object 
exp2_bact_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

#Apply Dr.Tipton's function to the 16S data
bact.nsampls <- apply(bact.sort2, 1, gt0)
bact.clean <- bact.sort2[which(bact.nsampls>1),]


#Setup to create the phyloseq object
OTU = otu_table(bact.clean, taxa_are_rows = TRUE)
TAX = tax_table(bact.taxa)
samples = sample_data(md.bact.sort2)
exp2_bact <- phyloseq(OTU, TAX, samples)

#Aggregate to Genus level
exp2_bact_gen <- tax_glom(exp2_bact, taxrank = "Genus")

# 5% prevalence filter to remove spurious ASVs
exp2_bact_prev <- filter_taxa(exp2_bact , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)
exp2_bact_gen_prev <- filter_taxa(exp2_bact_gen , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)

#Subset to BPD or non-BPD only
exp2_bact_prev_bpd <- subset_samples(exp2_bact_prev, BPD=="BPD")
exp2_bact_prev_pprd <- subset_samples(exp2_bact_prev, BPD=="No_BPD")


####MULTIKINGDOM DATA PREP####

#Combine 16S and ITS dataframes and fill in empty sections as needed
combined <- rbind.fill(bact,fung)
rownames(combined) <- c(row.names(bact),row.names(fung))
combined[is.na(combined)] <- 0

#Combine 16S and ITS taxonomy tables
combined.taxa <- rbind.fill(as.data.frame(bact.taxa),as.data.frame(fung.taxa))
rownames(combined.taxa) <- c(row.names(bact.taxa),row.names(fung.taxa))

roughmetadata_combined <- read_csv("Human_Metadata2.csv")

#Uncomment if using the complete metadata instead of the public version
#roughmetadata_combined <- roughmetadata_combined[!(is.na(roughmetadata_combined$BPD_Severity)),]

combined.taxa <- as.matrix(combined.taxa)

#Load metadata and the rest of this is all the same
combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined[which(unlist(roughmetadata_combined$Name) %in% colnames(combined.sort)),]
summary(colSums(combined.sort))

md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <100, 0 removed
combined.sort2 <- combined.sort[,c(which(colSums(combined.sort)>=100))]
md.combined.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.sort2)),]

#This is the same
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

combined.nsampls <- apply(combined.sort2, 1, gt0)
combined.clean <- combined.sort2[which(combined.nsampls>1),]

#Setup to create the phyloseq object
OTU = otu_table(combined.clean, taxa_are_rows = TRUE)
TAX = tax_table(combined.taxa)
samples = sample_data(md.combined.sort2)

#Convert!
exp2_combined <- phyloseq(OTU, TAX, samples)
exp2_combined_prev <- filter_taxa(exp2_combined , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)

