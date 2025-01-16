library(tidyverse)
library(PERFect)
library(phyloseq)
library(plyr)
library(RColorBrewer)


### Fig. S2A: ITS Contaminant Only Taxa Barplot ###

#Read in ASV and taxonomy table of only contaminants (called removed.taxa and fung_removed_taxa in this example, respectively)

fung_contam <- as.data.frame(t(removed.taxa))
fung.contam.taxa <- fung_removed_taxa
roughmetadata_fung_contam <- read_csv("Human_Metadata2.csv")

#This is only necessary with the public version of the metadata, since phyloseq does not allow metadata with only one column
roughmetadata_fung_contam$dummy <- 0

#Remove samples with unknown BPD severity since these samples were excluded from analysis (already removed in the public metadata, uncomment if using the full metadata)
#roughmetadata_fung <- roughmetadata_fung[!(is.na(roughmetadata_fung$BPD_Severity)),]

#### Contaminant ITS Data Prep ####

fung_contam <- fung_contam[rowSums(is.na(fung_contam)) == 0, ]
fung_contam <- fung_contam %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.contam.taxa <- fung.contam.taxa %>%
  tibble::column_to_rownames("Name")

fung.contam.taxa[is.na(fung.contam.taxa)] <- ""

for (i in 1:7){ fung.contam.taxa[,i] <- as.character(fung.contam.taxa[,i])}
####### Fills holes in the tax table
fung.contam.taxa[is.na(fung.contam.taxa)] <- ""
for (i in 1:nrow(fung.contam.taxa)){
  if (fung.contam.taxa[i,2] == ""){
    kingdom <- paste("Unidentified_", fung.contam.taxa[i,1], sep = "")
    fung.contam.taxa[i, 2:7] <- kingdom
  } else if (fung.contam.taxa[i,3] == ""){
    phylum <- paste("Unidentified_", fung.contam.taxa[i,2], sep = "")
    fung.contam.taxa[i, 3:7] <- phylum
  } else if (fung.contam.taxa[i,4] == ""){
    class <- paste("Unidentified_", fung.contam.taxa[i,3], sep = "")
    fung.contam.taxa[i, 4:7] <- class
  } else if (fung.contam.taxa[i,5] == ""){
    order <- paste("Unidentified_", fung.contam.taxa[i,4], sep = "")
    fung.contam.taxa[i, 5:7] <- order
  } else if (fung.contam.taxa[i,6] == ""){
    family <- paste("Unidentified_", fung.contam.taxa[i,5], sep = "")
    fung.contam.taxa[i, 6:7] <- family
  } else if (fung.contam.taxa[i,7] == ""){
    fung.contam.taxa$Species[i] <- paste("Unidentified_",fung.contam.taxa$Genus[i], sep = "")
  }
}

fung_contam <- fung_contam[which(row.names(fung_contam) %in% row.names(fung.contam.taxa)),]
fung.contam.taxa <- as.matrix(fung.contam.taxa)

#Load metadata and the rest of this is all the same
fung_contam.sort <- fung_contam[,order(colnames(fung_contam))]
md.fung.contam <- roughmetadata_fung_contam[which(unlist(roughmetadata_fung_contam$Name) %in% colnames(fung_contam.sort)),]
#Summarizes read counts of samples, useful sanity check
summary(colSums(fung_contam.sort))

md.fung.contam <- md.fung.contam %>%
  tibble::column_to_rownames("Name")


OTU_contam_rough = otu_table(fung_contam.sort, taxa_are_rows = TRUE)
TAX_contam_rough = tax_table(fung.contam.taxa)
samples_contam_rough = sample_data(md.fung.contam)

exp2_fung_contam_rough <- phyloseq(OTU_contam_rough, TAX_contam_rough, samples_contam_rough)
exp2_fung_gen_contam_rough <- tax_glom(exp2_fung_contam_rough, taxrank = "Genus")

#Identify top 49 genera and rename everything else to "Other"
top20_fung_contam_merged_list <- names(sort(taxa_sums(exp2_fung_gen_contam_rough), decreasing=TRUE)[1:49])
top20_merged_fung_contam_rel <- transform_sample_counts(exp2_fung_gen_contam_rough, function(x) x / sum(x) )
top20_merged_fung_contam_df <- psmelt(exp2_fung_gen_contam_rough)
top20_merged_fung_contam_df[!(top20_merged_fung_contam_df$OTU %in% top20_fung_contam_merged_list),]$Genus <- 'Other'


###Barplot of top 50 fungal genera###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(50)

barplot_gen_contam_its <- ggplot(top20_merged_fung_contam_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_contam_its + scale_fill_manual(values = barplot_colors)
#Save



### Fig. S2B: 16S Contaminant Only Taxa Barplot ###

#Read in ASV and taxonomy table of only contaminants (called removed.taxa and bact_removed_taxa in this example, respectively)

bact_contam <- as.data.frame(t(removed.taxa))
bact.contam.taxa <- bact_removed_taxa
roughmetadata_bact_contam <- read_csv("Human_Metadata2.csv")

#This is only necessary with the public version of the metadata, since phyloseq does not allow metadata with only one column
roughmetadata_bact_contam$dummy <- 0

#Remove samples with unknown BPD severity since these samples were excluded from analysis (already removed in the public metadata, uncomment if using the full metadata)
#roughmetadata_bact <- roughmetadata_bact[!(is.na(roughmetadata_bact$BPD_Severity)),]

#### Contaminant 16s Data Prep ####

bact_contam <- bact_contam[rowSums(is.na(bact_contam)) == 0, ]
bact_contam <- bact_contam %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
bact.contam.taxa <- bact.contam.taxa %>%
  tibble::column_to_rownames("Name")

bact.contam.taxa[is.na(bact.contam.taxa)] <- ""

for (i in 1:7){ bact.contam.taxa[,i] <- as.character(bact.contam.taxa[,i])}
####### Fills holes in the tax table
bact.contam.taxa[is.na(bact.contam.taxa)] <- ""
for (i in 1:nrow(bact.contam.taxa)){
  if (bact.contam.taxa[i,2] == ""){
    kingdom <- paste("Unidentified_", bact.contam.taxa[i,1], sep = "")
    bact.contam.taxa[i, 2:7] <- kingdom
  } else if (bact.contam.taxa[i,3] == ""){
    phylum <- paste("Unidentified_", bact.contam.taxa[i,2], sep = "")
    bact.contam.taxa[i, 3:7] <- phylum
  } else if (bact.contam.taxa[i,4] == ""){
    class <- paste("Unidentified_", bact.contam.taxa[i,3], sep = "")
    bact.contam.taxa[i, 4:7] <- class
  } else if (bact.contam.taxa[i,5] == ""){
    order <- paste("Unidentified_", bact.contam.taxa[i,4], sep = "")
    bact.contam.taxa[i, 5:7] <- order
  } else if (bact.contam.taxa[i,6] == ""){
    family <- paste("Unidentified_", bact.contam.taxa[i,5], sep = "")
    bact.contam.taxa[i, 6:7] <- family
  } else if (bact.contam.taxa[i,7] == ""){
    bact.contam.taxa$Species[i] <- paste("Unidentified_",bact.contam.taxa$Genus[i], sep = "")
  }
}

bact_contam <- bact_contam[which(row.names(bact_contam) %in% row.names(bact.contam.taxa)),]
bact.contam.taxa <- as.matrix(bact.contam.taxa)

#Load metadata and the rest of this is all the same
bact_contam.sort <- bact_contam[,order(colnames(bact_contam))]
md.bact.contam <- roughmetadata_bact_contam[which(unlist(roughmetadata_bact_contam$Name) %in% colnames(bact_contam.sort)),]
#Summarizes read counts of samples, useful sanity check
summary(colSums(bact_contam.sort))

md.bact.contam <- md.bact.contam %>%
  tibble::column_to_rownames("Name")

OTU_contam_rough = otu_table(bact_contam.sort, taxa_are_rows = TRUE)
TAX_contam_rough = tax_table(bact.contam.taxa)
samples_contam_rough = sample_data(md.bact.contam)

exp2_bact_contam_rough <- phyloseq(OTU_contam_rough, TAX_contam_rough, samples_contam_rough)
exp2_bact_gen_contam_rough <- tax_glom(exp2_bact_contam_rough, taxrank = "Genus")

#Identify top 49 genera and rename everything else to "Other"
top20_bact_contam_merged_list <- names(sort(taxa_sums(exp2_bact_gen_contam_rough), decreasing=TRUE)[1:49])
top20_merged_bact_contam_rel <- transform_sample_counts(exp2_bact_gen_contam_rough, function(x) x / sum(x) )
top20_merged_bact_contam_df <- psmelt(exp2_bact_gen_contam_rough)
top20_merged_bact_contam_df[!(top20_merged_bact_contam_df$OTU %in% top20_bact_contam_merged_list),]$Genus <- 'Other'


###Barplot of top 50 bacterial genera###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(50)

barplot_gen_contam_16s <- ggplot(top20_merged_bact_contam_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_contam_16s + scale_fill_manual(values = barplot_colors)
#Save



### Fig. S2C: ITS Decontaminated Taxa Barplot ###

#Read in filtered ASV table, table of unremoved taxa, and metadata
fung <- read_csv("Human_ITS_ASVs_Decontaminated.csv")
fung.taxa <- read_csv("Human_ITS_Taxa2.csv")
roughmetadata_fung <- read_csv("Human_Metadata2.csv")

#Remove samples with unknown BPD severity since these samples were excluded from analysis (already removed in the public metadata, uncomment if using the full metadata)
#roughmetadata_fung <- roughmetadata_fung[!(is.na(roughmetadata_fung$BPD_Severity)),]

#### Decontaminated ITS Data Prep ####

fung <- fung[rowSums(is.na(fung)) == 0, ]
fung <- fung %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")

fung.taxa[is.na(fung.taxa)] <- ""

for (i in 1:7){ fung.taxa[,i] <- as.character(fung.taxa[,i])}
####### Fills holes in the tax table
fung.taxa[is.na(fung.taxa)] <- ""
for (i in 1:nrow(fung.taxa)){
  if (fung.taxa[i,2] == ""){
    kingdom <- paste("Unidentified_", fung.taxa[i,1], sep = "")
    fung.taxa[i, 2:7] <- kingdom
  } else if (fung.taxa[i,3] == ""){
    phylum <- paste("Unidentified_", fung.taxa[i,2], sep = "")
    fung.taxa[i, 3:7] <- phylum
  } else if (fung.taxa[i,4] == ""){
    class <- paste("Unidentified_", fung.taxa[i,3], sep = "")
    fung.taxa[i, 4:7] <- class
  } else if (fung.taxa[i,5] == ""){
    order <- paste("Unidentified_", fung.taxa[i,4], sep = "")
    fung.taxa[i, 5:7] <- order
  } else if (fung.taxa[i,6] == ""){
    family <- paste("Unidentified_", fung.taxa[i,5], sep = "")
    fung.taxa[i, 6:7] <- family
  } else if (fung.taxa[i,7] == ""){
    fung.taxa$Species[i] <- paste("Unidentified_",fung.taxa$Genus[i], sep = "")
  }
}

fung <- fung[which(row.names(fung) %in% row.names(fung.taxa)),]
fung.taxa <- as.matrix(fung.taxa)

#Load metadata and the rest of this is all the same
fung.sort <- fung[,order(colnames(fung))]
md.fung <- roughmetadata_fung[which(unlist(roughmetadata_fung$Name) %in% colnames(fung.sort)),]
#Summarizes read counts of samples, useful sanity check
summary(colSums(fung.sort))

md.fung <- md.fung %>%
  tibble::column_to_rownames("Name")

#Create phyloseq object, there is no read count or prevalence filtering this time

OTU_decontam_rough = otu_table(fung.sort, taxa_are_rows = TRUE)
TAX_decontam_rough = tax_table(fung.taxa)
samples_decontam_rough = sample_data(md.fung)

exp2_fung_decontam_rough <- phyloseq(OTU_decontam_rough, TAX_decontam_rough, samples_decontam_rough)
exp2_fung_gen_decontam_rough <- tax_glom(exp2_fung_decontam_rough, taxrank = "Genus")


#Identify top 49 genera and rename everything else to "Other"
top20_fung_list <- names(sort(taxa_sums(exp2_fung_gen_decontam_rough), decreasing=TRUE)[1:49])
top20_fung_decontam_rel <- transform_sample_counts(exp2_fung_gen_decontam_rough, function(x) x / sum(x) )
top20_fung_decontam_df <- psmelt(top20_fung_decontam_rel)
top20_fung_decontam_df[!(top20_fung_decontam_df$OTU %in% top20_fung_decontam_list),]$Genus <- 'Other'


###Barplot of top 50 fungal genera###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(50)

barplot_gen_decontam_its <- ggplot(top20_fung_decontam_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_decontam_its + scale_fill_manual(values = barplot_colors)
#Save



### Fig. S2D: 16S Decontaminated Taxa Barplot ###

#Read in filtered ASV table, table of removed taxa, and metadata
bact <- read_csv("Human_16S_ASVs_Decontaminated.csv")
bact.taxa <- read_csv("Human_16S_Taxa2.csv")
roughmetadata_bact <- read_csv("Human_Metadata2.csv")

#Remove samples with unknown BPD severity if not using the metadata from GitHub
#roughmetadata_bact <- roughmetadata_bact[!(is.na(roughmetadata_bact$BPD_Severity)),]

#### Decontaminated 16S Data Prep ####

#Load ASV table
bact <- bact[rowSums(is.na(bact)) == 0, ]
bact <- bact %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("Name")

bact.taxa[is.na(bact.taxa)] <- ""

for (i in 1:6){ bact.taxa[,i] <- as.character(bact.taxa[,i])}
####### Fille holes in the tax table
bact.taxa[is.na(bact.taxa)] <- ""
for (i in 1:nrow(bact.taxa)){
  if (bact.taxa[i,2] == ""){
    kingdom <- paste("Unidentified_", bact.taxa[i,1], sep = "")
    bact.taxa[i, 2:7] <- kingdom
  } else if (bact.taxa[i,3] == ""){
    phylum <- paste("Unidentified_", bact.taxa[i,2], sep = "")
    bact.taxa[i, 3:7] <- phylum
  } else if (bact.taxa[i,4] == ""){
    class <- paste("Unidentified_", bact.taxa[i,3], sep = "")
    bact.taxa[i, 4:7] <- class
  } else if (bact.taxa[i,5] == ""){
    order <- paste("Unidentified_", bact.taxa[i,4], sep = "")
    bact.taxa[i, 5:7] <- order
  } else if (bact.taxa[i,6] == ""){
    family <- paste("Unidentified_", bact.taxa[i,5], sep = "")
    bact.taxa[i, 6:7] <- family
  } 
}

bact.taxa <- as.matrix(bact.taxa[,c(1:6)])
bact <- bact[which(row.names(bact) %in% row.names(bact.taxa)),]

#Load metadata and the rest of this is all the same
bact.sort <- bact[,order(colnames(bact))]
md.bact <- roughmetadata_bact[which(unlist(roughmetadata_bact$Name) %in% colnames(bact.sort)),]
#Summarizes read counts of samples, useful sanity check
summary(colSums(bact.sort))

md.bact <- md.bact %>%
  tibble::column_to_rownames("Name")

#Create phyloseq object, there is no read count or prevalence filtering this time

OTU_decontam_rough = otu_table(bact.sort, taxa_are_rows = TRUE)
TAX_decontam_rough = tax_table(bact.taxa)
samples_decontam_rough = sample_data(md.bact)

exp2_bact_decontam_rough <- phyloseq(OTU_decontam_rough, TAX_decontam_rough, samples_decontam_rough)
exp2_bact_gen_decontam_rough <- tax_glom(exp2_bact_decontam_rough, taxrank = "Genus")


#Identify top 49 genera and rename everything else to "Other"
top20_bact_list <- names(sort(taxa_sums(exp2_bact_gen_decontam_rough), decreasing=TRUE)[1:49])
top20_bact_decontam_rel <- transform_sample_counts(exp2_bact_gen_decontam_rough, function(x) x / sum(x) )
top20_bact_decontam_df <- psmelt(top20_bact_decontam_rel)
top20_bact_decontam_df[!(top20_bact_decontam_df$OTU %in% top20_bact_decontam_list),]$Genus <- 'Other'


###Barplot of top 50 bacterial genera###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(50)

barplot_gen_decontam_16s <- ggplot(top20_bact_decontam_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_decontam_16s + scale_fill_manual(values = barplot_colors)
#Save
