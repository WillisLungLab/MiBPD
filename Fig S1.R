library(tidyverse)
library(PERFect)
library(phyloseq)
library(plyr)
library(RColorBrewer)

###Fig. S1A: ITS Filtering Loss Plot###

#Read in non-filtered ASVs

fung <- fung %>%
  tibble::column_to_rownames("Name")
fung <- as.data.frame(t(fung))

fung_nozero <- fung[,colSums(fung)>0]

#Get p-values for permutation filtering
set.seed(1312)
res_sim <- PERFect_sim(X=fung_nozero)
fung_filt1 <- as.data.frame(t(res_sim[["filtX"]]))
new_order1 = sort(row.names(fung_filt1))

#Permutation filter
set.seed(1312)
res_perm2 <- PERFect_perm(X = fung_nozero, Order = "pvals",pvals_sim = res_sim, algorithm = "fast", alpha = 0.05, rollmean = FALSE)
fung_filt2 <- as.data.frame(t(res_perm2[["filtX"]]))
new_order2 = sort(row.names(fung_filt2))
fung_filt2 <- as.data.frame(fung_filt2[new_order2,])

#Filtering loss plot
pvals_Plots(PERFect = res_perm2, X = fung_nozero, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
#Save

#Clean up and save results
fung_filt <- as.data.frame(t(res_perm2[["filtX"]]))
new_order = sort(row.names(fung_filt))
fung_filt <- as.data.frame(fung_filt[new_order,])
#Save as .csv

#Get table of filtered taxa
unremoved <- colnames(fung_nozero[,colnames(fung_nozero) %in% row.names(fung_filt)])
#Read in unfiltered taxa
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")
unremoved.taxa <- fung.taxa[row.names(fung.taxa) %in% unremoved,]
#Save as .csv


###Fig. S1B: 16S Filtering Loss Plot###

#Read in non-filtered ASVs, taxa, and metadata
meta <- meta %>%
  tibble::column_to_rownames("Name")

bact_taxa <- bact_taxa %>%
  tibble::column_to_rownames("...1")
bact_taxa <- as.matrix(bact_taxa)

bact <- bact %>%
  tibble::column_to_rownames("...1")
bact <- as.data.frame(t(bact))
colnames(bact) %in% row.names(meta)

OTU = otu_table(bact_nozero, taxa_are_rows = TRUE)
TAX = tax_table(bact_taxa)
samples = sample_data(meta)
rough <- phyloseq(OTU, TAX, samples)

#Removes anything unassigned at the Phylum level and some especially low abundance taxa, this was required for PERFect to run without crashing
rough_phy = subset_taxa(rough, !(is.na(Phylum)))
rough_prev <- filter_taxa(rough_phy , function(x) sum(x >= 1) > (0.02*length(x)), TRUE)
otu_for_perfect <- as.data.frame(t(otu_table(rough_prev)))

set.seed(1312)
res_sim <- PERFect_sim(X=otu_for_perfect)

set.seed(1312)
res_perm2 <- PERFect_perm(X = otu_for_perfect, Order = "pvals",pvals_sim = res_sim, algorithm = "fast", alpha = 0.05, rollmean = FALSE)
pvals_Plots(PERFect = res_perm2, X = otu_for_perfect, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

bact_filt2 <- as.data.frame(t(res_perm2[["filtX"]]))
new_order2 = sort(row.names(bact_filt2))
bact_filt2 <- as.data.frame(bact_filt2[new_order2,])
#Save as .csv

#Get table of filtered taxa
unremoved <- row.names(bact_nozero[row.names(bact_nozero) %in% row.names(bact_filt),])
bact.taxa <- read_csv("Human_16s_Taxa.csv")
bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("Name")
unremoved.taxa <- bact.taxa[row.names(bact.taxa) %in% unremoved,]
#Save as .csv



###Fig. S1C: ITS Filtered Taxa Barplot###

#Read in filtered ASV table, table of unremoved taxa, and metadata
fung <- read_csv("Human_ITS_ASVs_Decontaminated.csv")
fung.taxa <- read_csv("Human_ITS_Taxa2.csv")
roughmetadata_fung <- read_csv("Human_Metadata2.csv")

#Remove samples with unknown BPD severity since these samples were excluded from analysis (already removed in the public metadata, uncomment if using the full metadata)
#roughmetadata_fung <- roughmetadata_fung[!(is.na(roughmetadata_fung$BPD_Severity)),]

####ITS Data Prep####

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


OTU_contam_rough = otu_table(fung.sort, taxa_are_rows = TRUE)
TAX_contam_rough = tax_table(fung.taxa)
samples_contam_rough = sample_data(md.fung)

exp2_fung_contam_rough <- phyloseq(OTU_contam_rough, TAX_contam_rough, samples_contam_rough)
exp2_fung_gen_contam_rough <- tax_glom(exp2_fung_contam_rough, taxrank = "Genus")

#Identify top 49 genera and rename everything else to "Other"
top20_fung_merged_list <- names(sort(taxa_sums(exp2_fung_gen_contam_rough), decreasing=TRUE)[1:49])
top20_merged_fung_rel <- transform_sample_counts(exp2_fung_gen_contam_rough, function(x) x / sum(x) )
top20_merged_fung_df <- psmelt(top20_merged_fung_rel)
top20_merged_fung_df[!(top20_merged_fung_df$OTU %in% top20_fung_merged_list),]$Genus <- 'Other'


###Barplot of top 50 fungal genera###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(50)

barplot_gen_contam_its <- ggplot(top20_merged_fung_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_contam_its + scale_fill_manual(values = barplot_colors)
#Save


###Fig. S1D: 16S Filtering Loss Plot###

#Read in filtered ASV table, table of removed taxa, and metadata
bact <- read_csv("Human_16S_ASVs_Decontaminated.csv")
bact.taxa <- read_csv("Human_16S_Taxa2.csv")
roughmetadata_bact <- read_csv("Human_Metadata2.csv")

#Remove samples with unknown BPD severity if not using the metadata from GitHub
#roughmetadata_bact <- roughmetadata_bact[!(is.na(roughmetadata_bact$BPD_Severity)),]

####16S Data Prep####

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


OTU_contam_rough = otu_table(bact.sort, taxa_are_rows = TRUE)
TAX_contam_rough = tax_table(bact.taxa)
samples_contam_rough = sample_data(md.bact)

exp2_bact_contam_rough <- phyloseq(OTU_contam_rough, TAX_contam_rough, samples_contam_rough)
exp2_bact_gen_contam_rough <- tax_glom(exp2_bact_contam_rough, taxrank = "Genus")

#Identify top 49 genera and rename everything else to "Other"
top20_bact_merged_list <- names(sort(taxa_sums(exp2_bact_gen_contam_rough), decreasing=TRUE)[1:49])
top20_merged_bact_rel <- transform_sample_counts(exp2_bact_gen_contam_rough, function(x) x / sum(x) )
top20_merged_bact_df <- psmelt(top20_merged_bact_rel)
top20_merged_bact_df[!(top20_merged_bact_df$OTU %in% top20_bact_merged_list),]$Genus <- 'Other'


###Barplot of top 50 bacterial genera###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(50)

barplot_gen_contam_16s <- ggplot(top20_merged_bact_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_contam_16s + scale_fill_manual(values = barplot_colors)
#Save
