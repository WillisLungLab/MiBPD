library(tidyverse)
library(PERFect)
library(phyloseq)

#Read in ASV table, taxonomy, and metadata
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
write.csv(bact_filt2, "MiBPD_16S_PERFect_ASVs_Perm2.csv")

removed <- row.names(bact_nozero[!(row.names(bact_nozero) %in% row.names(bact_filt)),])
bact.taxa <- read_csv("Human_16s_Taxa.csv")
bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("Name")
removed.taxa <- bact.taxa[row.names(bact.taxa) %in% removed,]
write.csv(removed.taxa, "PERFect_16s_Removed_Taxa.csv")


###Contamination Barplot###

#Read in unfiltered ASV table, table of removed taxa, and metadata

#Remove samples with unknown BPD severity since these samples were excluded from analysis
roughmetadata_bact <- roughmetadata_bact[!(is.na(roughmetadata_bact$BPD_Severity)),]

####16S Data Prep####

#Load ASV table
bact <- bact[rowSums(is.na(bact)) == 0, ]
bact <- bact %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("...1")

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

