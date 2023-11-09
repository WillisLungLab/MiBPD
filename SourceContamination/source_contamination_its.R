library(tidyverse)
library(PERFect)
library(phyloseq)
library(plyr)
library(RColorBrewer)

fung <- read_csv("Human_ITS_ASVs.csv")

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

#Get table of removed taxa
removed <- colnames(fung_nozero[,!(colnames(fung_nozero) %in% row.names(fung_filt))])
fung.taxa <- read_csv("Human_ITS_Taxa.csv")
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")
removed.taxa <- fung.taxa[row.names(fung.taxa) %in% removed,]
#Save as .csv


###Contamination Barplot###

#Read in unfiltered ASV table, table of removed taxa, and metadata

#Remove samples with unknown BPD severity since these samples were excluded from analysis
roughmetadata_fung <- roughmetadata_fung[!(is.na(roughmetadata_fung$BPD_Severity)),]

####ITS Data Prep####

#Load ASV table
fung <- fung[rowSums(is.na(fung)) == 0, ]
fung <- fung %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("...1")

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
