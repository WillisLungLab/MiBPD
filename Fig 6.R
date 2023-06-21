library(tidyverse)
library(readxl)
library(vegan)
library(phyloseq)
library(DESeq2)
library(indicspecies)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(factoextra)
library(ecole)

#############Fig. 6A: Gut ITS#############

###Data Prep###
#Load ASV table
fung <- fung %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")
fung.taxa <- as.matrix(fung.taxa)

#Load metadata
fung.sort <- fung[,order(colnames(fung))]
md.fung <- roughmetadata_its[which(unlist(roughmetadata_its$Name) %in% colnames(fung.sort)),]
summary(colSums(fung.sort[,-1]))

#Create phyloseq object with unfiltered data for alpha diversity
OTU_rough = otu_table(fung.sort, taxa_are_rows = TRUE)
TAX_rough = tax_table(fung.taxa)
md.fung <- md.fung %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.fung)
exp2_fung_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

fung.sort2 <- fung.sort[,-c(1,which(colSums(fung.sort[,-1])<50))]
md.fung.sort2 <- md.fung[which(row.names(md.fung) %in% colnames(fung.sort2)),]

#This is the same
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

fung.nsampls <- apply(fung.sort, 1, gt0)
fung.clean <- fung.sort2[which(fung.nsampls>1),]


#Setup to create the phyloseq object
OTU = otu_table(fung.clean, taxa_are_rows = TRUE)
TAX = tax_table(fung.taxa)
samples = sample_data(md.fung.sort2)

#Convert!
exp2_fung <- phyloseq(OTU, TAX, samples)
exp2_fung_prev <- filter_taxa(exp2_fung , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)

###Alpha Diversity###
fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
richness_est_fung <- estimate_richness(fr_fung, measures = c("Chao1", "Shannon"))

richness_est_fung <- richness_est_fung %>%
  mutate(
    Organ = md.fung$Organ,
    Oxygen = md.fung$Oxygen,
    FMT = md.fung$FMT,
    Condition  = md.fung$Condition
  )
#Save as .csv

###Beta Diversity###

exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)

exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
exp2_fung_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_fung_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Mouse Gut ITS PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.6), c(-0.6, 0.6), font = 2, font.lab = 2, xlab="PC1(20.23% Explained)", ylab="PC2(14.35% Explained)", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], factor2_exp2, kind = "se", conf = 0.95, col = pca_colors)
dev.off()

permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_fung_pcoa <- permanova_pairwise(exp2_fung_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(pairwise_permanova_fung_pcoa)
print(permdisp_fung_pcoa)


#############Fig. 6B: Gut Multikingdom#############

###Data Prep###

#Load ASV table
combined <- combined %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
combined.taxa <- combined.taxa %>%
  tibble::column_to_rownames("Name")
combined.taxa <- as.matrix(combined.taxa)

#Load metadata file
combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined[which(unlist(roughmetadata_combined$Name) %in% colnames(combined.sort)),]
summary(colSums(combined.sort[,-1]))

OTU_rough = otu_table(combined.sort, taxa_are_rows = TRUE)
TAX_rough = tax_table(combined.taxa)
md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.combined)
exp2_combined_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

combined.sort2 <- combined.sort[,-c(1,which(colSums(combined.sort[,-1])<50))]
md.combined.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.sort2)),]

combined.nsampls <- apply(combined.sort, 1, gt0)
combined.clean <- combined.sort2[which(combined.nsampls>1),]

#Setup to create the phyloseq object
OTU = otu_table(combined.clean, taxa_are_rows = TRUE)
TAX = tax_table(combined.taxa)
samples = sample_data(md.combined.sort2)

#Convert!
exp2_combined <- phyloseq(OTU, TAX, samples)
exp2_combined_prev <- filter_taxa(exp2_combined , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)


###Alpha Diversity###

fr_combined <- prune_taxa(taxa_sums(exp2_combined_rough) > 0, exp2_combined_rough)
richness_est_combined <- estimate_richness(fr_combined, measures = c("Chao1", "Shannon"))

richness_est_combined <- richness_est_combined %>%
  mutate(
    Organ = md.combined$Organ,
    Oxygen = md.combined$Oxygen,
    FMT = md.combined$FMT,
    Condition  = md.combined$Condition
  )
#Save as .csv

###Beta Diversity###

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)

exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
exp2_combined_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_combined_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Mouse Gut combined PCoA New CI.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.6), c(-0.7, 0.4), font = 2, font.lab = 2, xlab="PC1(17.36% Explained)", ylab="PC2(12.60% Explained)", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], factor2_exp2, kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_combined_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

permanova_pcoa_df <- data.frame(exp2_combined_meta_rel)
set.seed(1312)
permanova_combined_pcoa <- vegan::adonis2(exp2_combined_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_pcoa <- permanova_pairwise(exp2_combined_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

combined_pcoa_dist <- vegdist(exp2_combined_otu_rel, method = "bray")
disp_combined_pcoa <- betadisper(combined_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_combined_pcoa <- permutest(disp_combined_pcoa, permutations = 10000)

print(permanova_combined_pcoa)
print(pairwise_permanova_combined_pcoa)
print(permdisp_combined_pcoa)


#############Fig. 6C: Lung ITS#############

###Data Prep###

#Load ASV table
fung <- fung %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")
fung.taxa <- as.matrix(fung.taxa)

#Load metadata
fung.sort <- fung[,order(colnames(fung))]
md.fung <- roughmetadata_its[which(unlist(roughmetadata_its$Name) %in% colnames(fung.sort)),]
summary(colSums(fung.sort[,-1]))


OTU_rough = otu_table(fung.sort, taxa_are_rows = TRUE)

TAX_rough = tax_table(fung.taxa)

md.fung <- md.fung %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.fung)

exp2_fung_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

fung.sort2 <- fung.sort[,-c(1,which(colSums(fung.sort[,-1])<0))]
md.fung.sort2 <- md.fung[which(row.names(md.fung) %in% colnames(fung.sort2)),]

fung.nsampls <- apply(fung.sort, 1, gt0)
fung.clean <- fung.sort2[which(fung.nsampls>1),]
fung.clean.colsum <- fung.clean[,colSums(fung.clean)>0]


#Setup to create the phyloseq object
OTU = otu_table(fung.clean.colsum, taxa_are_rows = TRUE)
TAX = tax_table(fung.taxa)
samples = sample_data(md.fung.sort2)

#Convert!
exp2_fung <- phyloseq(OTU, TAX, samples)
exp2_fung_prev <- filter_taxa(exp2_fung , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)


###Alpha Diversity###

fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
richness_est_fung <- estimate_richness(fr_fung, measures = c("Chao1", "Shannon"))

richness_est_fung <- richness_est_fung %>%
  mutate(
    Organ = md.fung$Organ,
    Oxygen = md.fung$Oxygen,
    FMT = md.fung$FMT,
    Condition  = md.fung$Condition
  )
#Save

###Beta Diversity###

exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)

exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
exp2_fung_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_fung_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Mouse Lung ITS PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.7, 0.5), c(-0.6, 0.6), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95,factor2_exp2, col = pca_colors)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_fung_pcoa <- permanova_pairwise(exp2_fung_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(pairwise_permanova_fung_pcoa)
print(permdisp_fung_pcoa)


#############Fig. 6D: Lung Multikingdom#############

###Data Prep###

#Load ASV table
combined <- combined %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
combined.taxa <- combined.taxa %>%
  tibble::column_to_rownames("Name")
combined.taxa <- as.matrix(combined.taxa)

#Load metadata
combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined[which(unlist(roughmetadata_combined$Name) %in% colnames(combined.sort)),]
summary(colSums(combined.sort[,-1]))

OTU_rough = otu_table(combined.sort, taxa_are_rows = TRUE)
TAX_rough = tax_table(combined.taxa)
md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.combined)
exp2_combined_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

combined.sort2 <- combined.sort[,-c(1,which(colSums(combined.sort[,-1])<0))]
md.combined.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.sort2)),]

combined.nsampls <- apply(combined.sort, 1, gt0)
combined.clean <- combined.sort2[which(combined.nsampls>1),]

#Setup to create the phyloseq object
OTU = otu_table(combined.clean, taxa_are_rows = TRUE)
TAX = tax_table(combined.taxa)
samples = sample_data(md.combined.sort2)

#Convert!
exp2_combined <- phyloseq(OTU, TAX, samples)
exp2_combined_prev <- filter_taxa(exp2_combined , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)


###Alpha Diversity###

fr_combined <- prune_taxa(taxa_sums(exp2_combined_rough) > 0, exp2_combined_rough)
richness_est_combined <- estimate_richness(fr_combined, measures = c("Chao1", "Shannon"))

richness_est_combined <- richness_est_combined %>%
  mutate(
    Organ = md.combined$Organ,
    Oxygen = md.combined$Oxygen,
    FMT = md.combined$FMT,
    Condition  = md.combined$Condition
  )
#Save

###Beta Diversity###

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)

exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
exp2_combined_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_combined_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Mouse lung combined PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.7, 0.5), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1(17.36% Explained)", ylab="PC2(12.60% Explained)", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], factor2_exp2,kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_combined_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

permanova_pcoa_df <- data.frame(exp2_combined_meta_rel)
set.seed(1312)
permanova_combined_pcoa <- vegan::adonis2(exp2_combined_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_pcoa <- permanova_pairwise(exp2_combined_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

combined_pcoa_dist <- vegdist(exp2_combined_otu_rel, method = "bray")
disp_combined_pcoa <- betadisper(combined_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_combined_pcoa <- permutest(disp_combined_pcoa, permutations = 10000)

print(permanova_combined_pcoa)
print(pairwise_permanova_combined_pcoa)
print(permdisp_combined_pcoa)
