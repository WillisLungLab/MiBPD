library(tidyverse)
library(readxl)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(factoextra)
library(ecole)

#############Fig. 6A: Gut ITS#############

fung_mouse <- read_csv("Mouse_Gut_ITS_ASVs.csv")
fung.mouse.taxa <- read_csv("Mouse_Gut_ITS_Taxa.csv")
roughmetadata_fung_mouse <- read_csv("Mouse_Gut_Metadata.csv")
###Data Prep###
#Load ASV table
fung_mouse <- fung_mouse %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.mouse.taxa <- fung.mouse.taxa %>%
  tibble::column_to_rownames("Name")
fung.mouse.taxa <- as.matrix(fung.mouse.taxa)

#Load metadata
fung.mouse.sort <- fung_mouse[,order(colnames(fung_mouse))]
md.fung.mouse <- roughmetadata_fung_mouse[which(unlist(roughmetadata_fung_mouse$Name) %in% colnames(fung.mouse.sort)),]
summary(colSums(fung.mouse.sort))

#Create phyloseq object with unfiltered data for alpha diversity
OTU_rough = otu_table(fung.mouse.sort, taxa_are_rows = TRUE)
TAX_rough = tax_table(fung.mouse.taxa)
md.fung.mouse <- md.fung.mouse %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.fung.mouse)
exp2_fung_mouse_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

fung.mouse.sort2 <- fung.mouse.sort[,c(which(colSums(fung.mouse.sort)>=50))]
md.fung.mouse.sort2 <- md.fung.mouse[which(row.names(md.fung.mouse) %in% colnames(fung.mouse.sort2)),]

#This is the same
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

fung_mouse.nsampls <- apply(fung.mouse.sort2, 1, gt0)
fung_mouse.clean <- fung.mouse.sort2[which(fung_mouse.nsampls>1),]


#Setup to create the phyloseq object
OTU = otu_table(fung_mouse.clean, taxa_are_rows = TRUE)
TAX = tax_table(fung.mouse.taxa)
samples = sample_data(md.fung.mouse.sort2)

#Convert!
exp2_fung_mouse <- phyloseq(OTU, TAX, samples)
exp2_fung_mouse_prev <- filter_taxa(exp2_fung_mouse , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)

###Alpha Diversity###
fr_fung_mouse <- prune_taxa(taxa_sums(exp2_fung_mouse_rough) > 0, exp2_fung_mouse_rough)
richness_est_fung_mouse <- estimate_richness(fr_fung_mouse, measures = c("Simpson", "Shannon"))

wilcox_alpha_fung_mouse <- t(sapply(richness_est_fung_mouse, function(x) unlist(kruskal.test(x~sample_data(fr_fung_mouse)$Condition)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_fung_mouse

richness_est_fung_mouse <- richness_est_fung_mouse %>%
  mutate(
    Organ = md.fung.mouse$Organ,
    Oxygen = md.fung.mouse$Oxygen,
    FMT = md.fung.mouse$FMT,
    Condition  = md.fung.mouse$Condition
  )

#Save as .csv

###Beta Diversity###

exp2_fung_mouse_rel_prev <- transform_sample_counts(exp2_fung_mouse_prev, function(x) x / sum(x) )

exp2_fung_mouse_otu_rel <- as.data.frame(t(exp2_fung_mouse_rel_prev@otu_table))
exp2_fung_mouse_tax_rel <- as.data.frame(exp2_fung_mouse_rel_prev@tax_table)
exp2_fung_mouse_meta_rel <- as.data.frame(exp2_fung_mouse_rel_prev@sam_data)

exp2_fung_mouse_rel_bray = vegdist(exp2_fung_mouse_otu_rel, method='bray', na.rm = TRUE)
exp2_fung_mouse_rel_pcoa <- ape::pcoa(exp2_fung_mouse_rel_bray)
exp2_fung_mouse_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_fung_mouse_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Fig 6A PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.6), c(-0.6, 0.6), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_mouse_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordispider(exp2_fung_mouse_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
ordiellipse(exp2_fung_mouse_rel_pcoa$vectors[,1:2], factor2_exp2, kind = "se", conf = 0.95, col = pca_colors)
dev.off()

permanova_pcoa_df <- data.frame(exp2_fung_mouse_meta_rel)
set.seed(1312)
permanova_fung_mouse_pcoa <- vegan::adonis2(exp2_fung_mouse_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_fung_mouse_pcoa <- permanova_pairwise(exp2_fung_mouse_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

fung_mouse_pcoa_dist <- vegdist(exp2_fung_mouse_otu_rel, method = "bray")
disp_fung_mouse_pcoa <- betadisper(fung_mouse_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_fung_mouse_pcoa <- permutest(disp_fung_mouse_pcoa, permutations = 10000)

print(permanova_fung_mouse_pcoa)
print(pairwise_permanova_fung_mouse_pcoa)
print(permdisp_fung_mouse_pcoa)


#############Fig. 6B: Gut 16S#############

###Data Prep###

bact_mouse <- read_csv("Mouse_Gut_16s_ASVs.csv")
bact.mouse.taxa <- read_csv("Mouse_Gut_16s_Taxa.csv")
roughmetadata_bact_mouse <- read_csv("Mouse_Gut_Metadata.csv")

#Load ASV table
bact_mouse <- bact_mouse %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
bact.mouse.taxa <- bact.mouse.taxa %>%
  tibble::column_to_rownames("Name")
bact.mouse.taxa <- as.matrix(bact.mouse.taxa)

#Load metadata file
bact_mouse.sort <- bact_mouse[,order(colnames(bact_mouse))]
md.bact_mouse <- roughmetadata_bact_mouse[which(unlist(roughmetadata_bact_mouse$Name) %in% colnames(bact_mouse.sort)),]
summary(colSums(bact_mouse.sort))

OTU_rough = otu_table(bact_mouse.sort, taxa_are_rows = TRUE)
TAX_rough = tax_table(bact.mouse.taxa)
md.bact_mouse <- md.bact_mouse %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.bact_mouse)
exp2_bact_mouse_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

bact_mouse.sort2 <- bact_mouse.sort[,c(which(colSums(bact_mouse.sort)>=1000))]
md.bact_mouse.sort2 <- md.bact_mouse[which(row.names(md.bact_mouse) %in% colnames(bact_mouse.sort2)),]

bact_mouse.nsampls <- apply(bact_mouse.sort, 1, gt0)
bact_mouse.clean <- bact_mouse.sort2[which(bact_mouse.nsampls>1),]

#Setup to create the phyloseq object
OTU = otu_table(bact_mouse.clean, taxa_are_rows = TRUE)
TAX = tax_table(bact.mouse.taxa)
samples = sample_data(md.bact_mouse.sort2)

#Convert!
exp2_bact_mouse <- phyloseq(OTU, TAX, samples)
exp2_bact_mouse_prev <- filter_taxa(exp2_bact_mouse , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)


###Alpha Diversity###

fr_bact_mouse <- prune_taxa(taxa_sums(exp2_bact_mouse_rough) > 0, exp2_bact_mouse_rough)
richness_est_bact_mouse <- estimate_richness(fr_bact_mouse, measures = c("Simpson", "Shannon"))

wilcox_alpha_bact_mouse <- t(sapply(richness_est_bact_mouse, function(x) unlist(kruskal.test(x~sample_data(fr_bact_mouse)$Condition)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_bact_mouse

richness_est_bact_mouse <- richness_est_bact_mouse %>%
  mutate(
    Organ = md.bact_mouse$Organ,
    Oxygen = md.bact_mouse$Oxygen,
    FMT = md.bact_mouse$FMT,
    Condition  = md.bact_mouse$Condition
  )
#Save as .csv
write.csv(richness_est_bact_mouse, "Fig 6A Gut 16S Alpha Diversity.csv")

###Beta Diversity###

exp2_bact_mouse_rel_prev <- transform_sample_counts(exp2_bact_mouse_prev, function(x) x / sum(x) )

exp2_bact_mouse_otu_rel <- as.data.frame(t(exp2_bact_mouse_rel_prev@otu_table))
exp2_bact_mouse_tax_rel <- as.data.frame(exp2_bact_mouse_rel_prev@tax_table)
exp2_bact_mouse_meta_rel <- as.data.frame(exp2_bact_mouse_rel_prev@sam_data)

exp2_bact_mouse_rel_bray = vegdist(exp2_bact_mouse_otu_rel, method='bray')
exp2_bact_mouse_rel_pcoa <- ape::pcoa(exp2_bact_mouse_rel_bray)
exp2_bact_mouse_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_bact_mouse_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Fig 6B PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.75, 0.4), c(-0.6, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_bact_mouse_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordiellipse(exp2_bact_mouse_rel_pcoa$vectors[,1:2], factor2_exp2, kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_bact_mouse_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

permanova_pcoa_df <- data.frame(exp2_bact_mouse_meta_rel)
set.seed(1312)
permanova_bact_mouse_pcoa <- vegan::adonis2(exp2_bact_mouse_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_bact_mouse_pcoa <- permanova_pairwise(exp2_bact_mouse_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

bact_mouse_pcoa_dist <- vegdist(exp2_bact_mouse_otu_rel, method = "bray")
disp_bact_mouse_pcoa <- betadisper(bact_mouse_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_bact_mouse_pcoa <- permutest(disp_bact_mouse_pcoa, permutations = 10000)

print(permanova_bact_mouse_pcoa)
print(pairwise_permanova_bact_mouse_pcoa)
print(permdisp_bact_mouse_pcoa)


#############Fig. 6C: Lung ITS#############

fung_mouse <- read_csv("Mouse_Lung_ITS_ASVs.csv")
fung.mouse.taxa <- read_csv("Mouse_Lung_ITS_Taxa.csv")
roughmetadata_fung_mouse <- read_csv("Mouse_Lung_Metadata.csv")

###Data Prep###

#Load ASV table
fung_mouse <- fung_mouse %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
fung.mouse.taxa <- fung.mouse.taxa %>%
  tibble::column_to_rownames("Name")
fung.mouse.taxa <- as.matrix(fung.mouse.taxa)

#Load metadata
fung.mouse.sort <- fung_mouse[,order(colnames(fung_mouse))]
md.fung.mouse <- roughmetadata_fung_mouse[which(unlist(roughmetadata_fung_mouse$Name) %in% colnames(fung.mouse.sort)),]
summary(colSums(fung.mouse.sort))


OTU_rough = otu_table(fung.mouse.sort, taxa_are_rows = TRUE)

TAX_rough = tax_table(fung.mouse.taxa)

md.fung.mouse <- md.fung.mouse %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.fung.mouse)

exp2_fung_mouse_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

fung.mouse.sort2 <- fung.mouse.sort[,c(which(colSums(fung.mouse.sort)>=0))]
md.fung.mouse.sort2 <- md.fung.mouse[which(row.names(md.fung.mouse) %in% colnames(fung.mouse.sort2)),]

fung_mouse.nsampls <- apply(fung.mouse.sort, 1, gt0)
fung_mouse.clean <- fung.mouse.sort2[which(fung_mouse.nsampls>1),]
fung_mouse.clean.colsum <- fung_mouse.clean[,colSums(fung_mouse.clean)>0]


#Setup to create the phyloseq object
OTU = otu_table(fung_mouse.clean.colsum, taxa_are_rows = TRUE)
TAX = tax_table(fung.mouse.taxa)
samples = sample_data(md.fung.mouse.sort2)

#Convert!
exp2_fung_mouse <- phyloseq(OTU, TAX, samples)
exp2_fung_mouse_prev <- filter_taxa(exp2_fung_mouse , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)


###Alpha Diversity###

fr_fung_mouse <- prune_taxa(taxa_sums(exp2_fung_mouse_rough) > 0, exp2_fung_mouse_rough)
richness_est_fung_mouse <- estimate_richness(fr_fung_mouse, measures = c("Simpson", "Shannon"))

richness_est_fung_mouse <- richness_est_fung_mouse %>%
  mutate(
    Organ = md.fung.mouse$Organ,
    Oxygen = md.fung.mouse$Oxygen,
    FMT = md.fung.mouse$FMT,
    Condition  = md.fung.mouse$Condition
  )
#Save

###Beta Diversity###

exp2_fung_mouse_rel_prev <- transform_sample_counts(exp2_fung_mouse_prev, function(x) x / sum(x) )

exp2_fung_mouse_otu_rel <- as.data.frame(t(exp2_fung_mouse_rel_prev@otu_table))
exp2_fung_mouse_tax_rel <- as.data.frame(exp2_fung_mouse_rel_prev@tax_table)
exp2_fung_mouse_meta_rel <- as.data.frame(exp2_fung_mouse_rel_prev@sam_data)

exp2_fung_mouse_rel_bray = vegdist(exp2_fung_mouse_otu_rel, method='bray')
exp2_fung_mouse_rel_pcoa <- ape::pcoa(exp2_fung_mouse_rel_bray)
exp2_fung_mouse_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_fung_mouse_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Fig 6C PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.7), c(-0.6, 0.6), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_mouse_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordiellipse(exp2_fung_mouse_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95,factor2_exp2, col = pca_colors)
ordispider(exp2_fung_mouse_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

permanova_pcoa_df <- data.frame(exp2_fung_mouse_meta_rel)
set.seed(1312)
permanova_fung_mouse_pcoa <- vegan::adonis2(exp2_fung_mouse_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_fung_mouse_pcoa <- permanova_pairwise(exp2_fung_mouse_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

fung_mouse_pcoa_dist <- vegdist(exp2_fung_mouse_otu_rel, method = "bray")
disp_fung_mouse_pcoa <- betadisper(fung_mouse_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_fung_mouse_pcoa <- permutest(disp_fung_mouse_pcoa, permutations = 10000)

print(permanova_fung_mouse_pcoa)
print(pairwise_permanova_fung_mouse_pcoa)
print(permdisp_fung_mouse_pcoa)


#############Fig. 6D: Lung 16S#############

###Data Prep###
bact_mouse <- read_csv("Mouse_Lung_16s_ASVs.csv")
bact.mouse.taxa <- read_csv("Mouse_Lung_16s_Taxa.csv")
roughmetadata_bact_mouse <- read_csv("Mouse_Lung_Metadata.csv")

#Load ASV table
bact_mouse <- bact_mouse %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
bact.mouse.taxa <- bact.mouse.taxa %>%
  tibble::column_to_rownames("Name")
bact.mouse.taxa <- as.matrix(bact.mouse.taxa)

#Load metadata
bact_mouse.sort <- bact_mouse[,order(colnames(bact_mouse))]
md.bact_mouse <- roughmetadata_bact_mouse[which(unlist(roughmetadata_bact_mouse$Name) %in% colnames(bact_mouse.sort)),]
summary(colSums(bact_mouse.sort))

OTU_rough = otu_table(bact_mouse.sort, taxa_are_rows = TRUE)
TAX_rough = tax_table(bact.mouse.taxa)
md.bact_mouse <- md.bact_mouse %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.bact_mouse)
exp2_bact_mouse_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

bact_mouse.sort2 <- bact_mouse.sort[,c(which(colSums(bact_mouse.sort)>=0))]
md.bact_mouse.sort2 <- md.bact_mouse[which(row.names(md.bact_mouse) %in% colnames(bact_mouse.sort2)),]

bact_mouse.nsampls <- apply(bact_mouse.sort2,1, gt0)
bact_mouse.clean <- bact_mouse.sort2[which(bact_mouse.nsampls>1),]

#Setup to create the phyloseq object
OTU = otu_table(bact_mouse.clean, taxa_are_rows = TRUE)
TAX = tax_table(bact.mouse.taxa)
samples = sample_data(md.bact_mouse.sort2)

#Convert!
exp2_bact_mouse <- phyloseq(OTU, TAX, samples)
exp2_bact_mouse_prev <- filter_taxa(exp2_bact_mouse , function(x) sum(x >= 1) > (0.05*length(x)), TRUE)


###Alpha Diversity###

fr_bact_mouse <- prune_taxa(taxa_sums(exp2_bact_mouse_rough) > 0, exp2_bact_mouse_rough)
richness_est_bact_mouse <- estimate_richness(fr_bact_mouse, measures = c("Simpson", "Shannon"))

richness_est_bact_mouse <- richness_est_bact_mouse %>%
  mutate(
    Organ = md.bact_mouse$Organ,
    Oxygen = md.bact_mouse$Oxygen,
    FMT = md.bact_mouse$FMT,
    Condition  = md.bact_mouse$Condition
  )
#Save

###Beta Diversity###

exp2_bact_mouse_rel_prev <- transform_sample_counts(exp2_bact_mouse_prev, function(x) x / sum(x) )

exp2_bact_mouse_otu_rel <- as.data.frame(t(exp2_bact_mouse_rel_prev@otu_table))
exp2_bact_mouse_tax_rel <- as.data.frame(exp2_bact_mouse_rel_prev@tax_table)
exp2_bact_mouse_meta_rel <- as.data.frame(exp2_bact_mouse_rel_prev@sam_data)

exp2_bact_mouse_rel_bray = vegdist(exp2_bact_mouse_otu_rel, method='bray')
exp2_bact_mouse_rel_pcoa <- ape::pcoa(exp2_bact_mouse_rel_bray)
exp2_bact_mouse_rel_pcoa$values

factor2_exp2 <- as.factor(exp2_bact_mouse_meta_rel$Condition)
type2_exp2 <- as.numeric(factor2_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Mouse lung 16s PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.5), c(-0.7, 0.3), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_bact_mouse_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type2_exp2], lwd = 1)
ordiellipse(exp2_bact_mouse_rel_pcoa$vectors[,1:2], factor2_exp2,kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_bact_mouse_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

permanova_pcoa_df <- data.frame(exp2_bact_mouse_meta_rel)
set.seed(1312)
permanova_bact_mouse_pcoa <- vegan::adonis2(exp2_bact_mouse_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_bact_mouse_pcoa <- permanova_pairwise(exp2_bact_mouse_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

bact_mouse_pcoa_dist <- vegdist(exp2_bact_mouse_otu_rel, method = "bray")
disp_bact_mouse_pcoa <- betadisper(bact_mouse_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_bact_mouse_pcoa <- permutest(disp_bact_mouse_pcoa, permutations = 10000)

print(permanova_bact_mouse_pcoa)
print(pairwise_permanova_bact_mouse_pcoa)
print(permdisp_bact_mouse_pcoa)
