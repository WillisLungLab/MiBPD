library(tidyverse)
library(phyloseq)
library(vegan)
library(ecole)

###FIG. 4B: LUNG MULTIKINGDOM###

###DATA PREP###

#Load data
fung_mouse_lung <- read_csv("Mouse_Lung_ITS_ASVs.csv")
fung.mouse.lung.taxa <- read_csv("Mouse_Lung_ITS_Taxa.csv")

fung_mouse_lung <- fung_mouse_lung %>%
  tibble::column_to_rownames("Name")

fung.mouse.lung.taxa <- fung.mouse.lung.taxa %>%
  tibble::column_to_rownames("Name")

bact_mouse_lung <- read_csv("Mouse_Lung_16s_ASVs.csv")
bact.mouse.lung.taxa <- read_csv("Mouse_Lung_16s_Taxa.csv")

bact_mouse_lung <- bact_mouse_lung %>%
  tibble::column_to_rownames("Name")

bact.mouse.lung.taxa <- bact.mouse.lung.taxa %>%
  tibble::column_to_rownames("Name")

roughmetadata_combined_mouse_lung <- read_csv("Mouse_Lung_Metadata.csv")


#Combine 16S and ITS dataframes and fill in empty sections as needed
combined_mouse_lung <- rbind.fill(bact_mouse_lung,fung_mouse_lung)
rownames(combined_mouse_lung) <- c(row.names(bact_mouse_lung),row.names(fung_mouse_lung))
combined_mouse_lung[is.na(combined_mouse_lung)] <- 0

#Combine 16S and ITS taxonomy tables
combined.mouse.lung.taxa <- rbind.fill(as.data.frame(bact.mouse.lung.taxa),as.data.frame(fung.mouse.lung.taxa))
rownames(combined.mouse.lung.taxa) <- c(row.names(bact.mouse.lung.taxa),row.names(fung.mouse.lung.taxa))

combined.mouse.lung.taxa <- as.matrix(combined.mouse.lung.taxa)

#Load metadata and the rest of this is all the same

combined.mouse.lung.sort <- combined_mouse_lung[,order(colnames(combined_mouse_lung))]
md.combined <- roughmetadata_combined_mouse_lung[which(unlist(roughmetadata_combined_mouse_lung$Name) %in% colnames(combined.mouse.lung.sort)),]
summary(colSums(combined.mouse.lung.sort))

md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")

#Filters out samples with no reads
combined.mouse.lung.sort2 <- combined.mouse.lung.sort[,c(which(colSums(combined.mouse.lung.sort)>=0))]
md.combined.mouse.lung.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.mouse.lung.sort2)),]

#Phyloseq object for alpha diversity
OTU_rough = otu_table(combined.mouse.lung.sort2, taxa_are_rows = TRUE)
TAX_rough = tax_table(combined.mouse.lung.taxa)
samples_rough = sample_data(md.combined.mouse.lung.sort2)

#Get phyloseq object for alpha diversity
exp2_combined_mouse_lung_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

#Remove singletons (function written by Dr. Laura Tipton)
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

combined.mouse.lung.nsampls <- apply(combined.mouse.lung.sort2, 1, gt0)
combined.mouse.lung.clean <- combined.mouse.lung.sort2[which(combined.mouse.lung.nsampls>1),]

#Phyloseq object for beta diversity and supplemental analyses
OTU = otu_table(combined.mouse.lung.clean, taxa_are_rows = TRUE)
TAX = tax_table(combined.mouse.lung.taxa)
samples = sample_data(md.combined.mouse.lung.sort2)

exp2_combined_mouse_lung <- phyloseq(OTU, TAX, samples)


###FIG 4B: LUNG SHANNON ALPHA DIVERSITY###

fr_combined_mouse_lung <- prune_taxa(taxa_sums(exp2_combined_mouse_lung_rough) > 0, exp2_combined_mouse_lung_rough)
richness_est_combined_mouse_lung <- estimate_richness(fr_combined_mouse_lung, measures = c("Shannon"))

wilcox_alpha_combined_mouse_lung <- t(sapply(richness_est_combined_mouse_lung, function(x) unlist(kruskal.test(x~sample_data(fr_combined_mouse_lung)$Condition)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined_mouse_lung

richness_est_combined_mouse_lung <- richness_est_combined_mouse_lung %>%
  mutate(
    Organ = md.combined.mouse$Organ,
    Oxygen = md.combined.mouse$Oxygen,
    FMT = md.combined.mouse$FMT,
    Condition  = md.combined.mouse$Condition
  )

#Save data table and open it in Graphpad to create plot


###FIG 4B: LUNG BRAY-CURTIS BETA DIVERSITY PCOA###

exp2_combined_mouse_lung_rel <- transform_sample_counts(exp2_combined_mouse_lung, function(x) x / sum(x) )

exp2_combined_mouse_lung_otu_rel <- as.data.frame(t(exp2_combined_mouse_lung_rel@otu_table))
exp2_combined_mouse_lung_tax_rel <- as.data.frame(exp2_combined_mouse_lung_rel@tax_table)
exp2_combined_mouse_lung_meta_rel <- as.data.frame(exp2_combined_mouse_lung_rel@sam_data)

#Get Bray-Curtis dissimilarity matrix
exp2_combined_mouse_lung_rel_bray = vegdist(exp2_combined_mouse_lung_otu_rel, method='bray')

exp2_combined_mouse_lung_rel_pcoa <- ape::pcoa(exp2_combined_mouse_lung_rel_bray)
exp2_combined_mouse_lung_rel_pcoa$values

#This sets the color and factor order for the PCoA and PCA
factor_exp2 <- as.factor(exp2_combined_mouse_lung_meta_rel$Condition)
factor_exp2 <- ordered(factor_exp2, levels = c("BPD_HO","No_BPD_HO","BPD_NO", "No_BPD_NO"))
type_exp2 <- as.numeric(factor2_exp2)
pca_colors1 <- c("#007016","#94D57F","#FFFFFF","#FFFFFF")
pca_colors2 <- c("#000000","#000000","#007016","#94D57F")
ellipse_colors <- c("#007016","#94D57F","#007016","#94D57F")

#Plot PCoA
dpi=600
tiff("Fig 4B PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_combined_mouse_lung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors1[type_exp2], lwd = 1)
points(exp2_combined_mouse_lung_rel_pcoa$vectors[,1:2], pch = 1, cex = 1.3, col = pca_colors2[type_exp2], lwd = 1)
ordiellipse(exp2_combined_mouse_lung_rel_pcoa$vectors[,1:2], factor2_exp2,kind = "se", conf = 0.95, col = ellipse_colors)
ordispider(exp2_combined_mouse_lung_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

#Significance test with PERMANOVA, pairwise PERMANOVA, and PERMDISP
permanova_pcoa_df <- data.frame(exp2_combined_mouse_lung_meta_rel)
set.seed(1312)
permanova_combined_mouse_lung_pcoa <- vegan::adonis2(exp2_combined_mouse_lung_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_mouse_lung_pcoa <- permanova_pairwise(exp2_combined_mouse_lung_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

combined_mouse_lung_pcoa_dist <- vegdist(exp2_combined_mouse_lung_otu_rel, method = "bray")
disp_combined_mouse_lung_pcoa <- betadisper(combined_mouse_lung_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_combined_mouse_lung_pcoa <- permutest(disp_combined_mouse_lung_pcoa, permutations = 10000)

print(permanova_combined_mouse_lung_pcoa)
print(pairwise_permanova_combined_mouse_lung_pcoa)
print(permdisp_combined_mouse_lung_pcoa)


###FIG 4B: LUNG BETA DIVERSITY PCA###

exp2_combined_mouse_lung_rel <- transform_sample_counts(exp2_combined_mouse_lung, function(x) x / sum(x) )

exp2_combined_mouse_lung_otu <- as.data.frame(t(exp2_combined_mouse_lung_rel@otu_table))
exp2_combined_mouse_lung_tax <- as.data.frame(exp2_combined_mouse_lung_rel@tax_table)
exp2_combined_mouse_lung_meta <- as.data.frame(exp2_combined_mouse_lung_rel@sam_data)

#Hellinger transform and ordinate
exp2_combined_mouse_lung_otu_hel <- decostand(exp2_combined_mouse_lung_otu, "hellinger")
exp2_combined_mouse_lung_otu_pca <- rda(exp2_combined_mouse_lung_otu_hel)
summary(exp2_combined_mouse_lung_otu_pca)$cont


#Sets which taxa labels should be shown in the loading plot and which should be hidden, to prevent overcrowding
priority_exp2_combined_mouse_lung <- colSums(exp2_combined_mouse_lung_otu_hel)
labels_exp2_combined_mouse_lung <- orditorp(exp2_combined_mouse_lung_otu_pca, "sp", label = exp2_combined_mouse_lung_tax$Genus, priority=priority_exp2_combined_mouse_lung)

#PCA Loading Plot
dpi = 600
tiff("Fig 4B PCA Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp2_combined_mouse_lung_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (14.38% Explained)", ylab="PC2 (12.06% Explained)", ylim = c(-0.4,0.5),xlim = c(-0.6,0.5))
orditorp(exp2_combined_mouse_lung_otu_pca, "sp", label = exp2_combined_mouse_lung_tax$Genus, priority=priority_exp2_combined_mouse_lung, select = (labels_exp2_combined_mouse_lung == TRUE), cex = 0.7)
dev.off()

#PCA sample plot
#Color and factor mappings are the same as the PCoA
tiff("Fig 4B PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.8,0.8),c(-0.9,0.8), font = 2, font.lab = 2, xlab="PC1 (14.38% Explained)", ylab="PC2 (12.06% Explained)", type="n")
points(exp2_combined_mouse_lung_otu_pca, pch = 21, cex = 1.3, bg = pca_colors1[factor_exp2], lwd = 1)
points(exp2_combined_mouse_lung_otu_pca, pch = 1, cex = 1.3, col = pca_colors2[factor_exp2], lwd = 1)
ordiellipse(exp2_combined_mouse_lung_otu_pca, factor_exp2, col = ellipse_colors)
ordispider(exp2_combined_mouse_lung_otu_pca, factor_exp2, label = TRUE)
dev.off()

#PERMANOVA, pairwise PERMANOVA, and PERMDISP based on Euclidean distances
permanova_mouse_lung_df <- data.frame(exp2_combined_mouse_lung_meta)
set.seed(1312)
permanova_mouse_lung <- vegan::adonis2(exp2_combined_mouse_lung_otu_hel ~ Condition, data = permanova_mouse_lung_df, method="euclidean", permutations = 10000)
pairwise_permanova_mouse_lung <- permanova_pairwise(exp2_combined_mouse_lung_otu_hel, grp = permanova_mouse_lung_df$Condition, permutations = 10000, method = "euclidean", padj = "fdr")

mouse_lung_hel_dist <- vegdist(exp2_combined_mouse_lung_otu_hel, method = "euclidean")
disp_mouse_lung <- betadisper(mouse_lung_hel_dist, permanova_mouse_lung_df$Condition)
set.seed(1312)
permdisp_mouse_lung <- permutest(disp_mouse_lung, permutations = 10000)

print(permanova_mouse_lung)
print(pairwise_permanova_mouse_lung)
print(permdisp_mouse_lung)


###FIG 4C: GUT MULTIKINGDOM DATA PREP###

#Load data
fung_mouse_gut <- read_csv("Mouse_Gut_ITS_ASVs.csv")
fung.mouse.gut.taxa <- read_csv("Mouse_Gut_ITS_Taxa.csv")

fung_mouse_gut <- fung_mouse_gut %>%
  tibble::column_to_rownames("Name")

fung.mouse.gut.taxa <- fung.mouse.gut.taxa %>%
  tibble::column_to_rownames("Name")

bact_mouse_gut <- read_csv("Mouse_Gut_16s_ASVs.csv")
bact.mouse.gut.taxa <- read_csv("Mouse_Gut_16s_Taxa_noNA.csv")

bact_mouse_gut <- bact_mouse_gut %>%
  tibble::column_to_rownames("Name")

bact.mouse.gut.taxa <- bact.mouse.gut.taxa %>%
  tibble::column_to_rownames("Name")

roughmetadata_combined_mouse_gut <- read_csv("Mouse_Gut_Metadata.csv")


#Combine 16S and ITS dataframes and fill in empty sections as needed
combined_mouse_gut <- rbind.fill(bact_mouse_gut,fung_mouse_gut)
rownames(combined_mouse_gut) <- c(row.names(bact_mouse_gut),row.names(fung_mouse_gut))
combined_mouse_gut[is.na(combined_mouse_gut)] <- 0

#Combine 16S and ITS taxonomy tables
combined.mouse.gut.taxa <- rbind.fill(as.data.frame(bact.mouse.gut.taxa),as.data.frame(fung.mouse.gut.taxa))
rownames(combined.mouse.gut.taxa) <- c(row.names(bact.mouse.gut.taxa),row.names(fung.mouse.gut.taxa))

combined.mouse.gut.taxa <- as.matrix(combined.mouse.gut.taxa)

#Load metadata and the rest of this is all the same
combined.mouse.gut.sort <- combined_mouse_gut[,order(colnames(combined_mouse_gut))]
md.combined <- roughmetadata_combined_mouse_gut[which(unlist(roughmetadata_combined_mouse_gut$Name) %in% colnames(combined.mouse.gut.sort)),]
summary(colSums(combined.mouse.gut.sort))

md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth of 0
combined.mouse.gut.sort2 <- combined.mouse.gut.sort[,c(which(colSums(combined.mouse.gut.sort)>=0))]
md.combined.mouse.gut.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.mouse.gut.sort2)),]

OTU_rough = otu_table(combined.mouse.gut.sort2, taxa_are_rows = TRUE)
TAX_rough = tax_table(combined.mouse.gut.taxa)
samples_rough = sample_data(md.combined.mouse.gut.sort2)

#Get phyloseq object for alpha diversity
exp2_combined_mouse_gut_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

combined.mouse.gut.nsampls <- apply(combined.mouse.gut.sort2, 1, gt0)
combined.mouse.gut.clean <- combined.mouse.gut.sort2[which(combined.mouse.gut.nsampls>1),]

#Setup to create the main phyloseq object
OTU = otu_table(combined.mouse.gut.clean, taxa_are_rows = TRUE)
TAX = tax_table(combined.mouse.gut.taxa)
samples = sample_data(md.combined.mouse.gut.sort2)

#Convert!
exp2_combined_mouse_gut <- phyloseq(OTU, TAX, samples)


###FIG 4C: GUT SHANNON ALPHA DIVERSITY###

fr_combined_mouse_gut <- prune_taxa(taxa_sums(exp2_combined_mouse_gut_rough) > 0, exp2_combined_mouse_gut_rough)
richness_est_combined_mouse_gut <- estimate_richness(fr_combined_mouse_gut, measures = c("Shannon"))

wilcox_alpha_combined_mouse_gut <- t(sapply(richness_est_combined_mouse_gut, function(x) unlist(kruskal.test(x~sample_data(fr_combined_mouse_gut)$Condition)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined_mouse_gut

richness_est_combined_mouse_gut <- richness_est_combined_mouse_gut %>%
  mutate(
    Organ = md.combined.mouse$Organ,
    Oxygen = md.combined.mouse$Oxygen,
    FMT = md.combined.mouse$FMT,
    Condition  = md.combined.mouse$Condition
  )

#Save and open in GraphPad


###FIG 4C: GUT BRAY-CURTIS BETA DIVERSITY PCOA###

exp2_combined_mouse_gut_rel <- transform_sample_counts(exp2_combined_mouse_gut, function(x) x / sum(x) )

exp2_combined_mouse_gut_otu_rel <- as.data.frame(t(exp2_combined_mouse_gut_rel@otu_table))
exp2_combined_mouse_gut_tax_rel <- as.data.frame(exp2_combined_mouse_gut_rel@tax_table)
exp2_combined_mouse_gut_meta_rel <- as.data.frame(exp2_combined_mouse_gut_rel@sam_data)

#Create Bray-Curtis dissimilarity matrix
exp2_combined_mouse_gut_rel_bray = vegdist(exp2_combined_mouse_gut_otu_rel, method='bray')
exp2_combined_mouse_gut_rel_pcoa <- ape::pcoa(exp2_combined_mouse_gut_rel_bray)
exp2_combined_mouse_gut_rel_pcoa$values

#This sets the color and factor order for the PCoA and PCA
factor_exp2 <- as.factor(exp2_combined_mouse_gut_meta_rel$Condition)
factor_exp2 <- ordered(factor_exp2, levels = c("BPD_HO","NoBPD_HO","BPD_NO", "NoBPD_NO"))
type_exp2 <- as.numeric(factor_exp2)
pca_colors1 <- c("#007016","#94D57F","#FFFFFF","#FFFFFF")
pca_colors2 <- c("#000000","#000000","#007016","#94D57F")
ellipse_colors <- c("#007016","#94D57F","#007016","#94D57F")

#Plot PCoA
dpi=600
tiff("Fig 4C PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.7, 0.4), c(-0.4, 0.7), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_combined_mouse_gut_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors1[type_exp2], lwd = 1)
points(exp2_combined_mouse_gut_rel_pcoa$vectors[,1:2], pch = 1, cex = 1.3, col = pca_colors2[type_exp2], lwd = 1)
ordiellipse(exp2_combined_mouse_gut_rel_pcoa$vectors[,1:2], factor2_exp2,kind = "se", conf = 0.95, col = ellipse_colors)
ordispider(exp2_combined_mouse_gut_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

#PERMANOVA, pairwise PERMANOVA, and PERMDISP for Bray-Curtis distance
permanova_pcoa_df <- data.frame(exp2_combined_mouse_gut_meta_rel)
set.seed(1312)
permanova_combined_mouse_gut_pcoa <- vegan::adonis2(exp2_combined_mouse_gut_otu_rel ~ Condition, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_mouse_gut_pcoa <- permanova_pairwise(exp2_combined_mouse_gut_otu_rel, grp = permanova_pcoa_df$Condition, permutations = 10000, method = "bray", padj = "fdr")

combined_mouse_gut_pcoa_dist <- vegdist(exp2_combined_mouse_gut_otu_rel, method = "bray")
disp_combined_mouse_gut_pcoa <- betadisper(combined_mouse_gut_pcoa_dist, permanova_pcoa_df$Condition)
set.seed(1312)
permdisp_combined_mouse_gut_pcoa <- permutest(disp_combined_mouse_gut_pcoa, permutations = 10000)

print(permanova_combined_mouse_gut_pcoa)
print(pairwise_permanova_combined_mouse_gut_pcoa)
print(permdisp_combined_mouse_gut_pcoa)


###FIG 4C: GUT BETA DIVERSITY PCA###

exp2_combined_mouse_gut_rel <- transform_sample_counts(exp2_combined_mouse_gut, function(x) x / sum(x) )

exp2_combined_mouse_gut_otu <- as.data.frame(t(exp2_combined_mouse_gut_rel@otu_table))
exp2_combined_mouse_gut_tax <- as.data.frame(exp2_combined_mouse_gut_rel@tax_table)
exp2_combined_mouse_gut_meta <- as.data.frame(exp2_combined_mouse_gut_rel@sam_data)

#Hellinger transform and ordinate
exp2_combined_mouse_gut_otu_hel <- decostand(exp2_combined_mouse_gut_otu, "hellinger")
exp2_combined_mouse_gut_otu_pca <- rda(exp2_combined_mouse_gut_otu_hel)
summary(exp2_combined_mouse_gut_otu_pca)$cont


#Sets taxa label priority for loading plot
priority_exp2_combined_mouse_gut <- colSums(exp2_combined_mouse_gut_otu_hel)
labels_exp2_combined_mouse_gut <- orditorp(exp2_combined_mouse_gut_otu_pca, "sp", label = exp2_combined_mouse_gut_tax$Genus, priority=priority_exp2_combined_mouse_gut)
dpi = 600

#PCA loading plot
tiff("Fig 4C PCA Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp2_combined_mouse_gut_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (38.54% Explained)", ylab="PC2 (21.89% Explained)", ylim = c(-0.5,0.8),xlim = c(-1.1,1.2))
orditorp(exp2_combined_mouse_gut_otu_pca, "sp", label = exp2_combined_mouse_gut_tax$Genus, priority=priority_exp2_combined_mouse_gut, select = (labels_exp2_combined_mouse_gut == TRUE), cex = 0.7)
dev.off()

#PCA sample plot
#Color and factor mappings are the same as the PCoA
tiff("Fig 4C PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5,0.8),c(-0.4,1.0), font = 2, font.lab = 2, xlab="PC1 (38.54% Explained)", ylab="PC2 (21.89% Explained)", type="n")
points(exp2_combined_mouse_gut_otu_pca, pch = 21, cex = 1.3, bg = pca_colors1[factor_exp2], lwd = 1)
points(exp2_combined_mouse_gut_otu_pca, pch = 1, cex = 1.3, col = pca_colors2[factor_exp2], lwd = 1)
ordiellipse(exp2_combined_mouse_gut_otu_pca, factor_exp2, col = ellipse_colors)
ordispider(exp2_combined_mouse_gut_otu_pca, factor_exp2, label = TRUE)
dev.off()

####PERMANOVA, pairwise PERMANOVA, and PERMDISP based on Euclidean distances####
permanova_mouse_gut_df <- data.frame(exp2_combined_mouse_gut_meta)
set.seed(1312)
permanova_mouse_gut <- vegan::adonis2(exp2_combined_mouse_gut_otu_hel ~ Condition, data = permanova_mouse_gut_df, method="euclidean", permutations = 10000)
pairwise_permanova_mouse_gut <- permanova_pairwise(exp2_combined_mouse_gut_otu_hel, grp = permanova_mouse_gut_df$Condition, permutations = 10000, method = "euclidean", padj = "fdr")

mouse_gut_hel_dist <- vegdist(exp2_combined_mouse_gut_otu_hel, method = "euclidean")
disp_mouse_gut <- betadisper(mouse_gut_hel_dist, permanova_mouse_gut_df$Condition)
set.seed(1312)
permdisp_mouse_gut <- permutest(disp_mouse_gut, permutations = 10000)

print(permanova_mouse_gut)
print(pairwise_permanova_mouse_gut)
print(permdisp_mouse_gut)