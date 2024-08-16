library(tidyverse)
library(phyloseq)
library(plyr)
library(vegan)
library(ecole)

###Fig. 5B: LUNG MULTIKINGDOM###

###DATA PREP###

#Load data
fung_fluco_lung <- read_csv("Fluco_CPZ_Lung_ITS_ASV.csv")
fung.fluco.lung.taxa <- read_csv("Fluco_CPZ_Lung_ITS_Taxa.csv")

fung_fluco_lung <- fung_fluco_lung %>%
  tibble::column_to_rownames("Name")

fung.fluco.lung.taxa <- fung.fluco.lung.taxa %>%
  tibble::column_to_rownames("Name")

bact_fluco_lung <- read_csv("Fluco_CPZ_Lung_16S_ASV.csv")
bact.fluco.lung.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

bact_fluco_lung <- bact_fluco_lung %>%
  tibble::column_to_rownames("Name")

bact.fluco.lung.taxa <- bact.fluco.lung.taxa %>%
  tibble::column_to_rownames("Name")

roughmetadata_combined_fluco_lung <- read_csv("Fluco_Lung_Metadata.csv")

#Renaming the ITS ASVs so they don't overlap with 16S ASV names
row.names(fung_fluco_lung) <- paste("ASV",5386:5658,sep="")
bact_fluco_lung <-as.data.frame(t(bact_fluco_lung))

combined <- rbind.fill(bact_fluco_lung,fung_fluco_lung)
rownames(combined) <- c(row.names(bact_fluco_lung),row.names(fung_fluco_lung))
combined[is.na(combined)] <- 0

#Combine 16S and ITS taxonomy tables, rename the ASVs in the taxonomy table to match the new ITS names
row.names(fung.fluco.lung.taxa) <- paste("ASV",5386:5658,sep="")
combined.taxa <- rbind.fill(as.data.frame(bact.fluco.lung.taxa),as.data.frame(fung.fluco.lung.taxa))
rownames(combined.taxa) <- c(row.names(bact.fluco.lung.taxa),row.names(fung.fluco.lung.taxa))

combined.taxa <- as.matrix(combined.taxa)

combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined_fluco_lung[which(unlist(roughmetadata_combined_fluco_lung$Name) %in% colnames(combined.sort)),]

#Shows distribution of read depth
summary(colSums(combined.sort))

md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <50, 0 removed
combined.sort2 <- combined.sort[,c(which(colSums(combined.sort)>=50))]
combined.taxa2 <- combined.taxa[row.names(combined.taxa) %in% row.names(combined.sort2),]
md.combined.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.sort2)),]

OTU_rough = otu_table(combined.sort2, taxa_are_rows = TRUE)
TAX_rough = tax_table(combined.taxa)
samples_rough = sample_data(md.combined.sort2)

#Convert ITS data to phyloseq object 
fig5_lung_combined_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)


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
fig5_combined_fluco_lung <- phyloseq(OTU, TAX, samples)



###Fig 5B: LUNG SHANNON ALPHA DIVERSITY###

fr_combined_fluco_lung <- prune_taxa(taxa_sums(fig5_lung_combined_rough) > 0, fig5_lung_combined_rough)
richness_est_combined_fluco_lung <- estimate_richness(fr_combined_fluco_lung, measures = c("Shannon"))

wilcox_alpha_combined_fluco_lung <- t(sapply(richness_est_combined_fluco_lung, function(x) unlist(kruskal.test(x~sample_data(fr_combined_fluco_lung)$Sample.Name)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined_fluco_lung

richness_est_combined_fluco_lung <- richness_est_combined_fluco_lung %>%
  mutate(
    Organ = md.combined.fluco$Organ,
    Oxygen = md.combined.fluco$Oxygen,
    FMT = md.combined.fluco$FMT,
    Sample.Name  = md.combined.fluco$Sample.Name
  )

#Save data table and open it in Graphpad to create plot


###Fig 5B: LUNG BRAY-CURTIS BETA DIVERSITY PCOA###

fig5_combined_fluco_lung_rel <- transform_sample_counts(fig5_combined_fluco_lung, function(x) x / sum(x) )

fig5_combined_fluco_lung_otu_rel <- as.data.frame(t(fig5_combined_fluco_lung_rel@otu_table))
fig5_combined_fluco_lung_tax_rel <- as.data.frame(fig5_combined_fluco_lung_rel@tax_table)
fig5_combined_fluco_lung_meta_rel <- as.data.frame(fig5_combined_fluco_lung_rel@sam_data)

#Get Bray-Curtis dissimilarity matrix
fig5_combined_fluco_lung_bray = vegdist(fig5_combined_fluco_lung_otu_rel, method='bray')

fig5_combined_fluco_lung_pcoa <- ape::pcoa(fig5_combined_fluco_lung_bray)
fig5_combined_fluco_lung_pcoa$values

#This sets the color and factor order for the PCoA and PCA
factor_fig5 <- as.factor(fig5_combined_fluco_lung_meta_rel$Sample.Name)
factor_fig5 <- ordered(factor_fig5, levels = c("Fluco HO","SPF HO","Fluco NO", "SPF NO"))
type_fig5 <- as.numeric(factor_fig5)
pca_colors1 <- c("#94D57F","#B7B7B7")
pca_colors2 <- c("#000000","#000000","#94D57F","#000000")
ellipse_colors <- c("#94D57F","#000000","#94D57F","#000000")

#Plot PCoA
dpi=600
tiff("Fig 5B PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(fig5_combined_fluco_lung_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors1[type_fig5], lwd = 1)
points(fig5_combined_fluco_lung_pcoa$vectors[,1:2], pch = 1, cex = 1.3, col = pca_colors2[type_fig5], lwd = 1)
ordiellipse(fig5_combined_fluco_lung_pcoa$vectors[,1:2], factor_fig5,kind = "se", conf = 0.95, col = ellipse_colors)
ordispider(fig5_combined_fluco_lung_pcoa$vectors[,1:2], factor_fig5, label = TRUE)
dev.off()

#Significance test with PERMANOVA, pairwise PERMANOVA, and PERMDISP
permanova_pcoa_df <- data.frame(fig5_combined_fluco_lung_meta_rel)
set.seed(1312)
permanova_combined_fluco_lung_pcoa <- vegan::adonis2(fig5_combined_fluco_lung_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_fluco_lung_pcoa <- permanova_pairwise(fig5_combined_fluco_lung_otu_rel, grp = permanova_pcoa_df$Sample.Name, permutations = 10000, method = "bray", padj = "fdr")

combined_fluco_lung_pcoa_dist <- vegdist(fig5_combined_fluco_lung_otu_rel, method = "bray")
disp_combined_fluco_lung_pcoa <- betadisper(combined_fluco_lung_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_combined_fluco_lung_pcoa <- permutest(disp_combined_fluco_lung_pcoa, permutations = 10000)

print(permanova_combined_fluco_lung_pcoa)
print(pairwise_permanova_combined_fluco_lung_pcoa)
print(permdisp_combined_fluco_lung_pcoa)


###Fig 5B: LUNG BETA DIVERSITY PCA###

fig5_combined_fluco_lung_rel <- transform_sample_counts(fig5_combined_fluco_lung, function(x) x / sum(x) )

fig5_combined_fluco_lung_otu <- as.data.frame(t(fig5_combined_fluco_lung_rel@otu_table))
fig5_combined_fluco_lung_tax <- as.data.frame(fig5_combined_fluco_lung_rel@tax_table)
fig5_combined_fluco_lung_meta <- as.data.frame(fig5_combined_fluco_lung_rel@sam_data)

#Hellinger transform and ordinate
fig5_combined_fluco_lung_otu_hel <- decostand(fig5_combined_fluco_lung_otu, "hellinger")
fig5_combined_fluco_lung_otu_pca <- rda(fig5_combined_fluco_lung_otu_hel)
summary(fig5_combined_fluco_lung_otu_pca)$cont


#Sets which taxa labels should be shown in the loading plot and which should be hidden, to prevent overcrowding
priority_fig5_combined_fluco_lung <- colSums(fig5_combined_fluco_lung_otu_hel)
labels_fig5_combined_fluco_lung <- orditorp(fig5_combined_fluco_lung_otu_pca, "sp", label = fig5_combined_fluco_lung_tax$Genus, priority=priority_fig5_combined_fluco_lung)

#PCA Loading Plot
dpi = 600
tiff("Fig 5B PCA Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(fig5_combined_fluco_lung_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (18.46% Explained)", ylab="PC2 (11.20% Explained)", ylim = c(-0.5,0.4),xlim = c(-0.6,0.6))
orditorp(fig5_combined_fluco_lung_otu_pca, "sp", label = fig5_combined_fluco_lung_tax$Genus, priority=priority_fig5_combined_fluco_lung, select = (labels_fig5_combined_fluco_lung == TRUE), cex = 0.7)
dev.off()

#PCA sample plot
#Color and factor mappings are the same as the PCoA
tiff("Fig 5B PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.9,0.8),c(-0.9,0.8), font = 2, font.lab = 2, xlab="PC1 (18.46% Explained)", ylab="PC2 (11.20% Explained)", type="n")
points(fig5_combined_fluco_lung_otu_pca, pch = 21, cex = 1.3, bg = pca_colors1[factor_fig5], lwd = 1)
points(fig5_combined_fluco_lung_otu_pca, pch = 1, cex = 1.3, col = pca_colors2[factor_fig5], lwd = 1)
ordiellipse(fig5_combined_fluco_lung_otu_pca, factor_fig5, col = ellipse_colors)
ordispider(fig5_combined_fluco_lung_otu_pca, factor_fig5, label = TRUE)
dev.off()

#PERMANOVA, pairwise PERMANOVA, and PERMDISP based on Euclidean distances
permanova_fluco_lung_df <- data.frame(fig5_combined_fluco_lung_meta)
set.seed(1312)
permanova_fluco_lung <- vegan::adonis2(fig5_combined_fluco_lung_otu_hel ~ Sample.Name, data = permanova_fluco_lung_df, method="euclidean", permutations = 10000)
pairwise_permanova_fluco_lung <- permanova_pairwise(fig5_combined_fluco_lung_otu_hel, grp = permanova_fluco_lung_df$Sample.Name, permutations = 10000, method = "euclidean", padj = "fdr")

fluco_lung_hel_dist <- vegdist(fig5_combined_fluco_lung_otu_hel, method = "euclidean")
disp_fluco_lung <- betadisper(fluco_lung_hel_dist, permanova_fluco_lung_df$Sample.Name)
set.seed(1312)
permdisp_fluco_lung <- permutest(disp_fluco_lung, permutations = 10000)

print(permanova_fluco_lung)
print(pairwise_permanova_fluco_lung)
print(permdisp_fluco_lung)


###Fig 5C: GUT MULTIKINGDOM###

###DATA PREP###

#Load data
fung_fluco_gut <- read_csv("Fluco_CPZ_Gut_ITS_ASV.csv")
fung.fluco.gut.taxa <- read_csv("Fluco_CPZ_Gut_ITS_Taxa.csv")

fung_fluco_gut <- fung_fluco_gut %>%
  tibble::column_to_rownames("Name")

fung.fluco.gut.taxa <- fung.fluco.gut.taxa %>%
  tibble::column_to_rownames("Name")

bact_fluco_gut <- read_csv("Fluco_CPZ_Gut_16S_ASV.csv")
bact.fluco.gut.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

bact_fluco_gut <- bact_fluco_gut %>%
  tibble::column_to_rownames("Name")

bact.fluco.gut.taxa <- bact.fluco.gut.taxa %>%
  tibble::column_to_rownames("Name")

roughmetadata_combined_fluco_gut <- read_csv("Fluco_Gut_Metadata.csv")

#Renaming the ITS ASVs so they don't overlap with 16S ASV names
row.names(fung_fluco_gut) <- paste("ASV",5386:5434,sep="")
bact_fluco_gut <-as.data.frame(t(bact_fluco_gut))

combined <- rbind.fill(bact_fluco_gut,fung_fluco_gut)
rownames(combined) <- c(row.names(bact_fluco_gut),row.names(fung_fluco_gut))
combined[is.na(combined)] <- 0

#Combine 16S and ITS taxonomy tables, rename the ASVs in the taxonomy table to match the new ITS names
row.names(fung.fluco.gut.taxa) <- paste("ASV",5386:5434,sep="")
combined.taxa <- rbind.fill(as.data.frame(bact.fluco.gut.taxa),as.data.frame(fung.fluco.gut.taxa))
rownames(combined.taxa) <- c(row.names(bact.fluco.gut.taxa),row.names(fung.fluco.gut.taxa))

combined.taxa <- as.matrix(combined.taxa)

combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined_fluco_gut[which(unlist(roughmetadata_combined_fluco_gut$Name) %in% colnames(combined.sort)),]

#Shows distribution of read depth
summary(colSums(combined.sort))

md.combined <- md.combined %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <50, 0 removed
combined.sort2 <- combined.sort[,c(which(colSums(combined.sort)>=50))]
combined.taxa2 <- combined.taxa[row.names(combined.taxa) %in% row.names(combined.sort2),]
md.combined.sort2 <- md.combined[which(row.names(md.combined) %in% colnames(combined.sort2)),]

OTU_rough = otu_table(combined.sort2, taxa_are_rows = TRUE)
TAX_rough = tax_table(combined.taxa)
samples_rough = sample_data(md.combined.sort2)

#Convert ITS data to phyloseq object 
fig5_gut_combined_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)


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
fig5_combined_fluco_gut <- phyloseq(OTU, TAX, samples)



###Fig 5C: GUT SHANNON ALPHA DIVERSITY###

fr_combined_fluco_gut <- prune_taxa(taxa_sums(fig5_gut_combined_rough) > 0, fig5_gut_combined_rough)
richness_est_combined_fluco_gut <- estimate_richness(fr_combined_fluco_gut, measures = c("Shannon"))

wilcox_alpha_combined_fluco_gut <- t(sapply(richness_est_combined_fluco_gut, function(x) unlist(kruskal.test(x~sample_data(fr_combined_fluco_gut)$Sample.Name)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined_fluco_gut

richness_est_combined_fluco_gut <- richness_est_combined_fluco_gut %>%
  mutate(
    Organ = md.combined.sort2$Organ,
    Oxygen = md.combined.sort2$Oxygen,
    FMT = md.combined.sort2$FMT,
    Sample.Name  = md.combined.sort2$Sample.Name
  )

#Save data table and open it in Graphpad to create plot


###Fig 5C: GUT BRAY-CURTIS BETA DIVERSITY PCOA###

fig5_combined_fluco_gut_rel <- transform_sample_counts(fig5_combined_fluco_gut, function(x) x / sum(x) )

fig5_combined_fluco_gut_otu_rel <- as.data.frame(t(fig5_combined_fluco_gut_rel@otu_table))
fig5_combined_fluco_gut_tax_rel <- as.data.frame(fig5_combined_fluco_gut_rel@tax_table)
fig5_combined_fluco_gut_meta_rel <- as.data.frame(fig5_combined_fluco_gut_rel@sam_data)

#Get Bray-Curtis dissimilarity matrix
fig5_combined_fluco_gut_bray = vegdist(fig5_combined_fluco_gut_otu_rel, method='bray')

fig5_combined_fluco_gut_pcoa <- ape::pcoa(fig5_combined_fluco_gut_bray)
fig5_combined_fluco_gut_pcoa$values

#This sets the color and factor order for the PCoA and PCA
factor_fig5 <- as.factor(fig5_combined_fluco_gut_meta_rel$Sample.Name)
factor_fig5 <- ordered(factor_fig5, levels = c("Fluco HO","SPF HO","Fluco NO", "SPF NO"))
type_fig5 <- as.numeric(factor_fig5)
pca_colors1 <- c("#94D57F","#B7B7B7")
pca_colors2 <- c("#000000","#000000","#94D57F","#000000")
ellipse_colors <- c("#94D57F","#000000","#94D57F","#000000")

#Plot PCoA
dpi=600
tiff("Fig 5C PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.3), c(-0.3, 0.2), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(fig5_combined_fluco_gut_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors1[type_fig5], lwd = 1)
points(fig5_combined_fluco_gut_pcoa$vectors[,1:2], pch = 1, cex = 1.3, col = pca_colors2[type_fig5], lwd = 1)
ordiellipse(fig5_combined_fluco_gut_pcoa$vectors[,1:2], factor_fig5,kind = "se", conf = 0.95, col = ellipse_colors)
ordispider(fig5_combined_fluco_gut_pcoa$vectors[,1:2], factor_fig5, label = TRUE)
dev.off()

#Significance test with PERMANOVA, pairwise PERMANOVA, and PERMDISP
permanova_pcoa_df <- data.frame(fig5_combined_fluco_gut_meta_rel)
set.seed(1312)
permanova_combined_fluco_gut_pcoa <- vegan::adonis2(fig5_combined_fluco_gut_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_fluco_gut_pcoa <- permanova_pairwise(fig5_combined_fluco_gut_otu_rel, grp = permanova_pcoa_df$Sample.Name, permutations = 10000, method = "bray", padj = "fdr")

combined_fluco_gut_pcoa_dist <- vegdist(fig5_combined_fluco_gut_otu_rel, method = "bray")
disp_combined_fluco_gut_pcoa <- betadisper(combined_fluco_gut_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_combined_fluco_gut_pcoa <- permutest(disp_combined_fluco_gut_pcoa, permutations = 10000)

print(permanova_combined_fluco_gut_pcoa)
print(pairwise_permanova_combined_fluco_gut_pcoa)
print(permdisp_combined_fluco_gut_pcoa)


###Fig 5C: GUT BETA DIVERSITY PCA###

fig5_combined_fluco_gut_rel <- transform_sample_counts(fig5_combined_fluco_gut, function(x) x / sum(x) )

fig5_combined_fluco_gut_otu <- as.data.frame(t(fig5_combined_fluco_gut_rel@otu_table))
fig5_combined_fluco_gut_tax <- as.data.frame(fig5_combined_fluco_gut_rel@tax_table)
fig5_combined_fluco_gut_meta <- as.data.frame(fig5_combined_fluco_gut_rel@sam_data)

#Hellinger transform and ordinate
fig5_combined_fluco_gut_otu_hel <- decostand(fig5_combined_fluco_gut_otu, "hellinger")
fig5_combined_fluco_gut_otu_pca <- rda(fig5_combined_fluco_gut_otu_hel)
summary(fig5_combined_fluco_gut_otu_pca)$cont


#Sets which taxa labels should be shown in the loading plot and which should be hidden, to prevent overcrowding
priority_fig5_combined_fluco_gut <- colSums(fig5_combined_fluco_gut_otu_hel)
labels_fig5_combined_fluco_gut <- orditorp(fig5_combined_fluco_gut_otu_pca, "sp", label = fig5_combined_fluco_gut_tax$Genus, priority=priority_fig5_combined_fluco_gut)

#PCA Loading Plot
dpi = 600
tiff("Fig 5C PCA Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(fig5_combined_fluco_gut_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (70.77% Explained)", ylab="PC2 (11.23% Explained)", ylim = c(-0.8,0.6),xlim = c(-0.6,0.5))
orditorp(fig5_combined_fluco_gut_otu_pca, "sp", label = fig5_combined_fluco_gut_tax$Genus, priority=priority_fig5_combined_fluco_gut, select = (labels_fig5_combined_fluco_gut == TRUE), cex = 0.7)
dev.off()

#PCA sample plot
#Color and factor mappings are the same as the PCoA
tiff("Fig 5C PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6,0.4),c(-0.5,0.8), font = 2, font.lab = 2, xlab="PC1 (70.77% Explained)", ylab="PC2 (11.23% Explained)", type="n")
points(fig5_combined_fluco_gut_otu_pca, pch = 21, cex = 1.3, bg = pca_colors1[factor_fig5], lwd = 1)
points(fig5_combined_fluco_gut_otu_pca, pch = 1, cex = 1.3, col = pca_colors2[factor_fig5], lwd = 1)
ordiellipse(fig5_combined_fluco_gut_otu_pca, factor_fig5, col = ellipse_colors)
ordispider(fig5_combined_fluco_gut_otu_pca, factor_fig5, label = TRUE)
dev.off()

#PERMANOVA, pairwise PERMANOVA, and PERMDISP based on Euclidean distances
permanova_fluco_gut_df <- data.frame(fig5_combined_fluco_gut_meta)
set.seed(1312)
permanova_fluco_gut <- vegan::adonis2(fig5_combined_fluco_gut_otu_hel ~ Sample.Name, data = permanova_fluco_gut_df, method="euclidean", permutations = 10000)
pairwise_permanova_fluco_gut <- permanova_pairwise(fig5_combined_fluco_gut_otu_hel, grp = permanova_fluco_gut_df$Sample.Name, permutations = 10000, method = "euclidean", padj = "fdr")

fluco_gut_hel_dist <- vegdist(fig5_combined_fluco_gut_otu_hel, method = "euclidean")
disp_fluco_gut <- betadisper(fluco_gut_hel_dist, permanova_fluco_gut_df$Sample.Name)
set.seed(1312)
permdisp_fluco_gut <- permutest(disp_fluco_gut, permutations = 10000)

print(permanova_fluco_gut)
print(pairwise_permanova_fluco_gut)
print(permdisp_fluco_gut)
