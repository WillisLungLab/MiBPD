library(tidyverse)
library(phyloseq)
library(plyr)
library(vegan)
library(ecole)

###Fig. 6B: LUNG MULTIKINGDOM###

###DATA PREP###

#Load data
fung_cpz_lung <- read_csv("Fluco_CPZ_Lung_ITS_ASV.csv")
fung.cpz.lung.taxa <- read_csv("Fluco_CPZ_Lung_ITS_Taxa.csv")

fung_cpz_lung <- fung_cpz_lung %>%
  tibble::column_to_rownames("Name")

fung.cpz.lung.taxa <- fung.cpz.lung.taxa %>%
  tibble::column_to_rownames("Name")

bact_cpz_lung <- read_csv("Fluco_CPZ_Lung_16S_ASV.csv")
bact.cpz.lung.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

bact_cpz_lung <- bact_cpz_lung %>%
  tibble::column_to_rownames("Name")

bact.cpz.lung.taxa <- bact.cpz.lung.taxa %>%
  tibble::column_to_rownames("Name")

roughmetadata_combined_cpz_lung <- read_csv("CPZ_Lung_Metadata.csv")

#Renaming the ITS ASVs so they don't overlap with 16S ASV names
row.names(fung_cpz_lung) <- paste("ASV",5386:5658,sep="")
bact_cpz_lung <-as.data.frame(t(bact_cpz_lung))

combined <- rbind.fill(bact_cpz_lung,fung_cpz_lung)
rownames(combined) <- c(row.names(bact_cpz_lung),row.names(fung_cpz_lung))
combined[is.na(combined)] <- 0

#Combine 16S and ITS taxonomy tables, rename the ASVs in the taxonomy table to match the new ITS names
row.names(fung.cpz.lung.taxa) <- paste("ASV",5386:5658,sep="")
combined.taxa <- rbind.fill(as.data.frame(bact.cpz.lung.taxa),as.data.frame(fung.cpz.lung.taxa))
rownames(combined.taxa) <- c(row.names(bact.cpz.lung.taxa),row.names(fung.cpz.lung.taxa))

combined.taxa <- as.matrix(combined.taxa)

combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined_cpz_lung[which(unlist(roughmetadata_combined_cpz_lung$Name) %in% colnames(combined.sort)),]

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
fig6_lung_combined_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)


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
fig6_combined_cpz_lung <- phyloseq(OTU, TAX, samples)



###Fig 6B: LUNG SHANNON ALPHA DIVERSITY###

fr_combined_cpz_lung <- prune_taxa(taxa_sums(fig6_lung_combined_rough) > 0, fig6_lung_combined_rough)
richness_est_combined_cpz_lung <- estimate_richness(fr_combined_cpz_lung, measures = c("Shannon"))

wilcox_alpha_combined_cpz_lung <- t(sapply(richness_est_combined_cpz_lung, function(x) unlist(kruskal.test(x~sample_data(fr_combined_cpz_lung)$Sample.Name)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined_cpz_lung

richness_est_combined_cpz_lung <- richness_est_combined_cpz_lung %>%
  mutate(
    Organ = md.combined.sort2$Organ,
    Oxygen = md.combined.sort2$Oxygen,
    FMT = md.combined.sort2$FMT,
    Sample.Name  = md.combined.sort2$Sample.Name
  )

#Save data table and open it in Graphpad to create plot


###Fig 6B: LUNG BRAY-CURTIS BETA DIVERSITY PCOA###

fig6_combined_cpz_lung_rel <- transform_sample_counts(fig6_combined_cpz_lung, function(x) x / sum(x) )

fig6_combined_cpz_lung_otu_rel <- as.data.frame(t(fig6_combined_cpz_lung_rel@otu_table))
fig6_combined_cpz_lung_tax_rel <- as.data.frame(fig6_combined_cpz_lung_rel@tax_table)
fig6_combined_cpz_lung_meta_rel <- as.data.frame(fig6_combined_cpz_lung_rel@sam_data)

#Get Bray-Curtis dissimilarity matrix
fig6_combined_cpz_lung_bray = vegdist(fig6_combined_cpz_lung_otu_rel, method='bray')

fig6_combined_cpz_lung_pcoa <- ape::pcoa(fig6_combined_cpz_lung_bray)
fig6_combined_cpz_lung_pcoa$values

#This sets the color and factor order for the PCoA and PCA
factor_fig6 <- as.factor(fig6_combined_cpz_lung_meta_rel$Sample.Name)
factor_fig6 <- ordered(factor_fig6, levels = c("CPZ C.trop HO","SPF HO","CPZ C.trop NO", "SPF NO"))
type_fig6 <- as.numeric(factor_fig6)
pca_colors1 <- c("#007016","#B7B7B7")
pca_colors2 <- c("#000000","#000000","#007016","#000000")
ellipse_colors <- c("#007016","#000000","#007016","#000000")

#Plot PCoA
dpi=600
tiff("Fig 6B PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.5), c(-0.4, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(fig6_combined_cpz_lung_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors1[type_fig6], lwd = 1)
points(fig6_combined_cpz_lung_pcoa$vectors[,1:2], pch = 1, cex = 1.3, col = pca_colors2[type_fig6], lwd = 1)
ordiellipse(fig6_combined_cpz_lung_pcoa$vectors[,1:2], factor_fig6,kind = "se", conf = 0.95, col = ellipse_colors)
ordispider(fig6_combined_cpz_lung_pcoa$vectors[,1:2], factor_fig6, label = TRUE)
dev.off()

#Significance test with PERMANOVA, pairwise PERMANOVA, and PERMDISP
permanova_pcoa_df <- data.frame(fig6_combined_cpz_lung_meta_rel)
set.seed(1312)
permanova_combined_cpz_lung_pcoa <- vegan::adonis2(fig6_combined_cpz_lung_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_cpz_lung_pcoa <- permanova_pairwise(fig6_combined_cpz_lung_otu_rel, grp = permanova_pcoa_df$Sample.Name, permutations = 10000, method = "bray", padj = "fdr")

combined_cpz_lung_pcoa_dist <- vegdist(fig6_combined_cpz_lung_otu_rel, method = "bray")
disp_combined_cpz_lung_pcoa <- betadisper(combined_cpz_lung_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_combined_cpz_lung_pcoa <- permutest(disp_combined_cpz_lung_pcoa, permutations = 10000)

print(permanova_combined_cpz_lung_pcoa)
print(pairwise_permanova_combined_cpz_lung_pcoa)
print(permdisp_combined_cpz_lung_pcoa)


###Fig 6B: LUNG BETA DIVERSITY PCA###

fig6_combined_cpz_lung_rel <- transform_sample_counts(fig6_combined_cpz_lung, function(x) x / sum(x) )

fig6_combined_cpz_lung_otu <- as.data.frame(t(fig6_combined_cpz_lung_rel@otu_table))
fig6_combined_cpz_lung_tax <- as.data.frame(fig6_combined_cpz_lung_rel@tax_table)
fig6_combined_cpz_lung_meta <- as.data.frame(fig6_combined_cpz_lung_rel@sam_data)

#Hellinger transform and ordinate
fig6_combined_cpz_lung_otu_hel <- decostand(fig6_combined_cpz_lung_otu, "hellinger")
fig6_combined_cpz_lung_otu_pca <- rda(fig6_combined_cpz_lung_otu_hel)
summary(fig6_combined_cpz_lung_otu_pca)$cont


#Sets which taxa labels should be shown in the loading plot and which should be hidden, to prevent overcrowding
priority_fig6_combined_cpz_lung <- colSums(fig6_combined_cpz_lung_otu_hel)
labels_fig6_combined_cpz_lung <- orditorp(fig6_combined_cpz_lung_otu_pca, "sp", label = fig6_combined_cpz_lung_tax$Genus, priority=priority_fig6_combined_cpz_lung)

#PCA Loading Plot
dpi = 600
tiff("Fig 6B PCA Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(fig6_combined_cpz_lung_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (16.84% Explained)", ylab="PC2 (12.87% Explained)", ylim = c(-0.4,0.5),xlim = c(-0.5,0.6))
orditorp(fig6_combined_cpz_lung_otu_pca, "sp", label = fig6_combined_cpz_lung_tax$Genus, priority=priority_fig6_combined_cpz_lung, select = (labels_fig6_combined_cpz_lung == TRUE), cex = 0.7)
dev.off()

#PCA sample plot
#Color and factor mappings are the same as the PCoA
tiff("Fig 6B PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.9,0.8),c(-0.9,0.8), font = 2, font.lab = 2, xlab="PC1 (16.84% Explained)", ylab="PC2 (12.87% Explained)", type="n")
points(fig6_combined_cpz_lung_otu_pca, pch = 21, cex = 1.3, bg = pca_colors1[factor_fig6], lwd = 1)
points(fig6_combined_cpz_lung_otu_pca, pch = 1, cex = 1.3, col = pca_colors2[factor_fig6], lwd = 1)
ordiellipse(fig6_combined_cpz_lung_otu_pca, factor_fig6, col = ellipse_colors)
ordispider(fig6_combined_cpz_lung_otu_pca, factor_fig6, label = TRUE)
dev.off()

#PERMANOVA, pairwise PERMANOVA, and PERMDISP based on Euclidean distances
permanova_cpz_lung_df <- data.frame(fig6_combined_cpz_lung_meta)
set.seed(1312)
permanova_cpz_lung <- vegan::adonis2(fig6_combined_cpz_lung_otu_hel ~ Sample.Name, data = permanova_cpz_lung_df, method="euclidean", permutations = 10000)
pairwise_permanova_cpz_lung <- permanova_pairwise(fig6_combined_cpz_lung_otu_hel, grp = permanova_cpz_lung_df$Sample.Name, permutations = 10000, method = "euclidean", padj = "fdr")

cpz_lung_hel_dist <- vegdist(fig6_combined_cpz_lung_otu_hel, method = "euclidean")
disp_cpz_lung <- betadisper(cpz_lung_hel_dist, permanova_cpz_lung_df$Sample.Name)
set.seed(1312)
permdisp_cpz_lung <- permutest(disp_cpz_lung, permutations = 10000)

print(permanova_cpz_lung)
print(pairwise_permanova_cpz_lung)
print(permdisp_cpz_lung)


###Fig 6C: GUT MULTIKINGDOM###

###DATA PREP###

#Load data
fung_cpz_gut <- read_csv("Fluco_CPZ_Gut_ITS_ASV.csv")
fung.cpz.gut.taxa <- read_csv("Fluco_CPZ_Gut_ITS_Taxa.csv")

fung_cpz_gut <- fung_cpz_gut %>%
  tibble::column_to_rownames("Name")

fung.cpz.gut.taxa <- fung.cpz.gut.taxa %>%
  tibble::column_to_rownames("Name")

bact_cpz_gut <- read_csv("Fluco_CPZ_Gut_16S_ASV.csv")
bact.cpz.gut.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

bact_cpz_gut <- bact_cpz_gut %>%
  tibble::column_to_rownames("Name")

bact.cpz.gut.taxa <- bact.cpz.gut.taxa %>%
  tibble::column_to_rownames("Name")

roughmetadata_combined_cpz_gut <- read_csv("CPZ_Gut_Metadata.csv")

#Renaming the ITS ASVs so they don't overlap with 16S ASV names
row.names(fung_cpz_gut) <- paste("ASV",5386:5434,sep="")
bact_cpz_gut <-as.data.frame(t(bact_cpz_gut))

combined <- rbind.fill(bact_cpz_gut,fung_cpz_gut)
rownames(combined) <- c(row.names(bact_cpz_gut),row.names(fung_cpz_gut))
combined[is.na(combined)] <- 0

#Combine 16S and ITS taxonomy tables, rename the ASVs in the taxonomy table to match the new ITS names
row.names(fung.cpz.gut.taxa) <- paste("ASV",5386:5434,sep="")
combined.taxa <- rbind.fill(as.data.frame(bact.cpz.gut.taxa),as.data.frame(fung.cpz.gut.taxa))
rownames(combined.taxa) <- c(row.names(bact.cpz.gut.taxa),row.names(fung.cpz.gut.taxa))

combined.taxa <- as.matrix(combined.taxa)

combined.sort <- combined[,order(colnames(combined))]
md.combined <- roughmetadata_combined_cpz_gut[which(unlist(roughmetadata_combined_cpz_gut$Name) %in% colnames(combined.sort)),]

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
fig6_gut_combined_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)


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
fig6_combined_cpz_gut <- phyloseq(OTU, TAX, samples)



###Fig 6C: GUT SHANNON ALPHA DIVERSITY###

fr_combined_cpz_gut <- prune_taxa(taxa_sums(fig6_gut_combined_rough) > 0, fig6_gut_combined_rough)
richness_est_combined_cpz_gut <- estimate_richness(fr_combined_cpz_gut, measures = c("Shannon"))

wilcox_alpha_combined_cpz_gut <- t(sapply(richness_est_combined_cpz_gut, function(x) unlist(kruskal.test(x~sample_data(fr_combined_cpz_gut)$Sample.Name)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined_cpz_gut

richness_est_combined_cpz_gut <- richness_est_combined_cpz_gut %>%
  mutate(
    Organ = md.combined.sort2$Organ,
    Oxygen = md.combined.sort2$Oxygen,
    FMT = md.combined.sort2$FMT,
    Sample.Name  = md.combined.sort2$Sample.Name
  )

#Save data table and open it in Graphpad to create plot


###Fig 6C: GUT BRAY-CURTIS BETA DIVERSITY PCOA###

fig6_combined_cpz_gut_rel <- transform_sample_counts(fig6_combined_cpz_gut, function(x) x / sum(x) )

fig6_combined_cpz_gut_otu_rel <- as.data.frame(t(fig6_combined_cpz_gut_rel@otu_table))
fig6_combined_cpz_gut_tax_rel <- as.data.frame(fig6_combined_cpz_gut_rel@tax_table)
fig6_combined_cpz_gut_meta_rel <- as.data.frame(fig6_combined_cpz_gut_rel@sam_data)

#Get Bray-Curtis dissimilarity matrix
fig6_combined_cpz_gut_bray = vegdist(fig6_combined_cpz_gut_otu_rel, method='bray')

fig6_combined_cpz_gut_pcoa <- ape::pcoa(fig6_combined_cpz_gut_bray)
fig6_combined_cpz_gut_pcoa$values

#This sets the color and factor order for the PCoA and PCA
factor_fig6 <- as.factor(fig6_combined_cpz_gut_meta_rel$Sample.Name)
factor_fig6 <- ordered(factor_fig6, levels = c("CPZ C.trop HO","SPF HO","CPZ C.trop NO", "SPF NO"))
type_fig6 <- as.numeric(factor_fig6)
pca_colors1 <- c("#007016","#B7B7B7")
pca_colors2 <- c("#000000","#000000","#007016","#000000")
ellipse_colors <- c("#007016","#000000","#007016","#000000")

#Plot PCoA
dpi=600
tiff("Fig 6C PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.8, 0.4), c(-0.6, 0.3), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(fig6_combined_cpz_gut_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors1[type_fig6], lwd = 1)
points(fig6_combined_cpz_gut_pcoa$vectors[,1:2], pch = 1, cex = 1.3, col = pca_colors2[type_fig6], lwd = 1)
ordiellipse(fig6_combined_cpz_gut_pcoa$vectors[,1:2], factor_fig6,kind = "se", conf = 0.95, col = ellipse_colors)
ordispider(fig6_combined_cpz_gut_pcoa$vectors[,1:2], factor_fig6, label = TRUE)
dev.off()

#Significance test with PERMANOVA, pairwise PERMANOVA, and PERMDISP
permanova_pcoa_df <- data.frame(fig6_combined_cpz_gut_meta_rel)
set.seed(1312)
permanova_combined_cpz_gut_pcoa <- vegan::adonis2(fig6_combined_cpz_gut_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_combined_cpz_gut_pcoa <- permanova_pairwise(fig6_combined_cpz_gut_otu_rel, grp = permanova_pcoa_df$Sample.Name, permutations = 10000, method = "bray", padj = "fdr")

combined_cpz_gut_pcoa_dist <- vegdist(fig6_combined_cpz_gut_otu_rel, method = "bray")
disp_combined_cpz_gut_pcoa <- betadisper(combined_cpz_gut_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_combined_cpz_gut_pcoa <- permutest(disp_combined_cpz_gut_pcoa, permutations = 10000)

print(permanova_combined_cpz_gut_pcoa)
print(pairwise_permanova_combined_cpz_gut_pcoa)
print(permdisp_combined_cpz_gut_pcoa)


###Fig 6C: GUT BETA DIVERSITY PCA###

fig6_combined_cpz_gut_rel <- transform_sample_counts(fig6_combined_cpz_gut, function(x) x / sum(x) )

fig6_combined_cpz_gut_otu <- as.data.frame(t(fig6_combined_cpz_gut_rel@otu_table))
fig6_combined_cpz_gut_tax <- as.data.frame(fig6_combined_cpz_gut_rel@tax_table)
fig6_combined_cpz_gut_meta <- as.data.frame(fig6_combined_cpz_gut_rel@sam_data)

#Hellinger transform and ordinate
fig6_combined_cpz_gut_otu_hel <- decostand(fig6_combined_cpz_gut_otu, "hellinger")
fig6_combined_cpz_gut_otu_pca <- rda(fig6_combined_cpz_gut_otu_hel)
summary(fig6_combined_cpz_gut_otu_pca)$cont


#Sets which taxa labels should be shown in the loading plot and which should be hidden, to prevent overcrowding
priority_fig6_combined_cpz_gut <- colSums(fig6_combined_cpz_gut_otu_hel)
labels_fig6_combined_cpz_gut <- orditorp(fig6_combined_cpz_gut_otu_pca, "sp", label = fig6_combined_cpz_gut_tax$Genus, priority=priority_fig6_combined_cpz_gut)

#PCA Loading Plot
dpi = 600
tiff("Fig 6C PCA Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(fig6_combined_cpz_gut_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (53.78% Explained)", ylab="PC2 (22.75% Explained)", ylim = c(-0.4,0.8),xlim = c(-0.9,0.8))
orditorp(fig6_combined_cpz_gut_otu_pca, "sp", label = fig6_combined_cpz_gut_tax$Genus, priority=priority_fig6_combined_cpz_gut, select = (labels_fig6_combined_cpz_gut == TRUE), cex = 0.7)
dev.off()

#PCA sample plot
#Color and factor mappings are the same as the PCoA
tiff("Fig 6C PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4,0.9),c(-0.6,0.7), font = 2, font.lab = 2, xlab="PC1 (53.78% Explained)", ylab="PC2 (22.75% Explained)", type="n")
points(fig6_combined_cpz_gut_otu_pca, pch = 21, cex = 1.3, bg = pca_colors1[factor_fig6], lwd = 1)
points(fig6_combined_cpz_gut_otu_pca, pch = 1, cex = 1.3, col = pca_colors2[factor_fig6], lwd = 1)
ordiellipse(fig6_combined_cpz_gut_otu_pca, factor_fig6, col = ellipse_colors)
ordispider(fig6_combined_cpz_gut_otu_pca, factor_fig6, label = TRUE)
dev.off()

#PERMANOVA, pairwise PERMANOVA, and PERMDISP based on Euclidean distances
permanova_cpz_gut_df <- data.frame(fig6_combined_cpz_gut_meta)
set.seed(1312)
permanova_cpz_gut <- vegan::adonis2(fig6_combined_cpz_gut_otu_hel ~ Sample.Name, data = permanova_cpz_gut_df, method="euclidean", permutations = 10000)
pairwise_permanova_cpz_gut <- permanova_pairwise(fig6_combined_cpz_gut_otu_hel, grp = permanova_cpz_gut_df$Sample.Name, permutations = 10000, method = "euclidean", padj = "fdr")

cpz_gut_hel_dist <- vegdist(fig6_combined_cpz_gut_otu_hel, method = "euclidean")
disp_cpz_gut <- betadisper(cpz_gut_hel_dist, permanova_cpz_gut_df$Sample.Name)
set.seed(1312)
permdisp_cpz_gut <- permutest(disp_cpz_gut, permutations = 10000)

print(permanova_cpz_gut)
print(pairwise_permanova_cpz_gut)
print(permdisp_cpz_gut)
print(permanova_bact_mouse_pcoa)
print(pairwise_permanova_bact_mouse_pcoa)
print(permdisp_bact_mouse_pcoa)
