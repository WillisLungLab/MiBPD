library(tidyverse)
library(phyloseq)
library(plyr)
library(vegan)
library(ecole)

###DEFINE FUNCTION###
dada2_to_phyloseq <- function(asv,taxa,meta, read_count=100,prev_cutoff=0.05,samples_are_cols=TRUE,rough=FALSE){
  asv <- asv %>%
    tibble::column_to_rownames("Name")
  if(samples_are_cols==FALSE){
    asv <- as.data.frame(t(asv))
  }
  
  #Load taxonomy table
  #Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
  taxa <- taxa %>%
    tibble::column_to_rownames("Name")
  taxa <- as.matrix(taxa)
  
  #Load metadata
  asv.sort <- asv[,order(colnames(asv))]
  md.sort <- meta[which(unlist(meta$Name) %in% colnames(asv.sort)),]
  print(summary(colSums(asv.sort)))
  
  #Create phyloseq object with unfiltered data for alpha diversity
  asv.sort2 <- asv.sort[,c(which(colSums(asv.sort)>= read_count))]
  md.sort <- md.sort %>%
    tibble::column_to_rownames("Name")
  md.sort2 <- md.sort[which(row.names(md.sort) %in% colnames(asv.sort2)),]
  
  if(rough==TRUE){
    OTU_rough = otu_table(asv.sort2, taxa_are_rows = TRUE)
    TAX_rough = tax_table(taxa)
    samples_rough = sample_data(md.sort)
    pseq_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)
    return(pseq_rough)
    
  } else if(rough==FALSE) {
    #This is the same
    gt0 <- function(vec){
      v <- as.numeric(vec)
      s <- sum(v>0)
      return(s)
    }
    
    asv.nsampls <- apply(asv.sort2, 1, gt0)
    asv.clean <- asv.sort2[which(asv.nsampls>1),]
    asv.clean.colsum <- asv.clean[,colSums(asv.clean)>0]
    
    
    #Setup to create the phyloseq object
    OTU = otu_table(asv.clean.colsum, taxa_are_rows = TRUE)
    TAX = tax_table(taxa)
    samples = sample_data(md.sort2)
    
    #Convert!
    pseq <- phyloseq(OTU, TAX, samples)
    return(pseq)
  } else {
    print("True/False value needed for variable 'rough'")
  }
}


###DATA PREP###

#Load data
fung_fluco_gut <- read_csv("Fluco_CPZ_Gut_ITS_ASV.csv")
fung.fluco.gut.taxa <- read_csv("Fluco_CPZ_Gut_ITS_Taxa.csv")

bact_fluco_gut <- read_csv("Fluco_CPZ_Gut_16S_ASV.csv")
bact.fluco.gut.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

roughmetadata_combined_fluco_gut <- read_csv("Fluco_Gut_Metadata.csv")

fung_fluco_lung <- read_csv("Fluco_CPZ_Lung_ITS_ASV.csv")
fung.fluco.lung.taxa <- read_csv("Fluco_CPZ_Lung_ITS_Taxa.csv")

bact_fluco_lung <- read_csv("Fluco_CPZ_Lung_16S_ASV.csv")
bact.fluco.lung.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

roughmetadata_combined_fluco_lung <- read_csv("Fluco_Lung_Metadata.csv")

#Gut ITS
exp3_gut_fung_rough <- dada2_to_phyloseq(fung_fluco_gut,fung.fluco.gut.taxa,roughmetadata_combined_fluco_gut, read_count = 50, rough = TRUE)
exp3_gut_fung <- dada2_to_phyloseq(fung_fluco_gut,fung.fluco.gut.taxa,roughmetadata_combined_fluco_gut, read_count =50, rough = FALSE)
exp3_gut_fung_prev <- filter_taxa(exp3_gut_fung, function(x) sum(x >= 1) > (0.05*length(x)), TRUE)

#Gut 16S
exp3_gut_bact_rough <- dada2_to_phyloseq(bact_fluco_gut,bact.fluco.gut.taxa,roughmetadata_combined_fluco_gut, read_count = 1000,samples_are_cols=FALSE, rough = TRUE)
exp3_gut_bact <- dada2_to_phyloseq(bact_fluco_gut,bact.fluco.gut.taxa,roughmetadata_combined_fluco_gut, read_count =1000, samples_are_cols=FALSE,rough = FALSE)
exp3_gut_bact_prev <- filter_taxa(exp3_gut_bact, function(x) sum(x >= 1) > (0.05*length(x)), TRUE)

#Lung ITS
exp3_lung_fung_rough <- dada2_to_phyloseq(fung_fluco_lung,fung.fluco.lung.taxa,roughmetadata_combined_fluco_lung, read_count = 50, rough = TRUE)
exp3_lung_fung <- dada2_to_phyloseq(fung_fluco_lung,fung.fluco.lung.taxa,roughmetadata_combined_fluco_lung, read_count =50, rough = FALSE)
exp3_lung_fung_prev <- filter_taxa(exp3_lung_fung, function(x) sum(x >= 1) > (0.05*length(x)), TRUE)

#Lung 16S
exp3_lung_bact_rough <- dada2_to_phyloseq(bact_fluco_lung,bact.fluco.lung.taxa,roughmetadata_combined_fluco_lung, read_count = 50,samples_are_cols=FALSE, rough = TRUE)
exp3_lung_bact <- dada2_to_phyloseq(bact_fluco_lung,bact.fluco.lung.taxa,roughmetadata_combined_fluco_lung, read_count =50, samples_are_cols=FALSE,rough = FALSE)
exp3_lung_bact_prev <- filter_taxa(exp3_lung_bact, function(x) sum(x >= 1) > (0.05*length(x)), TRUE)


###Fig. S13A: GUT ITS###

exp3_gut_fung_rel_prev <- transform_sample_counts(exp3_gut_fung, function(x) x / sum(x) )

exp3_gut_fung_otu_rel <- as.data.frame(t(exp3_gut_fung_rel_prev@otu_table))
exp3_gut_fung_tax_rel <- as.data.frame(exp3_gut_fung_rel_prev@tax_table)
exp3_gut_fung_meta_rel <- as.data.frame(exp3_gut_fung_rel_prev@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp3_gut_fung_rel_bray = vegdist(exp3_gut_fung_otu_rel, method='bray')
exp3_gut_fung_rel_pcoa <- ape::pcoa(exp3_gut_fung_rel_bray)
exp3_gut_fung_rel_pcoa$values

#Colors points according to BPD status
factor_exp3 <- as.factor(exp3_gut_fung_meta_rel$Sample.Name)
type_exp3 <- as.numeric(factor_exp3)
pca_colors <- c("#FFFFFF","#B3B3B3")
pca_colors2 <- c("#94D57F","#000000","#B3B3B3")
ellipse_colors <- c("#94D57F","#B3B3B3","#94D57F","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig S13A.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.5), c(-0.5, 0.6), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp3_gut_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp3], lwd = 1)
points(exp3_gut_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp3], lwd = 1)
ordiellipse(exp3_gut_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp3, col = pca_colors2)
ordispider(exp3_gut_fung_rel_pcoa$vectors[,1:2], factor_exp3, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df$ <- data.frame(exp3_gut_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp3_gut_fung_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp3_gut_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(permdisp_fung_pcoa)



###Fig. S13B: GUT 16S###

exp3_gut_bact_rel_prev <- transform_sample_counts(exp3_gut_bact_prev, function(x) x / sum(x) )

exp3_gut_bact_otu_rel <- as.data.frame(t(exp3_gut_bact_rel_prev@otu_table))
exp3_gut_bact_tax_rel <- as.data.frame(exp3_gut_bact_rel_prev@tax_table)
exp3_gut_bact_meta_rel <- as.data.frame(exp3_gut_bact_rel_prev@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp3_gut_bact_rel_bray = vegdist(exp3_gut_bact_otu_rel, method='bray')
exp3_gut_bact_rel_pcoa <- ape::pcoa(exp3_gut_bact_rel_bray)
exp3_gut_bact_rel_pcoa$values

#Colors points according to BPD status
factor_exp3 <- as.factor(exp3_gut_bact_meta_rel$Sample.Name)
type_exp3 <- as.numeric(factor_exp3)
pca_colors <- c("#94D57F","#FFFFFF","#B3B3B3","#FFFFFF")
pca_colors2 <- c("#000000","#94D57F","#000000","#B3B3B3")
ellipse_colors <- c("#94D57F","#94D57F","#B3B3B3","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig S13B.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.25, 0.15), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp3_gut_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp3], lwd = 1)
points(exp3_gut_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp3], lwd = 1)
ordiellipse(exp3_gut_bact_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp3, col = ellipse_colors)
ordispider(exp3_gut_bact_rel_pcoa$vectors[,1:2], factor_exp3, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp3_gut_bact_meta_rel)
set.seed(1312)
permanova_bact_pcoa <- vegan::adonis2(exp3_gut_bact_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
bact_pcoa_dist <- vegdist(exp3_gut_bact_otu_rel, method = "bray")
disp_bact_pcoa <- betadisper(bact_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_bact_pcoa <- permutest(disp_bact_pcoa, permutations = 10000)

print(permanova_bact_pcoa)
print(permdisp_bact_pcoa)



###Fig. S13C: LUNG ITS###

exp3_lung_fung_rel_prev <- transform_sample_counts(exp3_lung_fung_prev, function(x) x / sum(x) )

exp3_lung_fung_otu_rel <- as.data.frame(t(exp3_lung_fung_rel_prev@otu_table))
exp3_lung_fung_tax_rel <- as.data.frame(exp3_lung_fung_rel_prev@tax_table)
exp3_lung_fung_meta_rel <- as.data.frame(exp3_lung_fung_rel_prev@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp3_lung_fung_rel_bray = vegdist(exp3_lung_fung_otu_rel, method='bray')
exp3_lung_fung_rel_pcoa <- ape::pcoa(exp3_lung_fung_rel_bray)
exp3_lung_fung_rel_pcoa$values

#Colors points according to BPD status
factor_exp3 <- as.factor(exp3_lung_fung_meta_rel$Sample.Name)
type_exp3 <- as.numeric(factor_exp3)
pca_colors <- c("#94D57F","#B3B3B3")
pca_colors2 <- c("#000000","#000000","#94D57F","#000000")
ellipse_colors <- c("#94D57F","#B3B3B3","#94D57F","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig 13C.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.6), c(-0.5, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp3_lung_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp3], lwd = 1)
points(exp3_lung_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp3], lwd = 1)
ordiellipse(exp3_lung_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp3, col = ellipse_colors)
ordispider(exp3_lung_fung_rel_pcoa$vectors[,1:2], factor_exp3, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp3_lung_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp3_lung_fung_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp3_lung_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(permdisp_fung_pcoa)



###Fig. S13D: LUNG 16S###

exp3_lung_bact_rel_prev <- transform_sample_counts(exp3_lung_bact_prev, function(x) x / sum(x) )

exp3_lung_bact_otu_rel <- as.data.frame(t(exp3_lung_bact_rel_prev@otu_table))
exp3_lung_bact_tax_rel <- as.data.frame(exp3_lung_bact_rel_prev@tax_table)
exp3_lung_bact_meta_rel <- as.data.frame(exp3_lung_bact_rel_prev@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp3_lung_bact_rel_bray = vegdist(exp3_lung_bact_otu_rel, method='bray')
exp3_lung_bact_rel_pcoa <- ape::pcoa(exp3_lung_bact_rel_bray)
exp3_lung_bact_rel_pcoa$values

#Colors points according to BPD status
factor_exp3 <- as.factor(exp3_lung_bact_meta_rel$Sample.Name)
type_exp3 <- as.numeric(factor_exp3)
pca_colors <- c("#94D57F","#FFFFFF","#B3B3B3","#FFFFFF")
pca_colors2 <- c("#000000","#94D57F","#000000","#B3B3B3")
ellipse_colors <- c("#94D57F","#94D57F","#B3B3B3","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig S13D.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.45, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp3_lung_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp3], lwd = 1)
points(exp3_lung_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp3], lwd = 1)
ordiellipse(exp3_lung_bact_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp3, col = ellipse_colors)
ordispider(exp3_lung_bact_rel_pcoa$vectors[,1:2], factor_exp3, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp3_lung_bact_meta_rel)
set.seed(1312)
permanova_bact_pcoa <- vegan::adonis2(exp3_lung_bact_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
bact_pcoa_dist <- vegdist(exp3_lung_bact_otu_rel, method = "bray")
disp_bact_pcoa <- betadisper(bact_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_bact_pcoa <- permutest(disp_bact_pcoa, permutations = 10000)

print(permanova_bact_pcoa)
print(permdisp_bact_pcoa)
