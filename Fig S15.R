library(tidyverse)
library(phyloseq)
library(plyr)
library(vegan)
library(ecole)

###DEFINE FUNCTION###
#This function does some basic filtering and creates the phyloseq object
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

#This function performs alpha diversity significance testing
alphadiv <- function(physeq, vari){
  
  #Remove OTUs with all 0s, they mess with analysis and are not necessary
  physeq_zerorm <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  #Get Simpson and Shannon indices
  richness <- estimate_richness(physeq_zerorm, measures = c("Simpson", "Shannon"))
  kruskal <- t(sapply(richness, function(x) unlist(kruskal.test(x~physeq_zerorm@sam_data[[vari]])[c("estimate","p.value","statistic","conf.int")])))
  print(kruskal) 
  return(richness)
  
}

###DATA PREP###

#Load data
fung_cpz_gut <- read_csv("Fluco_CPZ_Gut_ITS_ASV.csv")
fung.cpz.gut.taxa <- read_csv("Fluco_CPZ_Gut_ITS_Taxa.csv")

bact_cpz_gut <- read_csv("Fluco_CPZ_Gut_16S_ASV.csv")
bact.cpz.gut.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

roughmetadata_combined_cpz_gut <- read_csv("CPZ_Gut_Metadata.csv")

fung_cpz_lung <- read_csv("Fluco_CPZ_Lung_ITS_ASV.csv")
fung.cpz.lung.taxa <- read_csv("Fluco_CPZ_Lung_ITS_Taxa.csv")

bact_cpz_lung <- read_csv("Fluco_CPZ_Lung_16S_ASV.csv")
bact.cpz.lung.taxa <- read_csv("Fluco_CPZ_16S_Taxa.csv")

roughmetadata_combined_cpz_lung <- read_csv("CPZ_Lung_Metadata.csv")

#Gut ITS
exp4_gut_fung_rough <- dada2_to_phyloseq(fung_cpz_gut,fung.cpz.gut.taxa,roughmetadata_combined_cpz_gut, read_count = 50, rough = TRUE)
exp4_gut_fung <- dada2_to_phyloseq(fung_cpz_gut,fung.cpz.gut.taxa,roughmetadata_combined_cpz_gut, read_count =50, rough = FALSE)

#Gut 16S
exp4_gut_bact_rough <- dada2_to_phyloseq(bact_cpz_gut,bact.cpz.gut.taxa,roughmetadata_combined_cpz_gut, read_count = 1000,samples_are_cols=FALSE, rough = TRUE)
exp4_gut_bact <- dada2_to_phyloseq(bact_cpz_gut,bact.cpz.gut.taxa,roughmetadata_combined_cpz_gut, read_count =1000, samples_are_cols=FALSE,rough = FALSE)

#Lung ITS
exp4_lung_fung_rough <- dada2_to_phyloseq(fung_cpz_lung,fung.cpz.lung.taxa,roughmetadata_combined_cpz_lung, read_count = 50, rough = TRUE)
exp4_lung_fung <- dada2_to_phyloseq(fung_cpz_lung,fung.cpz.lung.taxa,roughmetadata_combined_cpz_lung, read_count =50, rough = FALSE)

#Lung 16S
exp4_lung_bact_rough <- dada2_to_phyloseq(bact_cpz_lung,bact.cpz.lung.taxa,roughmetadata_combined_cpz_lung, read_count = 50,samples_are_cols=FALSE, rough = TRUE)
exp4_lung_bact <- dada2_to_phyloseq(bact_cpz_lung,bact.cpz.lung.taxa,roughmetadata_combined_cpz_lung, read_count =50, samples_are_cols=FALSE,rough = FALSE)


###Fig. S15A: CPZ C.TROP GUT ITS ALPHA AND BETA DIVERSITY###

###Alpha Diversity###

alphadiv_fung_cpz_gut <- alphadiv(exp4_gut_fung_rough, "Sample.Name")

alphadiv_fung_cpz_gut <- alphadiv_fung_cpz_gut %>%
  mutate(
    Organ = exp4_gut_fung_rough@sam_data[["Organ"]],
    Oxygen = exp4_gut_fung_rough@sam_data[["condition"]],
    Treatment = exp4_gut_fung_rough@sam_data[["treatment"]],
    Sample.Name = exp4_gut_fung_rough@sam_data[["Sample.Name"]]
  )

#Save

###Beta Diversity###

exp4_gut_fung_rel <- transform_sample_counts(exp4_gut_fung, function(x) x / sum(x) )

exp4_gut_fung_otu_rel <- as.data.frame(t(exp4_gut_fung_rel@otu_table))
exp4_gut_fung_tax_rel <- as.data.frame(exp4_gut_fung_rel@tax_table)
exp4_gut_fung_meta_rel <- as.data.frame(exp4_gut_fung_rel@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp4_gut_fung_rel_bray = vegdist(exp4_gut_fung_otu_rel, method='bray')
exp4_gut_fung_rel_pcoa <- ape::pcoa(exp4_gut_fung_rel_bray)
exp4_gut_fung_rel_pcoa$values

#Colors points according to BPD status
factor_exp4 <- as.factor(exp4_gut_fung_meta_rel$Sample.Name)
type_exp4 <- as.numeric(factor_exp4)
pca_colors <- c("#007016","#FFFFFF","#B3B3B3","#FFFFFF")
pca_colors2 <- c("#000000","#007016","#000000","#B3B3B3")
ellipse_colors <- c("#007016","#007016","#B3B3B3","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig S15A PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.7, 0.5), c(-0.8, 1), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp4_gut_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp4], lwd = 1)
points(exp4_gut_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp4], lwd = 1)
ordiellipse(exp4_gut_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp4, col = pca_colors2)
ordispider(exp4_gut_fung_rel_pcoa$vectors[,1:2], factor_exp4, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp4_gut_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp4_gut_fung_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp4_gut_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(permdisp_fung_pcoa)



###Fig. S15B: CPZ C.TROP GUT 16S ALPHA AND BETA DIVERSITY###

###Alpha Diversity###

alphadiv_bact_cpz_gut <- alphadiv(exp4_gut_bact_rough, "Sample.Name")

alphadiv_bact_cpz_gut <- alphadiv_bact_cpz_gut %>%
  mutate(
    Organ = exp4_gut_bact_rough@sam_data[["Organ"]],
    Oxygen = exp4_gut_bact_rough@sam_data[["condition"]],
    Treatment = exp4_gut_bact_rough@sam_data[["treatment"]],
    Sample.Name = exp4_gut_bact_rough@sam_data[["Sample.Name"]]
  )

#Save

###Beta Diversity###

exp4_gut_bact_rel <- transform_sample_counts(exp4_gut_bact, function(x) x / sum(x) )

exp4_gut_bact_otu_rel <- as.data.frame(t(exp4_gut_bact_rel@otu_table))
exp4_gut_bact_tax_rel <- as.data.frame(exp4_gut_bact_rel@tax_table)
exp4_gut_bact_meta_rel <- as.data.frame(exp4_gut_bact_rel@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp4_gut_bact_rel_bray = vegdist(exp4_gut_bact_otu_rel, method='bray')
exp4_gut_bact_rel_pcoa <- ape::pcoa(exp4_gut_bact_rel_bray)
exp4_gut_bact_rel_pcoa$values

#Colors points according to BPD status
factor_exp4 <- as.factor(exp4_gut_bact_meta_rel$Sample.Name)
type_exp4 <- as.numeric(factor_exp4)
pca_colors <- c("#007016","#FFFFFF","#B3B3B3","#FFFFFF")
pca_colors2 <- c("#000000","#007016","#000000","#B3B3B3")
ellipse_colors <- c("#007016","#007016","#B3B3B3","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig S15B PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.3, 0.8), c(-0.6, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp4_gut_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp4], lwd = 1)
points(exp4_gut_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp4], lwd = 1)
ordiellipse(exp4_gut_bact_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp4, col = ellipse_colors)
ordispider(exp4_gut_bact_rel_pcoa$vectors[,1:2], factor_exp4, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp4_gut_bact_meta_rel)
set.seed(1312)
permanova_bact_pcoa <- vegan::adonis2(exp4_gut_bact_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
bact_pcoa_dist <- vegdist(exp4_gut_bact_otu_rel, method = "bray")
disp_bact_pcoa <- betadisper(bact_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_bact_pcoa <- permutest(disp_bact_pcoa, permutations = 10000)

print(permanova_bact_pcoa)
print(permdisp_bact_pcoa)



###Fig. S15C: CPZ C.TROP LUNG ITS ALPHA AND BETA DIVERSITY###

###Alpha Diversity###

alphadiv_fung_cpz_lung <- alphadiv(exp4_lung_fung_rough, "Sample.Name")

alphadiv_fung_cpz_lung <- alphadiv_fung_cpz_lung %>%
  mutate(
    Organ = exp4_lung_fung_rough@sam_data[["Organ"]],
    Oxygen = exp4_lung_fung_rough@sam_data[["condition"]],
    Treatment = exp4_lung_fung_rough@sam_data[["treatment"]],
    Sample.Name = exp4_lung_fung_rough@sam_data[["Sample.Name"]]
  )

#Save

###Beta Diversity###

exp4_lung_fung_rel <- transform_sample_counts(exp4_lung_fung, function(x) x / sum(x) )

exp4_lung_fung_otu_rel <- as.data.frame(t(exp4_lung_fung_rel@otu_table))
exp4_lung_fung_tax_rel <- as.data.frame(exp4_lung_fung_rel@tax_table)
exp4_lung_fung_meta_rel <- as.data.frame(exp4_lung_fung_rel@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp4_lung_fung_rel_bray = vegdist(exp4_lung_fung_otu_rel, method='bray')
exp4_lung_fung_rel_pcoa <- ape::pcoa(exp4_lung_fung_rel_bray)
exp4_lung_fung_rel_pcoa$values

#Colors points according to BPD status
factor_exp4 <- as.factor(exp4_lung_fung_meta_rel$`Sample Name`)
type_exp4 <- as.numeric(factor_exp4)
pca_colors <- c("#007016","#FFFFFF","#B3B3B3","#FFFFFF")
pca_colors2 <- c("#000000","#007016","#000000","#B3B3B3")
ellipse_colors <- c("#007016","#007016","#B3B3B3","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig 13C PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.6, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp4_lung_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp4], lwd = 1)
points(exp4_lung_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp4], lwd = 1)
ordiellipse(exp4_lung_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp4, col = ellipse_colors)
ordispider(exp4_lung_fung_rel_pcoa$vectors[,1:2], factor_exp4, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp4_lung_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp4_lung_fung_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp4_lung_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(permdisp_fung_pcoa)



###Fig. S15D: CPZ C.TROP LUNG 16S ALPHA AND BETA DIVERSITY###

###Alpha Diversity###

alphadiv_bact_cpz_lung <- alphadiv(exp4_lung_bact_rough, "Sample.Name")

alphadiv_bact_cpz_lung <- alphadiv_bact_cpz_lung %>%
  mutate(
    Organ = exp4_lung_bact_rough@sam_data[["Organ"]],
    Oxygen = exp4_lung_bact_rough@sam_data[["condition"]],
    Treatment = exp4_lung_bact_rough@sam_data[["treatment"]],
    Sample.Name = exp4_lung_bact_rough@sam_data[["Sample.Name"]]
  )

#Save

###Beta Diversity###

exp4_lung_bact_rel <- transform_sample_counts(exp4_lung_bact, function(x) x / sum(x) )

exp4_lung_bact_otu_rel <- as.data.frame(t(exp4_lung_bact_rel@otu_table))
exp4_lung_bact_tax_rel <- as.data.frame(exp4_lung_bact_rel@tax_table)
exp4_lung_bact_meta_rel <- as.data.frame(exp4_lung_bact_rel@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp4_lung_bact_rel_bray = vegdist(exp4_lung_bact_otu_rel, method='bray')
exp4_lung_bact_rel_pcoa <- ape::pcoa(exp4_lung_bact_rel_bray)
exp4_lung_bact_rel_pcoa$values

#Colors points according to BPD status
factor_exp4 <- as.factor(exp4_lung_bact_meta_rel$Sample.Name)
type_exp4 <- as.numeric(factor_exp4)
pca_colors <- c("#007016","#FFFFFF","#B3B3B3","#FFFFFF")
pca_colors2 <- c("#000000","#007016","#000000","#B3B3B3")
ellipse_colors <- c("#007016","#007016","#B3B3B3","#B3B3B3")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig S15D PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.5), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp4_lung_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp4], lwd = 1)
points(exp4_lung_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors2[type_exp4], lwd = 1)
ordiellipse(exp4_lung_bact_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp4, col = ellipse_colors)
ordispider(exp4_lung_bact_rel_pcoa$vectors[,1:2], factor_exp4, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp4_lung_bact_meta_rel)
set.seed(1312)
permanova_bact_pcoa <- vegan::adonis2(exp4_lung_bact_otu_rel ~ Sample.Name, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
bact_pcoa_dist <- vegdist(exp4_lung_bact_otu_rel, method = "bray")
disp_bact_pcoa <- betadisper(bact_pcoa_dist, permanova_pcoa_df$Sample.Name)
set.seed(1312)
permdisp_bact_pcoa <- permutest(disp_bact_pcoa, permutations = 10000)

print(permanova_bact_pcoa)
print(permdisp_bact_pcoa)
