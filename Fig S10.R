###Define functions###

#Cleans up data and creates a phyloseq object
dada2_to_phyloseq <- function(asv,taxa,meta, read_count=100,prev_cutoff=0.05,rough=FALSE){
  asv <- asv %>%
    tibble::column_to_rownames("Name")
  
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
    pseq_prev <- filter_taxa(pseq, function(x) sum(x >= 1) > (prev_cutoff*length(x)), TRUE)
    return(pseq_prev)
  } else {
    print("True/False value needed for variable 'rough'")
  }
}


#Performs alpha diversity significance testing
alphadiv <- function(physeq, vari){
  
  #Remove OTUs with all 0s, they mess with analysis and are not necessary
  physeq_zerorm <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  #Get Simpson and Shannon indices
  richness <- estimate_richness(physeq_zerorm, measures = c("Simpson", "Shannon"))
  kruskal <- t(sapply(richness, function(x) unlist(kruskal.test(x~physeq_zerorm@sam_data[[vari]])[c("estimate","p.value","statistic","conf.int")])))
  print(kruskal) 
  return(richness)
  
}

#############Fig. S10A: Gut ITS#############

###Data Prep###
fung_mouse_gut <- read_csv("Mouse_Gut_ITS_ASVs.csv")
fung.mouse.gut.taxa <- read_csv("Mouse_Gut_ITS_Taxa.csv")
roughmetadata_fung_mouse_gut <- read_csv("Mouse_Gut_Metadata.csv")

###Alpha Diversity###
exp2_fung_mouse_gut_rough <- dada2_to_phyloseq(fung_mouse_gut,fung.mouse.gut.taxa,roughmetadata_fung_mouse_gut, read_count = 50, rough = TRUE)
alphadiv_fung_mouse_gut <- alphadiv(exp2_fung_mouse_gut_rough, "Condition")

alphadiv_fung_mouse_gut <- alphadiv_fung_mouse_gut %>%
  mutate(
    Organ = exp2_fung_mouse_gut_rough@sam_data[["Organ"]],
    Oxygen = exp2_fung_mouse_gut_rough@sam_data[["Oxygen"]],
    FMT = exp2_fung_mouse_gut_rough@sam_data[["FMT"]],
    Condition = exp2_fung_mouse_gut_rough@sam_data[["Condition"]]
  )

#Save


###Beta Diversity###
exp2_fung_mouse_gut_prev <- dada2_to_phyloseq(fung_mouse_gut,fung.mouse.gut.taxa,roughmetadata_fung_mouse_gut, read_count = 50, rough = FALSE)
exp2_fung_mouse_rel_prev <- transform_sample_counts(exp2_fung_mouse_gut_prev, function(x) x / sum(x) )

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
tiff("Fig S10A Mouse Gut ITS PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
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



#############Fig. S10B: Gut 16s#############

bact_mouse_gut <- read_csv("Mouse_Gut_16s_ASVs.csv")
bact.mouse.gut.taxa <- read_csv("Mouse_Gut_16s_Taxa.csv")
roughmetadata_bact_mouse_gut <- read_csv("Mouse_Gut_Metadata.csv")

###Alpha Diversity###
exp2_bact_mouse_gut_rough <- dada2_to_phyloseq(bact_mouse_gut,bact.mouse.gut.taxa,roughmetadata_bact_mouse_gut, read_count = 1000, rough = TRUE)
alphadiv_bact_mouse_gut <- alphadiv(exp2_bact_mouse_gut_rough, "Condition")

alphadiv_bact_mouse_gut <- alphadiv_bact_mouse_gut %>%
  mutate(
    Organ = exp2_bact_mouse_gut_rough@sam_data[["Organ"]],
    Oxygen = exp2_bact_mouse_gut_rough@sam_data[["Oxygen"]],
    FMT = exp2_bact_mouse_gut_rough@sam_data[["FMT"]],
    Condition = exp2_bact_mouse_gut_rough@sam_data[["Condition"]]
  )

#Save


###Beta Diversity###
exp2_bact_mouse_gut <- dada2_to_phyloseq(bact_mouse_gut,bact.mouse.gut.taxa,roughmetadata_bact_mouse_gut, read_count = 1000, rough = FALSE)
exp2_bact_mouse_rel_prev <- transform_sample_counts(exp2_bact_mouse_gut, function(x) x / sum(x) )

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
tiff("Fig S10B Mouse Gut 16S PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
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



#############Fig. S10C: Lung ITS#############

fung_mouse_lung <- read_csv("Mouse_Lung_ITS_ASVs.csv")
fung.mouse.lung.taxa <- read_csv("Mouse_Lung_ITS_Taxa.csv")
roughmetadata_fung_mouse_lung <- read_csv("Mouse_Lung_Metadata.csv")


###Alpha Diversity###

exp2_fung_mouse_lung_rough <- dada2_to_phyloseq(fung_mouse_lung,fung.mouse.lung.taxa,roughmetadata_fung_mouse_lung, read_count = 0, rough = TRUE)
alphadiv_fung_mouse_lung <- alphadiv(exp2_fung_mouse_lung_rough, "Condition")

alphadiv_fung_mouse_lung <- alphadiv_fung_mouse_lung %>%
  mutate(
    Organ = exp2_fung_mouse_lung_rough@sam_data[["Organ"]],
    Oxygen = exp2_fung_mouse_lung_rough@sam_data[["Oxygen"]],
    FMT = exp2_fung_mouse_lung_rough@sam_data[["FMT"]],
    Condition = exp2_fung_mouse_lung_rough@sam_data[["Condition"]]
  )

write.csv(alphadiv_fung_mouse_lung, "Fig S10C Lung ITS Alpha Diversity.csv")


###Beta Diversity###

exp2_fung_mouse_lung <- dada2_to_phyloseq(fung_mouse_lung,fung.mouse.lung.taxa,roughmetadata_fung_mouse_lung, read_count = 0, rough = FALSE)
exp2_fung_mouse_rel_prev <- transform_sample_counts(exp2_fung_mouse_lung, function(x) x / sum(x) )

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
tiff("Fig S10C Mouse Lung ITS PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
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



#############Fig. S10D: Lung 16S#############

###Data Prep###
bact_mouse_lung <- read_csv("Mouse_Lung_16s_ASVs.csv")
bact.mouse.lung.taxa <- read_csv("Mouse_Lung_16s_Taxa.csv")
roughmetadata_bact_mouse_lung <- read_csv("Mouse_Lung_Metadata.csv")


###Alpha Diversity###
exp2_bact_mouse_lung_rough <- dada2_to_phyloseq(bact_mouse_lung,bact.mouse.lung.taxa,roughmetadata_bact_mouse_lung, read_count = 0, rough = TRUE)
alphadiv_bact_mouse_lung <- alphadiv(exp2_bact_mouse_lung_rough, "Condition")

alphadiv_bact_mouse_lung <- alphadiv_bact_mouse_lung %>%
  mutate(
    Organ = exp2_bact_mouse_lung_rough@sam_data[["Organ"]],
    Oxygen = exp2_bact_mouse_lung_rough@sam_data[["Oxygen"]],
    FMT = exp2_bact_mouse_lung_rough@sam_data[["FMT"]],
    Condition = exp2_bact_mouse_lung_rough@sam_data[["Condition"]]
  )

#Save


###Beta Diversity###

exp2_bact_mouse_lung <- dada2_to_phyloseq(bact_mouse_lung,bact.mouse.lung.taxa,roughmetadata_bact_mouse_lung, read_count = 0, rough = FALSE)
exp2_bact_mouse_rel_prev <- transform_sample_counts(exp2_bact_mouse_lung, function(x) x / sum(x) )

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
tiff("Fig S10D Mouse lung 16s PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
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



