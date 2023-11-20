library(tidyverse)
library(readxl)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(factoextra)
library(plyr)


############################FIGURE S3A: Multikingdom Alpha Diversity############################

fr_combined <- prune_taxa(taxa_sums(exp2_combined_rough) > 0, exp2_combined_rough)
richness_est_combined <- estimate_richness(fr_combined, measures = c("Simpson", "Shannon"))

#Wilcox test for significance
wilcox_alpha_combined <- t(sapply(richness_est_combined, function(x) unlist(wilcox.test(x~sample_data(fr_combined)$BPD)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined

#Get dataframe with alpha diversity values and metadata values, to plot in GraphPad
richness_est_combined <- richness_est_combined %>%
  mutate(
    BPD = md.combined$BPD,
    Sex = md.combined$Sex,
    BPD_Sex = md.combined$BPD_Sex,
    Oxygen = md.combined$Oxygen_Days
  )
#Save as .csv


############################FIGURE S3B: Multikingdom Beta Diversity############################

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)

#PCoA using Bray-Curtis dissimilarity
exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
exp2_combined_rel_pcoa$values

factor_exp2 <- as.factor(exp2_combined_meta_rel$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Spider plot with confidence interval
dpi=600
tiff("Fig E3B.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.4), c(-0.4, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp2, col = pca_colors)
ordispider(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df <- data.frame(exp2_combined_meta_rel)
set.seed(1312)
permanova_combined_pcoa <- vegan::adonis2(exp2_combined_otu_rel ~ BPD, data = permanova_pcoa_df, method="bray", permutations = 10000)

#PERMDISP
combined_pcoa_dist <- vegdist(exp2_combined_otu_rel, method = "bray")
disp_combined_pcoa <- betadisper(combined_pcoa_dist, permanova_pcoa_df$BPD)
set.seed(1312)
permdisp_combined_pcoa <- permutest(disp_combined_pcoa, permutations = 10000)

print(permanova_combined_pcoa)
print(permdisp_combined_pcoa)
