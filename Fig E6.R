library(tidyverse)
library(phyloseq)
library(vegan)

############################FIGURE E6A: ITS PCOA BY OXYGEN CONCENTRATION############################

exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)

###Ordinate using Bray-Curtis dissimilarity###
exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
exp2_fung_rel_pcoa$values

###Set variable of interest and gradient color palette###
factor_exp2_ga <- ordered(exp2_fung_meta_rel$Oxygen_Weeks)
type_exp2_ga <- as.numeric(factor_exp2_ga)
pcoa_pal_ga <- viridis(19, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")


###Plot PCoA###
dpi=600
tiff("Oxygen Weeks PCoA ITS With CI.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.5, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pcoa_pal_ga[type_exp2_ga], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95,factor_exp2, label = TRUE,col = pca_colors)
dev.off()

###Create legend###
tiff("Oxygen Weeks ITS Legend.tif", width=5*dpi, height=5*dpi, res=dpi)
legend_image <- as.raster(matrix(rev(pcoa_pal_ga), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Weeks on Oxygen')
text(x=1.5, y = seq(0,1,l=5), labels = seq(1,19,l=5))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

###PERMANOVA, omit any NA value###
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_oxygen_weeks <- vegan::adonis2(exp2_fung_otu_rel ~ Oxygen_Weeks, data = permanova_pcoa_df, method="bray", na.action = na.omit,permutations = 10000)
print(permanova_fung_oxygen_weeks)

###Pairwise PERMANOVA, remove NA values from dataframe first###
permanova_pcoa_df_noNA <-permanova_pcoa_df[-c(1,18,41),]
exp2_fung_otu_rel_noNA <-exp2_fung_otu_rel[-c(1,18,41),]

pairwise_permanova_fung_oxygen_weeks <- permanova_pairwise(exp2_fung_otu_rel_noNA, grp = permanova_pcoa_df_noNA$Oxygen_Weeks, permutations = 10000, method = "bray", padj = "fdr")

###PERMDISP###
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Oxygen_Weeks)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_oxygen_weeks)
print(pairwise_permanova_fung_oxygen_weeks)
print(permdisp_fung_pcoa)


############################FIGURE E6B: MULTIKINGDOM PCOA BY OXYGEN CONCENTRATION############################

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)

###Ordinate using Bray-Curtis dissimilarity###
exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
exp2_combined_rel_pcoa$values

###Set variable of interest and gradient color palette###
factor_exp2_ga <- ordered(exp2_combined_meta_rel$Oxygen_Weeks)
type_exp2_ga <- as.numeric(factor_exp2_ga)
pcoa_pal_ga <- viridis(19, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")

###Plot PCoA, using same legend as ITS so no need to plot that twice###
dpi=600
tiff("Oxygen Weeks PCoA Multikingdom.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.4, 0.5), font = 2, font.lab = 2, xlab="PC1(5.41% Explained)", ylab="PC2(4.35% Explained)", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pcoa_pal_ga[type_exp2_ga], lwd = 1)
dev.off()

###PERMANOVA###
permanova_pcoa_df <- data.frame(exp2_combined_meta_rel)
set.seed(1312)
permanova_combined_oxygen_weeks <- vegan::adonis2(exp2_combined_otu_rel ~ Oxygen_Weeks, data = permanova_pcoa_df, method="bray", na.action = na.omit,permutations = 10000)
print(permanova_combined_oxygen_weeks)
oxygen_weeks_p <- permanova_combined_oxygen_weeks$`Pr(>F)`

###Pairwise PERMANOVA###
permanova_pcoa_df_noNA <-permanova_pcoa_df[-c(1,18,41),]
exp2_combined_otu_rel_noNA <-exp2_combined_otu_rel[-c(1,18,41),]

pairwise_permanova_combined_oxygen_weeks <- permanova_pairwise(exp2_combined_otu_rel_noNA, grp = permanova_pcoa_df_noNA$Oxygen_Weeks, permutations = 10000, method = "bray", padj = "fdr")

###PERMDISP###
combined_pcoa_dist <- vegdist(exp2_combined_otu_rel, method = "bray")
disp_combined_pcoa <- betadisper(combined_pcoa_dist, permanova_pcoa_df$Oxygen_Weeks)
set.seed(1312)
permdisp_combined_pcoa <- permutest(disp_combined_pcoa, permutations = 10000)

print(permanova_combined_oxygen_weeks)
print(pairwise_permanova_combined_oxygen_weeks)
print(permdisp_combined_pcoa)

