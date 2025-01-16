library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(factoextra)
library(SpiecEasi)
library(igraph)
library(plyr)

############################FIGURE E4A: SAAB ITS ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
#Get Chao1 and Shannon diversity
richness_est_fung <- estimate_richness(fr_fung, measures = c("Chao1", "Shannon"))

#Wilcoxon rank-sum test for significance
wilcox_alpha_fung <- t(sapply(richness_est_fung, function(x) unlist(wilcox.test(x~sample_data(fr_fung)$Sex)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_fung

#Create a data frame with alpha diversity and relevant metadata values for use in GraphPad
richness_est_fung <- richness_est_fung %>%
  mutate(
    BPD = md.fung$BPD,
    Sex = md.fung$Sex,
    BPD_Sex = md.fung$BPD_Sex,
    Oxygen = md.fung$Oxygen_Days
  )
#Save as .csv and use in GraphPad


############################FIGURE E4B: SAAB 16S ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_fung <- prune_taxa(taxa_sums(exp2_combined_rough) > 0, exp2_combined_rough)
#Get Chao1 and Shannon diversity
richness_est_combined <- estimate_richness(fr_combined, measures = c("Chao1", "Shannon"))

#Wilcoxon rank-sum test for significance
wilcox_alpha_combined <- t(sapply(richness_est_combined, function(x) unlist(wilcox.test(x~sample_data(fr_combined)$Sex)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_combined

#Create a data frame with alpha diversity and relevant metadata values for use in GraphPad
richness_est_combined <- richness_est_combined %>%
  mutate(
    BPD = md.combined$BPD,
    Sex = md.combined$Sex,
    BPD_Sex = md.combined$BPD_Sex,
    Oxygen = md.combined$Oxygen_Days
  )
#Save as .csv and use in GraphPad

############################FIGURE E4C: SAAB ITS BETA DIVERSITY PCoA############################

#Convert to relative abundance instead of counts
exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)

#Get Bray-Curtis dissimilarity and run PCoA
exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
exp2_fung_rel_pcoa$values

#Colors points according to BPD status
factor_exp2 <- as.factor(exp2_fung_meta_rel$Sex)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Beta Diversity ITS PCoA SAAB.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1(9.44% Explained)", ylab="PC2(7.31% Explained)", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp2, col = pca_colors)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ Sex, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$Sex)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(permdisp_fung_pcoa)


############################FIGURE E4D: SAAB 16S BETA DIVERSITY PCoA############################

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)

#Ordinate PCoA by Bray-Curtis dissimilarity
exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
exp2_combined_rel_pcoa$values

#Set variable of interest and color palette
factor_exp2 <- as.factor(exp2_combined_meta_rel$Sex)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Plot!
dpi=600
tiff("Beta Diversity 16S PCoA SAAB.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.4), c(-0.4, 0.5), font = 2, font.lab = 2, xlab="PC1(5.41% Explained)", ylab="PC2(4.35% Explained)", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2, kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df <- data.frame(exp2_combined_meta_rel)
set.seed(1312)
permanova_combined_pcoa <- vegan::adonis2(exp2_combined_otu_rel ~ Sex, data = permanova_pcoa_df, method="bray", permutations = 10000)

#PERMDISP
combined_pcoa_dist <- vegdist(exp2_combined_otu_rel, method = "bray")
disp_combined_pcoa <- betadisper(combined_pcoa_dist, permanova_pcoa_df$Sex)
set.seed(1312)
permdisp_combined_pcoa <- permutest(disp_combined_pcoa, permutations = 10000)

print(permanova_combined_pcoa)
print(permdisp_combined_pcoa)

############################FIGURE E4E: SAAB 16S DESeq2 ASV Boxplot############################

diff = phyloseq_to_deseq2(exp2_combined_gen_prev, ~ Sex)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")
deseq2results = results(diff, pAdjustMethod = "fdr")
deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp2_combined_gen_prev)[rownames(deseq2results), ], "matrix"))
sigtab = deseq2results[(deseq2results$padj < 0.05), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp2_combined_gen_prev)[rownames(sigtab), ], "matrix"))
sigtabgenus = subset(sigtab, !is.na(Genus))
print(sigtabgenus)
#Save table as .csv

exp2_combined_rel_gen_prev <- transform_sample_counts(exp2_combined_gen_prev, function(x) x / sum(x) )

exp2_combined_signif_genus <- subset_taxa(exp2_combined_rel_gen_prev, rownames(otu_table(exp2_combined_rel_gen_prev)) %in% rownames(sigtabgenus))
exp2_combined_signif_otu_genus <- as.data.frame(t(exp2_combined_signif_genus@otu_table))
exp2_combined_signif_meta_genus <- as.data.frame(exp2_combined_signif_genus@sam_data)
colnames(exp2_combined_signif_otu_genus) <- c("ASV7: Bifidobacterium","ASV18049: Candida")
exp2_combined_diffabund_signif_genus <- exp2_combined_signif_otu_genus %>%
  dplyr::mutate(Sex = exp2_combined_signif_meta_genus$Sex)
#Save table as .csv for use in GraphPad

############################FIGURE E4F: SAAB 16S DESeq2 Plot############################

deseq2_colors <- c("#007016","#94D57F")
deseq2plot <- function(sigtabgenus){
  diffabund <- ggplot(sigtabgenus,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
    geom_col() +
    coord_flip() +
    scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
  diffabund2 <- diffabund + theme_bw()
  diffabund2
}
deseq2plot(sigtabgenus)
#Save

