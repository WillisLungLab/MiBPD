library(tidyverse)
library(vegan)
library(phyloseq)
library(ggpubr)
library(DESeq2)
library(indicspecies)
library(data.table)
library(RColorBrewer)
library(factoextra)
library(SpiecEasi)
library(igraph)
library(plyr)
library(mia)
library(miaViz)

############### ALL SAMPLES DMM ITS (FIG. 3A-C) ###############

exp2_fung_tse <- makeTreeSEFromPhyloseq(exp2_fung_prev)

set.seed(1312)
tse_dmn <- mia::runDMN(exp2_fung_tse, name = "DMN", k = 1:7)
getDMN(tse_dmn)

plotDMNFit(tse_dmn, type = "laplace")
best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")

#######

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
fitted <- head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))

prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
# Add column names
colnames(prob) <- c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4", "Cluster 5")

#For each row, finds column that has the highest value. Then extract the column 
#names of highest values.
vec <- colnames(prob)[max.col(prob,ties.method = "first")]

#Convert to relative abundance, add corresponding component to metadata, and run PCoA using Bray-Curtis dissimilarity
exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)
exp2_fung_meta_rel <- cbind(exp2_fung_meta_rel,cluster = vec)

exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
bray_pcoa_df <- exp2_fung_rel_pcoa$vectors[,1:2]

bray_dmm_pcoa_df <- cbind(bray_pcoa_df,
                          cluster = vec)

### ITS PCoA BY DMM CLUSTER (Fig. 3A) ###
factor_exp2 <- as.factor(bray_dmm_pcoa_df[,3])
type_exp2 <- as.numeric(factor_exp2)
dmm_colors_all <- c("#8FD96C","#64BD4B","#409E2F","#228019","#096620")

dpi=600
tiff("Figure_3A.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.5), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21,cex = 1.3, bg = dmm_colors_all[type_exp2], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ cluster, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Pairwise PERMANOVA
pairwise_permanova_fung_pcoa <- permanova_pairwise(exp2_fung_otu_rel, grp = permanova_pcoa_df$cluster, permutations = 10000, method = "bray", padj = "fdr")

#PERMDISP
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$cluster)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(pairwise_permanova_fung_pcoa)
print(permdisp_fung_pcoa)

#Adds component to metadata of the main phyloseq object
sample_data(exp2_fung_prev)$cluster <- vec
sample_data(exp2_fung_gen_prev)$cluster <- vec

### ITS ALPHA DIVERSITY BY DMM CLUSTER (Fig. 3B) ###

fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
richness_est_fung <- estimate_richness(fr_fung, measures = c("Simpson", "Shannon"))


vec_df <- as.data.frame(vec)
row.names(vec_df) <- row.names(exp2_fung_meta_rel)
richness_est_fung_cluster <- richness_est_fung[which(unlist(row.names(richness_est_fung)) %in% row.names(vec_df)),]

wilcox_alpha_fung_cluster <- t(sapply(richness_est_fung_cluster, function(x) unlist(kruskal.test(x~vec_df$vec)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_fung_cluster

richness_est_fung_cluster <- cbind(richness_est_fung_cluster, vec_df)
#Save .csv for GraphPad

shannon_prelim <- ggboxplot(richness_est_fung_cluster, x = "vec", y = "Shannon", add = "jitter", color = "vec")
shannon_prelim + stat_compare_means()
#Save


simpson_prelim <- ggboxplot(richness_est_fung_cluster, x = "vec", y = "Simpson", add = "jitter", color = "vec")
simpson_prelim + stat_compare_means()
#Save



### ITS BARPLOT BY DMM CLUSTER (Fig. 3C) ###

exp2_fung_merged = merge_samples(exp2_fung_prev, "cluster")
exp2_fung_gen_merged <- tax_glom(exp2_fung_merged, taxrank = 'Genus')

top20_fung_prev_list <- names(sort(taxa_sums(exp2_fung_gen_merged), decreasing=TRUE)[1:19])
top20_prev_fung_rel <- transform_sample_counts(exp2_fung_gen_merged, function(x) x / sum(x) )
top20_prev_fung_df <- psmelt(top20_prev_fung_rel)
top20_prev_fung_df[!(top20_prev_fung_df$OTU %in% top20_fung_prev_list),]$Genus <- 'Other'


###Barplot of top 20 combined genera, and an entry for all others###

barplot_gen_bpd_its <- ggplot(top20_prev_fung_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_its + scale_fill_manual(values = barplot_colors)
#Save


############### ITS DMM NoBPD ONLY (FIG. 3D-F) ###############
exp2_fung_tse_pprd <- makeTreeSEFromPhyloseq(exp2_fung_prev_pprd)

set.seed(1312)
tse_dmn_pprd <- mia::runDMN(exp2_fung_tse_pprd, name = "DMN", k = 1:7)
getDMN(tse_dmn_pprd)

plotDMNFit(tse_dmn_pprd, type = "laplace")
best_dmn_pprd <- getBestDMNFit(tse_dmn_pprd, type = "laplace")

######

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn_pprd))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_pprd)))
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn_pprd)))

prob_pprd <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_pprd))
# Add column names
colnames(prob_pprd) <- c("Cluster 1", "Cluster 2")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec_pprd <- colnames(prob_pprd)[max.col(prob_pprd,ties.method = "first")]

#Set up PCoA
exp2_fung_rel_prev_pprd <- transform_sample_counts(exp2_fung_prev_pprd, function(x) x / sum(x) )

exp2_fung_otu_rel_pprd <- as.data.frame(t(exp2_fung_rel_prev_pprd@otu_table))
exp2_fung_tax_rel_pprd <- as.data.frame(exp2_fung_rel_prev_pprd@tax_table)
exp2_fung_meta_rel_pprd <- as.data.frame(exp2_fung_rel_prev_pprd@sam_data)
exp2_fung_meta_rel_pprd <- cbind(exp2_fung_meta_rel_pprd,
                                 cluster = vec_pprd)

exp2_fung_rel_bray_pprd = vegdist(exp2_fung_otu_rel_pprd, method='bray')
exp2_fung_rel_pcoa_pprd <- ape::pcoa(exp2_fung_rel_bray_pprd)
bray_pcoa_df_pprd <- exp2_fung_rel_pcoa_pprd$vectors[,1:2]

bray_dmm_pcoa_df_pprd <- cbind(bray_pcoa_df_pprd,
                               cluster = vec_pprd)

### ITS PCoA BY DMM CLUSTER (Fig. 3D) ###
factor_exp2_pprd <- as.factor(bray_dmm_pcoa_df_pprd[,3])
type_exp2_pprd <- as.numeric(factor_exp2_pprd)
dmm_colors_pprd <- c("#8FD96C","#228019")

dpi=600
tiff("Figure_3D.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.5), c(-0.4, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa_pprd$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_pprd[type_exp2_pprd], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa_pprd$vectors[,1:2], factor_exp2_pprd)
ordispider(exp2_fung_rel_pcoa_pprd$vectors[,1:2], factor_exp2_pprd, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df_pprd <- data.frame(exp2_fung_meta_rel_pprd)
set.seed(1312)
permanova_fung_pcoa_pprd <- vegan::adonis2(exp2_fung_otu_rel_pprd ~ cluster, data = permanova_pcoa_df_pprd, method="bray", permutations = 10000)

#PERMDISP
fung_pcoa_dist_pprd <- vegdist(exp2_fung_otu_rel_pprd, method = "bray")
disp_fung_pcoa_pprd <- betadisper(fung_pcoa_dist_pprd, permanova_pcoa_df_pprd$cluster)
set.seed(1312)
permdisp_fung_pcoa_pprd <- permutest(disp_fung_pcoa_pprd, permutations = 10000)

print(permanova_fung_pcoa_pprd)
print(permdisp_fung_pcoa_pprd)

sample_data(exp2_fung_prev_pprd)$cluster <- vec_pprd
sample_data(exp2_fung_gen_prev_pprd)$cluster <- vec_pprd

### ITS ALPHA DIVERSITY BY DMM CLUSTER (Fig. 3E) ###

vec_pprd_df <- as.data.frame(vec_pprd)
row.names(vec_pprd_df) <- row.names(exp2_fung_meta_rel_pprd)
richness_est_fung_cluster_pprd <- richness_est_fung[which(unlist(row.names(richness_est_fung)) %in% row.names(vec_pprd_df)),]

wilcox_alpha_bact <- t(sapply(richness_est_fung_cluster_pprd, function(x) unlist(kruskal.test(x~vec_pprd_df$vec)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_bact

richness_est_fung_cluster_pprd <- cbind(richness_est_fung_cluster_pprd, vec_pprd_df)
#Save as .csv

shannon_prelim_pprd <- ggboxplot(richness_est_fung_cluster_pprd, x = "vec_pprd", y = "Shannon", add = "jitter", color = "vec_pprd")
shannon_prelim_pprd + stat_compare_means()
#Save


simpson_prelim <- ggboxplot(richness_est_fung_cluster_pprd, x = "vec_pprd", y = "Simpson", add = "jitter", color = "vec_pprd")
simpson_prelim + stat_compare_means()
#Save


### ITS TAXA BARPLOT BY DMM CLUSTER (Fig. 3F) ###

exp2_fung_merged_pprd = merge_samples(exp2_fung_prev_pprd, "cluster")
exp2_fung_gen_merged_pprd <- tax_glom(exp2_fung_merged_pprd, taxrank = 'Genus')

top20_fung_prev_list_pprd <- names(sort(taxa_sums(exp2_fung_gen_merged_pprd), decreasing=TRUE)[1:19])
top20_prev_fung_rel_pprd <- transform_sample_counts(exp2_fung_gen_merged_pprd, function(x) x / sum(x) )
top20_prev_fung_df_pprd <- psmelt(top20_prev_fung_rel_pprd)
top20_prev_fung_df_pprd[!(top20_prev_fung_df_pprd$OTU %in% top20_fung_prev_list_pprd),]$Genus <- 'Other'


#Barplot of top 20 combined genera, and an entry for all others

barplot_gen_its_pprd <- ggplot(top20_prev_fung_df_pprd, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_its_pprd + scale_fill_manual(values = barplot_colors)
#Save
#Save



############### ITS DMM BPD ONLY (FIG. 3G-I) ###############

exp2_fung_tse_bpd <- makeTreeSEFromPhyloseq(exp2_fung_prev_bpd)

set.seed(1312)
tse_dmn_bpd <- mia::runDMN(exp2_fung_tse_bpd, name = "DMN", k = 1:7)
getDMN(tse_dmn_bpd)

plotDMNFit(tse_dmn_bpd, type = "laplace")
best_dmn_bpd <- getBestDMNFit(tse_dmn_bpd, type = "laplace")

######

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn_bpd))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_bpd)))
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn_bpd)))

prob_bpd <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_bpd))
# Add column names
colnames(prob_bpd) <- c("Cluster 1", "Cluster 2")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec_bpd <- colnames(prob_bpd)[max.col(prob_bpd,ties.method = "first")]

#Set up PCoA
exp2_fung_rel_prev_bpd <- transform_sample_counts(exp2_fung_prev_bpd, function(x) x / sum(x) )

exp2_fung_otu_rel_bpd <- as.data.frame(t(exp2_fung_rel_prev_bpd@otu_table))
exp2_fung_tax_rel_bpd <- as.data.frame(exp2_fung_rel_prev_bpd@tax_table)
exp2_fung_meta_rel_bpd <- as.data.frame(exp2_fung_rel_prev_bpd@sam_data)
exp2_fung_meta_rel_bpd <- cbind(exp2_fung_meta_rel_bpd,
                                cluster = vec_bpd)

### ITS PCoA BY DMM CLUSTER (Fig. 3G) ###
exp2_fung_rel_bray_bpd = vegdist(exp2_fung_otu_rel_bpd, method='bray')
exp2_fung_rel_pcoa_bpd <- ape::pcoa(exp2_fung_rel_bray_bpd)
bray_pcoa_df_bpd <- exp2_fung_rel_pcoa_bpd$vectors[,1:2]

bray_dmm_pcoa_df_bpd <- cbind(bray_pcoa_df_bpd,
                              cluster = vec_bpd)

factor_exp2_bpd <- as.factor(bray_dmm_pcoa_df_bpd[,3])
type_exp2_bpd <- as.numeric(factor_exp2_bpd)
dmm_colors_bpd <- c("#8FD96C","#228019")

dpi=600
tiff("Figure_3G.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa_bpd$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_bpd[type_exp2_bpd], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa_bpd$vectors[,1:2], factor_exp2_bpd)
ordispider(exp2_fung_rel_pcoa_bpd$vectors[,1:2], factor_exp2_bpd, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df_bpd <- data.frame(exp2_fung_meta_rel_bpd)
set.seed(1312)
permanova_fung_pcoa_bpd <- vegan::adonis2(exp2_fung_otu_rel_bpd ~ cluster, data = permanova_pcoa_df_bpd, method="bray", permutations = 10000)

#Pairwise PERMANOVA
pairwise_permanova_fung_pcoa_bpd <- permanova_pairwise(exp2_fung_otu_rel_bpd, grp = permanova_pcoa_df_bpd$cluster, permutations = 10000, method = "bray", padj = "fdr")

#PERMDISP
fung_pcoa_dist_bpd <- vegdist(exp2_fung_otu_rel_bpd, method = "bray")
disp_fung_pcoa_bpd <- betadisper(fung_pcoa_dist_bpd, permanova_pcoa_df_bpd$cluster)
set.seed(1312)
permdisp_fung_pcoa_bpd <- permutest(disp_fung_pcoa_bpd, permutations = 10000)

print(permanova_fung_pcoa_bpd)
print(pairwise_permanova_fung_pcoa_bpd)
print(permdisp_fung_pcoa_bpd)

sample_data(exp2_fung_prev_bpd)$cluster <- vec_bpd
sample_data(exp2_fung_gen_prev_bpd)$cluster <- vec_bpd

### ITS ALPHA DIVERSITY BY DMM CLUSTER (Fig. 3H) ###

vec_bpd_df <- as.data.frame(vec_bpd)
row.names(vec_bpd_df) <- row.names(exp2_fung_meta_rel_bpd)
richness_est_fung_cluster_bpd <- richness_est_fung[which(unlist(row.names(richness_est_fung)) %in% row.names(vec_bpd_df)),]

wilcox_alpha_bact <- t(sapply(richness_est_fung_cluster_bpd, function(x) unlist(kruskal.test(x~vec_bpd_df$vec)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_bact

richness_est_fung_cluster_bpd <- cbind(richness_est_fung_cluster_bpd, vec_bpd_df)
#Save as .csv for GraphPad


shannon_prelim_bpd <- ggboxplot(richness_est_fung_cluster_bpd, x = "vec_bpd", y = "Shannon", add = "jitter", color = "vec_bpd")
shannon_prelim_bpd + stat_compare_means()
#Save


simpson_prelim <- ggboxplot(richness_est_fung_cluster_bpd, x = "vec_bpd", y = "Simpson", add = "jitter", color = "vec_bpd")
simpson_prelim + stat_compare_means()
#Save


### ITS TAXA BARPLOT BY DMM CLUSTER (Fig. 3I) ###

exp2_fung_merged_bpd = merge_samples(exp2_fung_prev_bpd, "cluster")
exp2_fung_gen_merged_bpd <- tax_glom(exp2_fung_merged_bpd, taxrank = 'Genus')

top20_fung_prev_list_bpd <- names(sort(taxa_sums(exp2_fung_gen_merged_bpd), decreasing=TRUE)[1:19])
top20_prev_fung_rel_bpd <- transform_sample_counts(exp2_fung_gen_merged_bpd, function(x) x / sum(x) )
top20_prev_fung_df_bpd <- psmelt(top20_prev_fung_rel_bpd)
top20_prev_fung_df_bpd[!(top20_prev_fung_df_bpd$OTU %in% top20_fung_prev_list_bpd),]$Genus <- 'Other'


###Barplot of top 20 combined genera, and an entry for all others###

barplot_gen_bpd_its_bpd <- ggplot(top20_prev_fung_df_bpd, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_its_bpd + scale_fill_manual(values = barplot_colors)
#Save

