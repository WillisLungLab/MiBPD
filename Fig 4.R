library(tidyverse)
library(vegan)
library(phyloseq)
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

###############ALL SAMPLES###############
exp2_fung_tse <- makeTreeSEFromPhyloseq(exp2_fung_prev)

set.seed(1312)
tse_dmn <- mia::runDMN(exp2_fung_tse, name = "DMN", k = 1:7)
getDMN(tse_dmn)

plotDMNFit(tse_dmn, type = "laplace")
best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")


###############ALL SAMPLES, 2 COMPONENTS (FIG. 4A-C)###############
tse_dmn2 <- mia::runDMN(exp2_fung_tse, name = "DMN", k = 2)
getDMN(tse_dmn2)

plotDMNFit(tse_dmn2, type = "laplace")
best_dmn2 <- getBestDMNFit(tse_dmn, type = "laplace")

######

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn2))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn2)))
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn2)))

prob2 <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn2))
# Add column names
colnames(prob2) <- c("comp1", "comp2")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec2 <- colnames(prob2)[max.col(prob2,ties.method = "first")]

#Set up PCoA
exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel2 <- as.data.frame(exp2_fung_rel_prev@sam_data)
exp2_fung_meta_rel2 <- cbind(exp2_fung_meta_rel2,dmm_component = vec2)

exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
bray_pcoa_df <- exp2_fung_rel_pcoa$vectors[,1:2]

bray_dmm_pcoa_df2 <- cbind(bray_pcoa_df,
                           dmm_component = vec2)

###Plot PCoA (Fig. 4A)###
factor2_exp2 <- as.factor(bray_dmm_pcoa_df2[,3])
type2_exp2 <- as.numeric(factor2_exp2)
dmm_colors_all2 <- c("#8FD96C","#228019")

dpi=600
tiff("ITS all DMM 2 Clusters PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1(9.44% Explained)", ylab="PC2(7.31% Explained)", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_all2[type2_exp2], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], factor2_exp2)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor2_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df2 <- data.frame(exp2_fung_meta_rel2)
set.seed(1312)
permanova_fung_pcoa2 <- vegan::adonis2(exp2_fung_otu_rel2 ~ dmm_component, data = permanova_pcoa_df2, method="bray", permutations = 10000)

#PERMDISP
fung_pcoa_dist2 <- vegdist(exp2_fung_otu_rel2, method = "bray")
disp_fung_pcoa2 <- betadisper(fung_pcoa_dist2, permanova_pcoa_df2$dmm_component)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa2)
print(permdisp_fung_pcoa2)

#Add component to main phyloseq object
exp2_fung_prev2 <- exp2_fung_prev
sample_data(exp2_fung_prev2)$cluster <- vec2
exp2_fung_gen_prev2 <- exp2_fung_gen_prev
sample_data(exp2_fung_gen_prev)$cluster <- vec2

###ALPHA DIVERSITY (Fig. 4B)###

fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
richness_est_fung <- estimate_richness(fr_fung, measures = c("Chao1", "Shannon"))


vec_df2 <- as.data.frame(vec2)
row.names(vec_df2) <- row.names(exp2_fung_meta_rel)
richness_est_fung_cluster2 <- richness_est_fung[which(unlist(row.names(richness_est_fung)) %in% row.names(vec_df2)),]

richness_est_fung_cluster2 <- cbind(richness_est_fung_cluster2, vec_df2)
#Save .csv for GraphPad

###BARPLOT (Fig. 4C)###

exp2_fung_merged2 = merge_samples(exp2_fung_prev2, "cluster")
exp2_fung_gen_merged2 <- tax_glom(exp2_fung_merged2, taxrank = 'Genus')

top20_fung_prev_list2 <- names(sort(taxa_sums(exp2_fung_gen_merged2), decreasing=TRUE)[1:19])
top20_prev_fung_rel2 <- transform_sample_counts(exp2_fung_gen_merged2, function(x) x / sum(x) )
top20_prev_fung_df2 <- psmelt(top20_prev_fung_rel2)
top20_prev_fung_df2[!(top20_prev_fung_df2$OTU %in% top20_fung_prev_list2),]$Genus <- 'Other'


###Barplot of top 20 combined genera, and an entry for all others###

barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(20)

barplot_gen_its_2 <- ggplot(top20_prev_fung_df2, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_its_2 + scale_fill_manual(values = barplot_colors)
#Save


###############ALL SAMPLES, 5 COMPONENTS (FIG. 4D-F)###############

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
fitted <- head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))

prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
# Add column names
colnames(prob) <- c("comp1", "comp2","comp3", "comp4", "comp5")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec <- colnames(prob)[max.col(prob,ties.method = "first")]

#Convert to relative abundance, add corresponding component to metadata, and run PCoA using Bray-Curtis dissimilarity
exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)
exp2_fung_meta_rel <- cbind(exp2_fung_meta_rel,dmm_component = vec)

exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
bray_pcoa_df <- exp2_fung_rel_pcoa$vectors[,1:2]

bray_dmm_pcoa_df <- cbind(bray_pcoa_df,
                          dmm_component = vec)

###Plot PCoA (Fig. 4D)###
factor_exp2 <- as.factor(bray_dmm_pcoa_df[,3])
type_exp2 <- as.numeric(factor_exp2)
dmm_colors_all <- c("#8FD96C","#64BD4B","#409E2F","#228019","#096620")

dpi=600
tiff("ITS all DMM 5 Clusters PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1(9.44% Explained)", ylab="PC2(7.31% Explained)", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21,cex = 1.3, bg = dmm_colors_all[type_exp2], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ dmm_component, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Pairwise PERMANOVA
pairwise_permanova_fung_pcoa <- permanova_pairwise(exp2_fung_otu_rel, grp = permanova_pcoa_df$dmm_component, permutations = 10000, method = "bray", padj = "fdr")

#PERMDISP
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$dmm_component)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(pairwise_permanova_fung_pcoa)
print(permdisp_fung_pcoa)

#Adds component to metadata of the main phyloseq object
sample_data(exp2_fung_prev)$cluster <- vec
sample_data(exp2_fung_gen_prev)$cluster <- vec

###ALPHA DIVERSITY (Fig. 4E)###

fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
richness_est_fung <- estimate_richness(fr_fung, measures = c("Chao1", "Shannon"))


vec_df <- as.data.frame(vec)
row.names(vec_df) <- row.names(exp2_fung_meta_rel)
richness_est_fung_cluster <- richness_est_fung[which(unlist(row.names(richness_est_fung)) %in% row.names(vec_df)),]

richness_est_fung_cluster <- cbind(richness_est_fung_cluster, vec_df)
#Save .csv for GraphPad

###BARPLOT (Fig. 4F)###

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


###############BPD ONLY (FIG. 4G-I)###############

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
colnames(prob_bpd) <- c("comp1", "comp2","comp3")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec_bpd <- colnames(prob_bpd)[max.col(prob_bpd,ties.method = "first")]

#Set up PCoA
exp2_fung_rel_prev_bpd <- transform_sample_counts(exp2_fung_prev_bpd, function(x) x / sum(x) )

exp2_fung_otu_rel_bpd <- as.data.frame(t(exp2_fung_rel_prev_bpd@otu_table))
exp2_fung_tax_rel_bpd <- as.data.frame(exp2_fung_rel_prev_bpd@tax_table)
exp2_fung_meta_rel_bpd <- as.data.frame(exp2_fung_rel_prev_bpd@sam_data)
exp2_fung_meta_rel_bpd <- cbind(exp2_fung_meta_rel_bpd,
                                dmm_component = vec_bpd)

###Plot PCoA (Fig. 4G)###
exp2_fung_rel_bray_bpd = vegdist(exp2_fung_otu_rel_bpd, method='bray')
exp2_fung_rel_pcoa_bpd <- ape::pcoa(exp2_fung_rel_bray_bpd)
bray_pcoa_df_bpd <- exp2_fung_rel_pcoa_bpd$vectors[,1:2]

bray_dmm_pcoa_df_bpd <- cbind(bray_pcoa_df_bpd,
                              dmm_component = vec_bpd)

factor_exp2_bpd <- as.factor(bray_dmm_pcoa_df_bpd[,3])
type_exp2_bpd <- as.numeric(factor_exp2_bpd)
dmm_colors_bpd <- c("#8FD96C","#409E2F","#096620")

dpi=600
tiff("ITS all DMM 5 Clusters PCoA BPD Only.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1(19.4% Explained)", ylab="PC2(13.9% Explained)", type="n")
points(exp2_fung_rel_pcoa_bpd$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_bpd[type_exp2_bpd], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa_bpd$vectors[,1:2], factor_exp2_bpd)
ordispider(exp2_fung_rel_pcoa_bpd$vectors[,1:2], factor_exp2_bpd, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df_bpd <- data.frame(exp2_fung_meta_rel_bpd)
set.seed(1312)
permanova_fung_pcoa_bpd <- vegan::adonis2(exp2_fung_otu_rel_bpd ~ dmm_component, data = permanova_pcoa_df_bpd, method="bray", permutations = 10000)

#Pairwise PERMANOVA
pairwise_permanova_fung_pcoa_bpd <- permanova_pairwise(exp2_fung_otu_rel_bpd, grp = permanova_pcoa_df_bpd$dmm_component, permutations = 10000, method = "bray", padj = "fdr")

#PERMDISP
fung_pcoa_dist_bpd <- vegdist(exp2_fung_otu_rel_bpd, method = "bray")
disp_fung_pcoa_bpd <- betadisper(fung_pcoa_dist_bpd, permanova_pcoa_df_bpd$dmm_component)
set.seed(1312)
permdisp_fung_pcoa_bpd <- permutest(disp_fung_pcoa_bpd, permutations = 10000)

print(permanova_fung_pcoa_bpd)
print(pairwise_permanova_fung_pcoa_bpd)
print(permdisp_fung_pcoa_bpd)

sample_data(exp2_fung_prev_bpd)$cluster <- vec_bpd
sample_data(exp2_fung_gen_prev_bpd)$cluster <- vec_bpd

###ALPHA DIVERSITY (Fig. 4H)###

vec_bpd_df <- as.data.frame(vec_bpd)
row.names(vec_bpd_df) <- row.names(exp2_fung_meta_rel_bpd)
richness_est_fung_cluster_bpd <- richness_est_fung[which(unlist(row.names(richness_est_fung)) %in% row.names(vec_bpd_df)),]

richness_est_fung_cluster_bpd <- cbind(richness_est_fung_cluster_bpd, vec_bpd_df)
#Save as .csv

###BARPLOT (Fig. 4I)###

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

###############PPRD ONLY (FIG. 4J-L)###############
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
colnames(prob_pprd) <- c("comp1", "comp2","comp3")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec_pprd <- colnames(prob_pprd)[max.col(prob_pprd,ties.method = "first")]

#Set up PCoA
exp2_fung_rel_prev_pprd <- transform_sample_counts(exp2_fung_prev_pprd, function(x) x / sum(x) )

exp2_fung_otu_rel_pprd <- as.data.frame(t(exp2_fung_rel_prev_pprd@otu_table))
exp2_fung_tax_rel_pprd <- as.data.frame(exp2_fung_rel_prev_pprd@tax_table)
exp2_fung_meta_rel_pprd <- as.data.frame(exp2_fung_rel_prev_pprd@sam_data)
exp2_fung_meta_rel_pprd <- cbind(exp2_fung_meta_rel_pprd,
                                 dmm_component = vec_pprd)

exp2_fung_rel_bray_pprd = vegdist(exp2_fung_otu_rel_pprd, method='bray')
exp2_fung_rel_pcoa_pprd <- ape::pcoa(exp2_fung_rel_bray_pprd)
bray_pcoa_df_pprd <- exp2_fung_rel_pcoa_pprd$vectors[,1:2]

bray_dmm_pcoa_df_pprd <- cbind(bray_pcoa_df_pprd,
                               dmm_component = vec_pprd)

###Plot PCoA (Fig. 4J)###
factor_exp2_pprd <- as.factor(bray_dmm_pcoa_df_pprd[,3])
type_exp2_pprd <- as.numeric(factor_exp2_pprd)
dmm_colors_pprd <- c("#8FD96C","#409E2F","#096620")

dpi=600
tiff("ITS all DMM 5 Clusters PCoA PPRD Only.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.6, 0.4), font = 2, font.lab = 2, xlab="PC1(14.4% Explained)", ylab="PC2(10.7% Explained)", type="n")
points(exp2_fung_rel_pcoa_pprd$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_pprd[type_exp2_pprd], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa_pprd$vectors[,1:2], factor_exp2_pprd)
ordispider(exp2_fung_rel_pcoa_pprd$vectors[,1:2], factor_exp2_pprd, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df_pprd <- data.frame(exp2_fung_meta_rel_pprd)
set.seed(1312)
permanova_fung_pcoa_pprd <- vegan::adonis2(exp2_fung_otu_rel_pprd ~ dmm_component, data = permanova_pcoa_df_pprd, method="bray", permutations = 10000)

#Pairwise PERMANOVA
pairwise_permanova_fung_pcoa_pprd <- permanova_pairwise(exp2_fung_otu_rel_pprd, grp = permanova_pcoa_df_pprd$dmm_component, permutations = 10000, method = "bray", padj = "fdr")

#PERMDISP
fung_pcoa_dist_pprd <- vegdist(exp2_fung_otu_rel_pprd, method = "bray")
disp_fung_pcoa_pprd <- betadisper(fung_pcoa_dist_pprd, permanova_pcoa_df_pprd$dmm_component)
set.seed(1312)
permdisp_fung_pcoa_pprd <- permutest(disp_fung_pcoa_pprd, permutations = 10000)

print(permanova_fung_pcoa_pprd)
print(pairwise_permanova_fung_pcoa_pprd)
print(permdisp_fung_pcoa_pprd)

sample_data(exp2_fung_prev_pprd)$cluster <- vec_pprd
sample_data(exp2_fung_gen_prev_pprd)$cluster <- vec_pprd

###ALPHA DIVERSITY (Fig. 4K)###

vec_pprd_df <- as.data.frame(vec_pprd)
row.names(vec_pprd_df) <- row.names(exp2_fung_meta_rel_pprd)
richness_est_fung_cluster_pprd <- richness_est_fung[which(unlist(row.names(richness_est_fung)) %in% row.names(vec_pprd_df)),]

richness_est_fung_cluster_pprd <- cbind(richness_est_fung_cluster_pprd, vec_pprd_df)
#Save as .csv

###BARPLOT (Fig. 4L)###

exp2_fung_merged_pprd = merge_samples(exp2_fung_prev_pprd, "cluster")
exp2_fung_gen_merged_pprd <- tax_glom(exp2_fung_merged_pprd, taxrank = 'Genus')

top20_fung_prev_list_pprd <- names(sort(taxa_sums(exp2_fung_gen_merged_pprd), decreasing=TRUE)[1:19])
top20_prev_fung_rel_pprd <- transform_sample_counts(exp2_fung_gen_merged_pprd, function(x) x / sum(x) )
top20_prev_fung_df_pprd <- psmelt(top20_prev_fung_rel_pprd)
top20_prev_fung_df_pprd[!(top20_prev_fung_df_pprd$OTU %in% top20_fung_prev_list_pprd),]$Genus <- 'Other'


###Barplot of top 20 combined genera, and an entry for all others###

barplot_gen_pprd_its_pprd <- ggplot(top20_prev_fung_df_pprd, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_pprd_its_pprd + scale_fill_manual(values = barplot_colors)
#Save
