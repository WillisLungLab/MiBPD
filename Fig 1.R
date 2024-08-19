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
library(ggpubr)


############################FIGURE 1A: ITS ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
#Get Chao1 and Shannon diversity
richness_est_fung <- estimate_richness(fr_fung, measures = c("Simpson", "Shannon"))

#Wilcoxon rank-sum test for significance
wilcox_alpha_fung <- t(sapply(richness_est_fung, function(x) unlist(wilcox.test(x~sample_data(fr_fung)$BPD)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_fung

#Create a data frame with alpha diversity and relevant metadata values for use in GraphPad
richness_est_fung <- richness_est_fung %>%
  mutate(
    BPD = sample_data(fr_fung)$BPD,
    Sex = sample_data(fr_fung)$Sex,
    BPD_Sex = sample_data(fr_fung)$BPD_Sex,
    BPD_Severity = sample_data(fr_fung)$BPD_Severity,
    Oxygen = sample_data(fr_fung)$Oxygen
  )
#Save as .csv and use in GraphPad


############################FIGURE 1B: 16S ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_bact <- prune_taxa(taxa_sums(exp2_bact_rough) > 0, exp2_bact_rough)
#Get Chao1 and Shannon diversity
richness_est_bact <- estimate_richness(fr_bact, measures = c("Simpson", "Shannon"))

#Wilcoxon rank-sum test for significance
wilcox_alpha_bact <- t(sapply(richness_est_bact, function(x) unlist(wilcox.test(x~sample_data(fr_bact)$BPD)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_bact

#Create a data frame with alpha diversity and relevant metadata values for use in GraphPad
richness_est_bact <- richness_est_bact %>%
  mutate(
    BPD = sample_data(fr_bact)$BPD,
    Sex = sample_data(fr_bact)$Sex,
    BPD_Sex = sample_data(fr_bact)$BPD_Sex,
    BPD_Severity = sample_data(fr_bact)$BPD_Severity,
    Oxygen = sample_data(fr_bact)$Oxygen
  )
#Save as .csv and use in GraphPad

############################FIGURE 1C: ITS BETA DIVERSITY PCoA############################

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
factor_exp2 <- as.factor(exp2_fung_meta_rel$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Spider plot with the points and labeled centroids, with an ellipse representing the 95% confidence interval
dpi=600
tiff("Fig 1D.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.5), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp2, col = pca_colors)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ BPD, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$BPD)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(permdisp_fung_pcoa)

############################FIGURE 1D: ITS BETA DIVERSITY LOADING PLOT + PCA############################

exp2_fung_otu <- as.data.frame(t(exp2_fung_prev@otu_table))
exp2_fung_tax <- as.data.frame(exp2_fung_prev@tax_table)
exp2_fung_meta <- as.data.frame(exp2_fung_prev@sam_data)

#Hellinger transform and get PCA
exp2_fung_otu_hel <- decostand(exp2_fung_otu, "hellinger")
exp2_fung_otu_pca <- rda(exp2_fung_otu_hel)
summary(exp2_fung_otu_pca)$cont

#Sets priority for which labels will show if labels cover each other
priority_exp2_fung <- colSums(exp2_fung_otu)
labels_exp2_fung <- orditorp(exp2_fung_otu_pca, "sp", label = exp2_fung_tax$Genus, priority=priority_exp2_fung)

#Plot loading plot
dpi = 600
tiff("Fig 1E.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp2_fung_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (23.81% Explained)", ylab="PC2 (11.99% Explained)", ylim = c(-0.7,0.3),xlim = c(-0.8,0.5))
orditorp(exp2_fung_otu_pca, "sp", label = exp2_fung_tax$Genus, priority=priority_exp2_fung, select = (labels_exp2_fung == TRUE), cex = 0.7)
dev.off()

#Plot PCA, same method as the PCoA
factor_exp2 <- as.factor(exp2_fung_meta$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

dpi=600
tiff("Fig 1E PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.7, 0.7), font = 2, font.lab = 2, xlab="PC1 (23.81% Explained)", ylab="PC2 (11.99% Explained)", type="n")
points(exp2_fung_otu_pca, pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_fung_otu_pca, factor_exp2, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp2_fung_otu_pca, factor_exp2, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pca_df <- data.frame(exp2_fung_meta)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_hel ~ BPD, data = permanova_pca_df, method="euclidean", permutations = 10000)

#Run PERMDISP
fung_pca_dist <- vegdist(exp2_fung_otu_hel, method = "euclidean")
disp_fung_pca <- betadisper(fung_pca_dist, permanova_pca_df$BPD)
set.seed(1312)
permdisp_fung_pca <- permutest(disp_fung_pca, permutations = 10000)

print(permanova_fung_pcoa)
print(permdisp_fung_pca)


############################FIGURE 1E: 16S BETA DIVERSITY PCoA############################

exp2_bact_rel_prev <- transform_sample_counts(exp2_bact_prev, function(x) x / sum(x) )

exp2_bact_otu_rel <- as.data.frame(t(exp2_bact_rel_prev@otu_table))
exp2_bact_tax_rel <- as.data.frame(exp2_bact_rel_prev@tax_table)
exp2_bact_meta_rel <- as.data.frame(exp2_bact_rel_prev@sam_data)

#Ordinate PCoA by Bray-Curtis dissimilarity
exp2_bact_rel_bray = vegdist(exp2_bact_otu_rel, method='bray')
exp2_bact_rel_pcoa <- ape::pcoa(exp2_bact_rel_bray)
exp2_bact_rel_pcoa$values

#Set variable of interest and color palette
factor_exp2 <- as.factor(exp2_bact_meta_rel$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Plot!
dpi=600
tiff("Fig 1F.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.4), c(-0.3, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_bact_rel_pcoa$vectors[,1:2], factor_exp2, kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_bact_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df <- data.frame(exp2_bact_meta_rel)
set.seed(1312)
permanova_bact_pcoa <- vegan::adonis2(exp2_bact_otu_rel ~ BPD, data = permanova_pcoa_df, method="bray", permutations = 10000)

#PERMDISP
bact_pcoa_dist <- vegdist(exp2_bact_otu_rel, method = "bray")
disp_bact_pcoa <- betadisper(bact_pcoa_dist, permanova_pcoa_df$BPD)
set.seed(1312)
permdisp_bact_pcoa <- permutest(disp_bact_pcoa, permutations = 10000)

print(permanova_bact_pcoa)
print(permdisp_bact_pcoa)


############################FIGURE 1F: 16S BETA DIVERSITY LOADING PLOT + PCA############################

exp2_bact_otu <- as.data.frame(t(exp2_bact_prev@otu_table))
exp2_bact_tax <- as.data.frame(exp2_bact_prev@tax_table)
exp2_bact_meta <- as.data.frame(exp2_bact_prev@sam_data)

#Hellinger transform and get PCA
exp2_bact_otu_hel <- decostand(exp2_bact_otu, "hellinger")
exp2_bact_otu_pca <- rda(exp2_bact_otu_hel)
summary(exp2_bact_otu_pca)$cont

factor_exp2 <- as.factor(exp2_bact_meta$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Set label priority and get loading plot
priority_exp2_bact <- colSums(exp2_bact_otu)
labels_exp2_bact <- orditorp(exp2_bact_otu_pca, "sp", label = exp2_bact_tax$Genus, priority=priority_exp2_bact)

dpi = 600
tiff("Fig 1G.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp2_bact_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (14.03% Explained)", ylab="PC2 (9.36% Explained)", ylim = c(-0.3,0.4),xlim = c(-0.5,0.6))
orditorp(exp2_bact_otu_pca, "sp", label = exp2_bact_tax$Genus, priority=priority_exp2_bact, select = (labels_exp2_bact == TRUE), cex = 0.7)
dev.off()

#Plot PCA, same method as the PCoA
factor_exp2 <- as.factor(exp2_bact_meta$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

dpi=600
tiff("Fig 1G PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.6), c(-0.7, 0.5), font = 2, font.lab = 2, xlab="PC1 (14.03% Explained)", ylab="PC2 (9.36% Explained)", type="n")
points(exp2_bact_otu_pca, pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_bact_otu_pca, factor_exp2, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp2_bact_otu_pca, factor_exp2, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pca_df <- data.frame(exp2_bact_meta)
set.seed(1312)
permanova_bact_pca <- vegan::adonis2(exp2_bact_otu_hel ~ BPD, data = permanova_pca_df, method="euclidean", permutations = 10000)

#Run PERMDISP
bact_pca_dist <- vegdist(exp2_bact_otu_hel, method = "euclidean")
disp_bact_pca <- betadisper(bact_pca_dist, permanova_pca_df$BPD)
set.seed(1312)
permdisp_bact_pca <- permutest(disp_bact_pca, permutations = 10000)

print(permanova_bact_pca)
print(permdisp_bact_pca)



###SPIEC-EASI DATA PREP###
#Using the data aggregated to Genus level + 30% prevalence filter
exp2_combined_gen <- tax_glom(exp2_combined, taxrank = "Genus")
exp2_combined_30_gen <- filter_taxa(exp2_combined_gen, function(x) sum(x >= 1) > (0.30*length(x)), TRUE)
exp2_combined_30_gen_noNA = subset_taxa(exp2_combined_30_gen, Genus!="NA")
exp2_combined_30_gen_noNA_bpd <- subset_samples(exp2_combined_30_gen_noNA, BPD=="BPD")
exp2_combined_30_gen_noNA_pprd <- subset_samples(exp2_combined_30_gen_noNA, BPD=="No_BPD")

exp2_combined_16s = subset_taxa(exp2_combined_30_gen, Kingdom=="Bacteria")
exp2_combined_16s_bpd = subset_samples(exp2_combined_16s, BPD=="BPD")
exp2_combined_16s_pprd = subset_samples(exp2_combined_16s, BPD=="No_BPD")
exp2_combined_its = subset_taxa(exp2_combined_30_gen, Kingdom=="Fungi")
exp2_combined_its_bpd = subset_samples(exp2_combined_its, BPD=="BPD")
exp2_combined_its_pprd = subset_samples(exp2_combined_its, BPD=="No_BPD")


############################FIGURE 1G: NOBPD SPIEC-EASI NETWORK############################

###NOBPD SPIEC-EASI###
set.seed(1312)
se.exp2.gen.pprd <- spiec.easi(list(exp2_combined_16s_pprd, exp2_combined_its_pprd), method='mb', nlambda=99,
                               lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

#NoBPD only network plot with labels
spiec.graph.gen.pprd=adj2igraph(getRefit(se.exp2.gen.pprd), vertex.attr=list(name=taxa_names(exp2_combined_30_gen_noNA)))
prev_gen_spiec_plot_pprd <- plot_network(spiec.graph.gen.pprd, exp2_combined_30_gen_noNA, type='taxa',label=NULL,color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 2)
#Save

#NoBPD only network plot without labels
prev_gen_spiec_plot_pprd + scale_color_manual(values = spiec.colors)
#Save

#Get the nodes and edge counts
nodes_pprd <- gorder(spiec.graph.gen.pprd)
edges_pprd <- gsize(spiec.graph.gen.pprd)


############################FIGURE 1H: BPD SPIEC-EASI NETWORK############################

###BPD SPIEC-EASI###
set.seed(1312)
se.exp2.gen.bpd <- spiec.easi(list(exp2_combined_16s_bpd, exp2_combined_its_bpd), method='mb', nlambda=99,
                              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

#Color palette, bacteria in gray and fungi in green
spiec.colors <- c("#aaaaaa","#74c476")

###PLOT NETWORK###

#BPD only network plot with labels
spiec.graph.gen.bpd=adj2igraph(getRefit(se.exp2.gen.bpd), vertex.attr=list(name=taxa_names(exp2_combined_30_gen_noNA)))
prev_gen_spiec_plot_bpd <- plot_network(spiec.graph.gen.bpd, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 2)
#Save

#BPD only network plot with no labels
prev_gen_spiec_plot_bpd <- plot_network(spiec.graph.gen.bpd, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd + scale_color_manual(values = spiec.colors)
#Save

#Get number of nodes and edges
nodes_bpd <- gorder(spiec.graph.gen.bpd)
edges_bpd <- gsize(spiec.graph.gen.bpd)
