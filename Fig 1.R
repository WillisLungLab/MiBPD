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

############################FIGURE 1A: ITS TOP 20 BARPLOT############################

#Merge all samples together by BPD status and aggregate to Genus level
exp2_fung_merged = merge_samples(exp2_fung_prev, "BPD")
exp2_fung_gen_merged <- tax_glom(exp2_fung_merged, taxrank = 'Genus')

#Identify top 19 genera and rename everything else to "Other"
top20_fung_merged_list <- names(sort(taxa_sums(exp2_fung_gen_merged), decreasing=TRUE)[1:19])
top20_merged_fung_rel <- transform_sample_counts(exp2_fung_gen_merged, function(x) x / sum(x) )
top20_merged_fung_df <- psmelt(top20_merged_fung_rel)
top20_merged_fung_df[!(top20_merged_fung_df$OTU %in% top20_fung_merged_list),]$Genus <- 'Other'


###Barplot of top 20 fungal genera###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(20)

barplot_gen_bpd_its <- ggplot(top20_merged_fung_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_its + scale_fill_manual(values = barplot_colors)
#Save

############################FIGURE 1B: ITS ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
#Get Chao1 and Shannon diversity
richness_est_fung <- estimate_richness(fr_fung, measures = c("Chao1", "Shannon"))

#Wilcoxon rank-sum test for significance
wilcox_alpha_fung <- t(sapply(richness_est_fung, function(x) unlist(wilcox.test(x~sample_data(fr_fung)$BPD)[c("estimate","p.value","statistic","conf.int")])))
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


############################FIGURE 1C: MULTIKINGDOM ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_fung <- prune_taxa(taxa_sums(exp2_combined_rough) > 0, exp2_combined_rough)
#Get Chao1 and Shannon diversity
richness_est_combined <- estimate_richness(fr_combined, measures = c("Chao1", "Shannon"))

#Wilcoxon rank-sum test for significance
wilcox_alpha_combined <- t(sapply(richness_est_combined, function(x) unlist(wilcox.test(x~sample_data(fr_combined)$BPD)[c("estimate","p.value","statistic","conf.int")])))
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

############################FIGURE 1D: ITS BETA DIVERSITY PCoA############################

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
tiff("Beta Diversity PCoA BPD.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1(9.44% Explained)", ylab="PC2(7.31% Explained)", type="n")
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

############################FIGURE 1E: ITS BETA DIVERSITY LOADING PLOT + PCA############################

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
tiff("MiBPD ITS Loading Plot ASV.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp2_fung_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (22.73% Explained)", ylab="PC2 (14.66% Explained)", ylim = c(-0.7,0.3),xlim = c(-0.6,0.8))
orditorp(exp2_fung_otu_pca, "sp", label = exp2_fung_tax$Genus, priority=priority_exp2_fung, select = (labels_exp2_fung == TRUE), cex = 0.7)
dev.off()

#Plot PCA, same method as the PCoA
factor_exp2 <- as.factor(exp2_fung_meta$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

dpi=600
tiff("MiBPD ITS PCA BPD ASV.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.7, 0.5), c(-0.7, 0.7), font = 2, font.lab = 2, xlab="PC1 (22.73% Explained)", ylab="PC2 (14.66% Explained)", type="n")
points(exp2_fung_otu_pca, pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_fung_otu_pca, factor_exp2, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp2_fung_otu_pca, factor_exp2, label = TRUE)
dev.off()

############################FIGURE 1F: MULTIKINGDOM BETA DIVERSITY PCoA############################

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)

#Ordinate PCoA by Bray-Curtis dissimilarity
exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
exp2_combined_rel_pcoa$values

#Set variable of interest and color palette
factor_exp2 <- as.factor(exp2_combined_meta_rel$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Plot!
dpi=600
tiff("Beta Diversity PCoA BPD ASV New CI.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.4), c(-0.4, 0.5), font = 2, font.lab = 2, xlab="PC1(5.41% Explained)", ylab="PC2(4.35% Explained)", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2, kind = "se", conf = 0.95, col = pca_colors)
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

############################FIGURE 1G: MULTIKINGDOM BETA DIVERSITY LOADING PLOT + PCA############################

exp2_combined_otu <- as.data.frame(t(exp2_combined_prev@otu_table))
exp2_combined_tax <- as.data.frame(exp2_combined_prev@tax_table)
exp2_combined_meta <- as.data.frame(exp2_combined_prev@sam_data)

#Hellinger transform and get PCA
exp2_combined_otu_hel <- decostand(exp2_combined_otu, "hellinger")
exp2_combined_otu_pca <- rda(exp2_combined_otu_hel)
summary(exp2_combined_otu_pca)$cont

factor_exp2 <- as.factor(exp2_combined_meta$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Set label priority and get loading plot
priority_exp2_combined <- colSums(exp2_combined_otu)
labels_exp2_combined <- orditorp(exp2_combined_otu_pca, "sp", label = exp2_combined_tax$Genus, priority=priority_exp2_combined)

dpi = 600
tiff("MiBPD Combined Loading Plot ASV.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp2_combined_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (13.09% Explained)", ylab="PC2 (9.31% Explained)", ylim = c(-0.3,0.6),xlim = c(-0.6,0.6))
orditorp(exp2_combined_otu_pca, "sp", label = exp2_combined_tax$Genus, priority=priority_exp2_combined, select = (labels_exp2_combined == TRUE), cex = 0.7)
dev.off()


############################FIGURE 1H: BPD SPIEC-EASI NETWORK############################

###DATA PREP###
#Using the data aggregated to Genus level + 30% prevalence filter 
exp2_combined_30_gen <- filter_taxa(exp2_combined_gen, function(x) sum(x >= 1) > (0.30*length(x)), TRUE)
exp2_combined_gen_prev_noNA = subset_taxa(exp2_combined_30_gen, Genus!="NA")

#Pull individual components out of the phyloseq object
exp2_prev_gen_otu <- as.data.frame(exp2_combined_gen_prev_noNA@otu_table)
exp2_prev_gen_tax <- as.data.frame(exp2_combined_gen_prev_noNA@tax_table)
exp2_prev_gen_meta <- as.data.frame(exp2_combined_gen_prev_noNA@sam_data)

#Subset to just bacteria
exp2_prev_gen_otu_bact <- exp2_prev_gen_otu[c(1:60),]
exp2_prev_gen_tax_bact <- exp2_prev_gen_tax[c(1:60),]

#Subset to just fungi
exp2_prev_gen_otu_fung <- exp2_prev_gen_otu[c(61:80),]
exp2_prev_gen_tax_fung <- exp2_prev_gen_tax[c(61:80),]

###BPD BACTERIA DATA PREP###
exp2_prev_gen_tax_bact <- as.matrix(exp2_prev_gen_tax_bact)
exp2_prev_gen_bact_bpd <- exp2_prev_gen_otu_bact[,which(exp2_prev_gen_meta$BPD == "BPD")]

#Create phyloseq object for BPD 16s
OTU_prev_gen_bact_bpd = otu_table(exp2_prev_gen_bact_bpd, taxa_are_rows = TRUE)
TAX_prev_gen_bact_bpd = tax_table(exp2_prev_gen_tax_bact)
samples_prev_gen_bact_bpd = sample_data(exp2_prev_gen_meta)

exp2_prev_bact_bpd <- phyloseq(OTU_prev_gen_bact_bpd, TAX_prev_gen_bact_bpd, samples_prev_gen_bact_bpd)


###BPD FUNGI DATA PREP###
exp2_prev_gen_tax_fung <- as.matrix(exp2_prev_gen_tax_fung)
exp2_prev_gen_fung_bpd <- exp2_prev_gen_otu_fung[,which(exp2_prev_gen_meta$BPD == "BPD")]

#BPD ITS phyloseq object
OTU_prev_gen_fung_bpd = otu_table(exp2_prev_gen_fung_bpd, taxa_are_rows = TRUE)
TAX_prev_gen_fung_bpd = tax_table(exp2_prev_gen_tax_fung)
samples_prev_gen_fung_bpd = sample_data(exp2_prev_gen_meta)

exp2_prev_fung_bpd <- phyloseq(OTU_prev_gen_fung_bpd, TAX_prev_gen_fung_bpd, samples_prev_gen_fung_bpd)

###BPD SPIEC-EASI###
set.seed(1312)
se.exp2.gen.bpd <- spiec.easi(list(exp2_prev_bact_bpd, exp2_prev_fung_bpd), method='mb', nlambda=99,
                              lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

#Color palette, bacteria in gray and fungi in green
spiec.colors <- c("#aaaaaa","#74c476")

###PLOT NETWORK###
library(igraph)
##BPD ONLY##
#BPD only network plot with no taxa labels
spiec.graph.gen.bpd=adj2igraph(getRefit(se.exp2.gen.bpd), vertex.attr=list(name=taxa_names(exp2_combined_gen_prev_noNA)))
prev_gen_spiec_plot_bpd <- plot_network(spiec.graph.gen.bpd, exp2_combined_gen_prev_noNA, type='taxa',label="Genus", color="Kingdom") + theme(legend.position = "none")
prev_gen_spiec_plot_bpd + scale_color_manual(values = spiec.colors)
#Save

#BPD only network plot with taxa labels
prev_gen_spiec_plot_bpd_legend <- plot_network(spiec.graph.gen.bpd, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

#Get number of notes, number of edges, betweenness values, and degree values
nodes_bpd <- gorder(spiec.graph.gen.bpd)
edges_bpd <- gsize(spiec.graph.gen.bpd)
betweenness_bpd <- as.list(betweenness(spiec.graph.gen.bpd, normalized = TRUE))
degree_bpd <- as.data.frame(degree(spiec.graph.gen.bpd))

############################FIGURE 1I: PPRD SPIEC-EASI NETWORK############################

###PPRD BACTERIA DATA PREP###
exp2_prev_gen_bact_pprd <- exp2_prev_gen_otu_bact[,which(exp2_prev_gen_meta$BPD == "PPRD")]

#New phyloseq object for PPRD 16s
OTU_prev_gen_bact_pprd = otu_table(exp2_prev_gen_bact_pprd, taxa_are_rows = TRUE)
TAX_prev_gen_bact_pprd = tax_table(exp2_prev_gen_tax_bact)
samples_prev_gen_bact_pprd = sample_data(exp2_prev_gen_meta)

exp2_prev_bact_pprd <- phyloseq(OTU_prev_gen_bact_pprd, TAX_prev_gen_bact_pprd, samples_prev_gen_bact_pprd)

###PPRD FUNGI DATA PREP###
exp2_prev_gen_fung_pprd <- exp2_prev_gen_otu_fung[,which(exp2_prev_gen_meta$BPD == "PPRD")]

#PPRD ITS phyloseq object
OTU_prev_gen_fung_pprd = otu_table(exp2_prev_gen_fung_pprd, taxa_are_rows = TRUE)
TAX_prev_gen_fung_pprd = tax_table(exp2_prev_gen_tax_fung)
samples_prev_gen_fung_pprd = sample_data(exp2_prev_gen_meta)

exp2_prev_fung_pprd <- phyloseq(OTU_prev_gen_fung_pprd, TAX_prev_gen_fung_pprd, samples_prev_gen_fung_pprd)

###PPRD SPIEC-EASI###
set.seed(1312)
se.exp2.gen.pprd <- spiec.easi(list(exp2_prev_bact_pprd, exp2_prev_fung_pprd), method='mb', nlambda=99,
                               lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

#PPRD only network plot with no taxa labels
prev_gen_spiec_plot_pprd <- plot_network(spiec.graph.gen.pprd, exp2_combined_gen_prev_noNA, type='taxa',label=NULL,color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd + scale_color_manual(values = spiec.colors)
#Save
                                    
#PPRD only network plot with taxa labels
prev_gen_spiec_plot_pprd <- plot_network(spiec.graph.gen.pprd, exp2_combined_gen_prev_noNA, type='taxa',label=NULL,color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

#Get the nodes, edges, betweenness, and degree for the network
nodes_pprd <- gorder(spiec.graph.gen.pprd)
edges_pprd <- gsize(spiec.graph.gen.pprd)
betweenness_pprd <- as.list(betweenness(spiec.graph.gen.pprd, normalized = TRUE))
degree_pprd <- as.data.frame(degree(spiec.graph.gen.pprd))

