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
library(FSA)
library(ecole)

############################FIGURE 2A: ITS TOP 20 BARPLOT############################

#Merge all samples together by BPD status and aggregate to Genus level
exp2_fung_merged = merge_samples(exp2_fung_prev, "BPD_Sex")
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

############################FIGURE 2B: ITS ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
#Get Chao1 and Shannon diversity
richness_est_fung<- estimate_richness(fr_fung, measures = c("Chao1", "Shannon"))

#Kruskal-Wallis test for significance and Dunn test as a post-hoc test
kruskal_alpha_fung <- t(sapply(richness_est_fung, function(x) unlist(kruskal.test(x~sample_data(fr_fung)$BPD_Sex)[c("estimate","p.value","statistic","conf.int")])))
dunnTest(richness_est_fung$Chao1 ~ sample_data(fr_fung)$BPD_Sex,
         data=richness_est_fung,
         method="bh")
dunnTest(richness_est_fung$Shannon ~ sample_data(fr_fung)$BPD_Sex,
         data=richness_est_fung,
         method="bh")

richness_est_fung <- richness_est_fung %>%
  mutate(
    BPD = md.fung$BPD,
    Sex = md.fung$Sex,
    BPD_Sex = md.fung$BPD_Sex,
    Oxygen = md.fung$Oxygen_Days
  )
#Save as .csv and use in GraphPad


############################FIGURE 2C: MULTIKINGDOM ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_combined <- prune_taxa(taxa_sums(exp2_combined_rough) > 0, exp2_combined_rough)
#Get Chao1 and Shannon diversity
richness_est_combined<- estimate_richness(fr_combined, measures = c("Chao1", "Shannon"))

#Kruskal-Wallis test for significance and Dunn test as a post-hoc test
kruskal_alpha_combined <- t(sapply(richness_est_combined, function(x) unlist(kruskal.test(x~sample_data(fr_combined)$BPD_Sex)[c("estimate","p.value","statistic","conf.int")])))
dunnTest(richness_est_combined$Chao1 ~ sample_data(fr_combined)$BPD_Sex,
         data=richness_est_combined,
         method="bh")
dunnTest(richness_est_combined$Shannon ~ sample_data(fr_combined)$BPD_Sex,
         data=richness_est_combined,
         method="bh")

richness_est_combined <- richness_est_combined %>%
  mutate(
    BPD = md.combined$BPD,
    Sex = md.combined$Sex,
    BPD_Sex = md.combined$BPD_Sex,
    Oxygen = md.combined$Oxygen_Days
  )
#Save as .csv and use in GraphPad

############################FIGURE 2D: ITS BETA DIVERSITY PCoA############################

#Convert to relative abundance instead of counts
exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)

#Get Bray-Curtis dissimilarity and PCoA ordination
exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
exp2_fung_rel_pcoa$values

#Set factor of interest
factor_exp2 <- as.factor(exp2_fung_meta_rel$BPD_Sex)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

#Spider plot with ellipses for 95% confidence interval
dpi=600
tiff("Beta Diversity PCoA BPD_SAAB.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1(9.44% Explained)", ylab="PC2(7.31% Explained)", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2, kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ BPD_Sex, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run pairwise PERMANOVA and FDR adjust p values
pairwise_permanova_fung_pcoa <- permanova_pairwise(exp2_fung_otu_rel, grp = permanova_pcoa_df$BPD_Sex, permutations = 10000, method = "bray", padj = "fdr")

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$BPD_Sex)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(pairwise_permanova_fung_pcoa)
print(permdisp_fung_pcoa)

############################FIGURE 2E: MULTIKINGDOM BETA DIVERSITY PCoA############################

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)

#Bray-Curtis dissimilarity and PCoA ordination
exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
exp2_combined_rel_pcoa$values

#Set colors and plot spider plot
factor_exp2 <- as.factor(exp2_combined_meta_rel$BPD_Sex)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Beta Diversity PCoA BPD_SAAB.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.5), c(-0.4, 0.5), font = 2, font.lab = 2, xlab="PC1(5.41% Explained)", ylab="PC2(4.35% Explained)", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], factor2_exp2, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df <- data.frame(exp2_combined_meta_rel)
set.seed(1312)
permanova_combined_pcoa <- vegan::adonis2(exp2_combined_otu_rel ~ BPD_Sex, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Pairwise PERMANOVA
pairwise_permanova_combined_pcoa <- permanova_pairwise(exp2_combined_otu_rel, grp = permanova_pcoa_df$BPD_Sex, permutations = 10000, method = "bray", padj = "fdr")

#PERMDISP
combined_pcoa_dist <- vegdist(exp2_combined_otu_rel, method = "bray")
disp_combined_pcoa <- betadisper(combined_pcoa_dist, permanova_pcoa_df$BPD_Sex)
set.seed(1312)
permdisp_combined_pcoa <- permutest(disp_combined_pcoa, permutations = 10000)

print(permanova_combined_pcoa)
print(pairwise_permanova_combined_pcoa)
print(permdisp_combined_pcoa)

############################FIGURE 2F: PPRD AMAB SPIEC-EASI NETWORK############################

###DATA PREP###
exp2_combined_30_gen <- filter_taxa(exp2_combined_gen, function(x) sum(x >= 1) > (0.30*length(x)), TRUE)
exp2_combined_gen_prev_noNA = subset_taxa(exp2_combined_30_gen, Genus!="NA")

exp2_prev_gen_otu <- as.data.frame(exp2_combined_gen_prev_noNA@otu_table)
exp2_prev_gen_tax <- as.data.frame(exp2_combined_gen_prev_noNA@tax_table)
exp2_prev_gen_meta <- as.data.frame(exp2_combined_gen_prev_noNA@sam_data)

exp2_prev_gen_otu_bact <- exp2_prev_gen_otu[c(1:60),]
exp2_prev_gen_tax_bact <- exp2_prev_gen_tax[c(1:60),]

exp2_prev_gen_otu_fung <- exp2_prev_gen_otu[c(61:80),]
exp2_prev_gen_tax_fung <- exp2_prev_gen_tax[c(61:80),]

exp2_prev_gen_tax_bact <- as.matrix(exp2_prev_gen_tax_bact)
exp2_prev_gen_tax_fung <- as.matrix(exp2_prev_gen_tax_fung)

exp2_prev_gen_bact_pprd_amab <- exp2_prev_gen_otu_bact[,which(exp2_prev_gen_meta$BPD_Sex == "PPRD_AMAB")]
exp2_prev_gen_fung_pprd_amab <- exp2_prev_gen_otu_fung[,which(exp2_prev_gen_meta$BPD_Sex == "PPRD_AMAB")]

###PPRD AMAB Bacteria###
OTU_prev_gen_bact_pprd_amab = otu_table(exp2_prev_gen_bact_pprd_amab, taxa_are_rows = TRUE)
TAX_prev_gen_bact_pprd_amab = tax_table(exp2_prev_gen_tax_bact)
samples_prev_gen_bact_pprd_amab = sample_data(exp2_prev_gen_meta)

exp2_prev_bact_pprd_amab <- phyloseq(OTU_prev_gen_bact_pprd_amab, TAX_prev_gen_bact_pprd_amab, samples_prev_gen_bact_pprd_amab)

###PPRD AMAB Fungi###
OTU_prev_gen_fung_pprd_amab = otu_table(exp2_prev_gen_fung_pprd_amab, taxa_are_rows = TRUE)
TAX_prev_gen_fung_pprd_amab = tax_table(exp2_prev_gen_tax_fung)
samples_prev_gen_fung_pprd_amab = sample_data(exp2_prev_gen_meta)

exp2_prev_fung_pprd_amab <- phyloseq(OTU_prev_gen_fung_pprd_amab, TAX_prev_gen_fung_pprd_amab, samples_prev_gen_fung_pprd_amab)

###Run SPIEC-EASI###
set.seed(1312)
se.exp2.gen.pprd.amab <- spiec.easi(list(exp2_prev_bact_pprd_amab, exp2_prev_fung_pprd_amab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))
spiec.graph.gen.pprd.amab=adj2igraph(getRefit(se.exp2.gen.pprd.amab), vertex.attr=list(name=taxa_names(exp2_combined_gen_prev_noNA)))

#Plot network with and without labels
prev_gen_spiec_plot_pprd_amab_nolegend <- plot_network(spiec.graph.gen.pprd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom")
prev_gen_spiec_plot_pprd_amab_nolegend + scale_color_manual(values = spiec.colors)
#Save

prev_gen_spiec_plot_pprd_amab_legend <- plot_network(spiec.graph.gen.pprd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

#Get relevant stats
nodes_pprd_amab <- gorder(spiec.graph.gen.pprd.amab)
edges_pprd_amab <- gsize(spiec.graph.gen.pprd.amab)
betweenness_pprd_amab <- as.list(betweenness(spiec.graph.gen.pprd.amab, normalized = TRUE))
degree_pprd_amab <- as.data.frame(degree(spiec.graph.gen.pprd.amab))

############################FIGURE 2G: PPRD AFAB SPIEC-EASI NETWORK############################

exp2_prev_gen_bact_pprd_afab <- exp2_prev_gen_otu_bact[,which(exp2_prev_gen_meta$BPD_Sex == "PPRD_AFAB")]
exp2_prev_gen_fung_pprd_afab <- exp2_prev_gen_otu_fung[,which(exp2_prev_gen_meta$BPD_Sex == "PPRD_AFAB")]

###PPRD AFAB Bacteria###
OTU_prev_gen_bact_pprd_afab = otu_table(exp2_prev_gen_bact_pprd_afab, taxa_are_rows = TRUE)
TAX_prev_gen_bact_pprd_afab = tax_table(exp2_prev_gen_tax_bact)
samples_prev_gen_bact_pprd_afab = sample_data(exp2_prev_gen_meta)

exp2_prev_bact_pprd_afab <- phyloseq(OTU_prev_gen_bact_pprd_afab, TAX_prev_gen_bact_pprd_afab, samples_prev_gen_bact_pprd_afab)

###PPRD AFAB Fungi###
OTU_prev_gen_fung_pprd_afab = otu_table(exp2_prev_gen_fung_pprd_afab, taxa_are_rows = TRUE)
TAX_prev_gen_fung_pprd_afab = tax_table(exp2_prev_gen_tax_fung)
samples_prev_gen_fung_pprd_afab = sample_data(exp2_prev_gen_meta)

exp2_prev_fung_pprd_afab <- phyloseq(OTU_prev_gen_fung_pprd_afab, TAX_prev_gen_fung_pprd_afab, samples_prev_gen_fung_pprd_afab)

#Run SPIEC-EASI
se.exp2.gen.pprd.afab <- spiec.easi(list(exp2_prev_bact_pprd_afab, exp2_prev_fung_pprd_afab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

spiec.graph.gen.pprd.afab=adj2igraph(getRefit(se.exp2.gen.pprd.afab), vertex.attr=list(name=taxa_names(exp2_combined_gen_prev_noNA)))

#Plot network with and without labels
prev_gen_spiec_plot_pprd_afab_nolegend <- plot_network(spiec.graph.gen.pprd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom")
prev_gen_spiec_plot_pprd_afab_nolegend + scale_color_manual(values = spiec.colors)
#Save

prev_gen_spiec_plot_pprd_afab_legend <- plot_network(spiec.graph.gen.pprd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

nodes_pprd_afab <- gorder(spiec.graph.gen.pprd.afab)
edges_pprd_afab <- gsize(spiec.graph.gen.pprd.afab)
betweenness_pprd_afab <- as.list(betweenness(spiec.graph.gen.pprd.afab, normalized = TRUE))
degree_pprd_afab <- as.data.frame(degree(spiec.graph.gen.pprd.afab))


############################FIGURE 2H: BPD AMAB SPIEC-EASI NETWORK############################

exp2_prev_gen_bact_bpd_amab <- exp2_prev_gen_otu_bact[,which(exp2_prev_gen_meta$BPD_Sex == "BPD_AMAB")]
exp2_prev_gen_fung_bpd_amab <- exp2_prev_gen_otu_fung[,which(exp2_prev_gen_meta$BPD_Sex == "BPD_AMAB")]

###BPD AMAB Bacteria###
OTU_prev_gen_bact_bpd_amab = otu_table(exp2_prev_gen_bact_bpd_amab, taxa_are_rows = TRUE)
TAX_prev_gen_bact_bpd_amab = tax_table(exp2_prev_gen_tax_bact)
samples_prev_gen_bact_bpd_amab = sample_data(exp2_prev_gen_meta)

exp2_prev_bact_bpd_amab <- phyloseq(OTU_prev_gen_bact_bpd_amab, TAX_prev_gen_bact_bpd_amab, samples_prev_gen_bact_bpd_amab)

###BPD AMAB Fungi###
OTU_prev_gen_fung_bpd_amab = otu_table(exp2_prev_gen_fung_bpd_amab, taxa_are_rows = TRUE)
TAX_prev_gen_fung_bpd_amab = tax_table(exp2_prev_gen_tax_fung)
samples_prev_gen_fung_bpd_amab = sample_data(exp2_prev_gen_meta)

exp2_prev_fung_bpd_amab <- phyloseq(OTU_prev_gen_fung_bpd_amab, TAX_prev_gen_fung_bpd_amab, samples_prev_gen_fung_bpd_amab)

set.seed(1312)
se.exp2.gen.bpd.amab <- spiec.easi(list(exp2_prev_bact_bpd_amab, exp2_prev_fung_bpd_amab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))
spiec.graph.gen.bpd.amab=adj2igraph(getRefit(se.exp2.gen.bpd.amab), vertex.attr=list(name=taxa_names(exp2_combined_gen_prev_noNA)))

#Plot network with and without labels
prev_gen_spiec_plot_bpd_amab_nolegend <- plot_network(spiec.graph.gen.bpd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom")
prev_gen_spiec_plot_bpd_amab_nolegend + scale_color_manual(values = spiec.colors)
#Save

prev_gen_spiec_plot_bpd_amab_legend <- plot_network(spiec.graph.gen.bpd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save


nodes_bpd_amab <- gorder(spiec.graph.gen.bpd.amab)
edges_bpd_amab <- gsize(spiec.graph.gen.bpd.amab)
betweenness_bpd_amab <- as.list(betweenness(spiec.graph.gen.bpd.amab, normalized = TRUE))
degree_bpd_amab <- as.data.frame(degree(spiec.graph.gen.bpd.amab))

############################FIGURE 2I: BPD AFAB SPIEC-EASI NETWORK############################

exp2_prev_gen_bact_bpd_afab <- exp2_prev_gen_otu_bact[,which(exp2_prev_gen_meta$BPD_Sex == "bpd_AFAB")]
exp2_prev_gen_fung_bpd_afab <- exp2_prev_gen_otu_fung[,which(exp2_prev_gen_meta$BPD_Sex == "bpd_AFAB")]

###Bacteria###
OTU_prev_gen_bact_bpd_afab = otu_table(exp2_prev_gen_bact_bpd_afab, taxa_are_rows = TRUE)
TAX_prev_gen_bact_bpd_afab = tax_table(exp2_prev_gen_tax_bact)
samples_prev_gen_bact_bpd_afab = sample_data(exp2_prev_gen_meta)

exp2_prev_bact_bpd_afab <- phyloseq(OTU_prev_gen_bact_bpd_afab, TAX_prev_gen_bact_bpd_afab, samples_prev_gen_bact_bpd_afab)

###Fungi###
OTU_prev_gen_fung_bpd_afab = otu_table(exp2_prev_gen_fung_bpd_afab, taxa_are_rows = TRUE)
TAX_prev_gen_fung_bpd_afab = tax_table(exp2_prev_gen_tax_fung)
samples_prev_gen_fung_bpd_afab = sample_data(exp2_prev_gen_meta)

exp2_prev_fung_bpd_afab <- phyloseq(OTU_prev_gen_fung_bpd_afab, TAX_prev_gen_fung_bpd_afab, samples_prev_gen_fung_bpd_afab)

se.exp2.gen.bpd.afab <- spiec.easi(list(exp2_prev_bact_bpd_afab, exp2_prev_fung_bpd_afab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))


#Plot network with and without labels
prev_gen_spiec_plot_bpd_afab_nolegend <- plot_network(spiec.graph.gen.bpd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom")
prev_gen_spiec_plot_bpd_afab_nolegend + scale_color_manual(values = spiec.colors)
#Save

prev_gen_spiec_plot_bpd_afab_legend <- plot_network(spiec.graph.gen.bpd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

nodes_bpd_afab <- gorder(spiec.graph.gen.bpd.afab)
edges_bpd_afab <- gsize(spiec.graph.gen.bpd.afab)
betweenness_bpd_afab <- as.list(betweenness(spiec.graph.gen.bpd.afab, normalized = TRUE))
degree_bpd_afab <- as.data.frame(degree(spiec.graph.gen.bpd.afab))