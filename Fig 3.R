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

############################FIGURE 3A: ITS TOP 20 BARPLOT############################

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

############################FIGURE 3B: ITS ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_fung <- prune_taxa(taxa_sums(exp2_fung_rough) > 0, exp2_fung_rough)
#Get Simpson and Shannon diversity
richness_est_fung<- estimate_richness(fr_fung, measures = c("Simpson", "Shannon"))

#Kruskal-Wallis test for significance and Dunn test as a post-hoc test
set.seed(1312)
kruskal_alpha_fung <- t(sapply(richness_est_fung, function(x) unlist(kruskal.test(x~sample_data(fr_fung)$BPD_Sex)[c("estimate","p.value","statistic","conf.int")])))
kruskal_alpha_fung
set.seed(1312)
dunnTest(richness_est_fung$Simpson ~ sample_data(fr_fung)$BPD_Sex,
         data=richness_est_fung,
         method="bh")
set.seed(1312)
dunnTest(richness_est_fung$Shannon ~ sample_data(fr_fung)$BPD_Sex,
         data=richness_est_fung,
         method="bh")

richness_est_fung <- richness_est_fung %>%
  mutate(
    BPD = sample_data(fr_fung)$BPD,
    Sex = sample_data(fr_fung)$Sex,
    BPD_Sex = sample_data(fr_fung)$BPD_Sex,
    BPD_Severity = sample_data(fr_fung)$BPD_Severity,
    Oxygen = sample_data(fr_fung)$Oxygen
  )
#Save as .csv and use in GraphPad


############################FIGURE 3C: MULTIKINGDOM ALPHA DIVERSITY############################

#Remove taxa that aren't present in any sample
fr_bact <- prune_taxa(taxa_sums(exp2_bact_rough) > 0, exp2_bact_rough)
#Get Simpson and Shannon diversity
richness_est_bact<- estimate_richness(fr_bact, measures = c("Simpson", "Shannon"))

#Kruskal-Wallis test for significance and Dunn test as a post-hoc test
kruskal_alpha_bact <- t(sapply(richness_est_bact, function(x) unlist(kruskal.test(x~sample_data(fr_bact)$BPD_Sex)[c("estimate","p.value","statistic","conf.int")])))
kruskal_alpha_bact
dunnTest(richness_est_bact$Simpson ~ sample_data(fr_bact)$BPD_Sex,
         data=richness_est_bact,
         method="bh")
dunnTest(richness_est_bact$Shannon ~ sample_data(fr_bact)$BPD_Sex,
         data=richness_est_bact,
         method="bh")

richness_est_bact <- richness_est_bact %>%
  mutate(
    BPD = sample_data(fr_bact)$BPD,
    Sex = sample_data(fr_bact)$Sex,
    BPD_Sex = sample_data(fr_bact)$BPD_Sex,
    BPD_Severity = sample_data(fr_bact)$BPD_Severity,
    Oxygen = sample_data(fr_bact)$Oxygen
  )
#Save as .csv


############################FIGURE 3D: ITS BETA DIVERSITY PCoA############################

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
tiff("Fig 3D.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_fung_rel_pcoa$vectors[,1:2], type_exp2, kind = "se", conf = 0.95, col = pca_colors)
ordispider(exp2_fung_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#Run PERMANOVA
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel)
set.seed(1312)
permanova_fung_pcoa <- vegan::adonis2(exp2_fung_otu_rel ~ BPD_Sex, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Run pairwise PERMANOVA and FDR adjust p values
set.seed(1312)
pairwise_permanova_fung_pcoa <- permanova_pairwise(exp2_fung_otu_rel, grp = permanova_pcoa_df$BPD_Sex, permutations = 10000, method = "bray", padj = "fdr")

#Run PERMDISP
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df$BPD_Sex)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_pcoa)
print(pairwise_permanova_fung_pcoa)
print(permdisp_fung_pcoa)

############################FIGURE 3E: MULTIKINGDOM BETA DIVERSITY PCoA############################

exp2_bact_rel_prev <- transform_sample_counts(exp2_bact_prev, function(x) x / sum(x) )

exp2_bact_otu_rel <- as.data.frame(t(exp2_bact_rel_prev@otu_table))
exp2_bact_tax_rel <- as.data.frame(exp2_bact_rel_prev@tax_table)
exp2_bact_meta_rel <- as.data.frame(exp2_bact_rel_prev@sam_data)

#Bray-Curtis dissimilarity and PCoA ordination
exp2_bact_rel_bray = vegdist(exp2_bact_otu_rel, method='bray')
exp2_bact_rel_pcoa <- ape::pcoa(exp2_bact_rel_bray)
exp2_bact_rel_pcoa$values

#Set colors and plot spider plot
factor_exp2 <- as.factor(exp2_bact_meta_rel$BPD_Sex)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#007016","#94D57F","#94D57F")

dpi=600
tiff("Fig 3E.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.4), c(-0.3, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_bact_rel_pcoa$vectors[,1:2], factor_exp2, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp2_bact_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()

#PERMANOVA
permanova_pcoa_df <- data.frame(exp2_bact_meta_rel)
set.seed(1312)
permanova_bact_pcoa <- vegan::adonis2(exp2_bact_otu_rel ~ BPD_Sex, data = permanova_pcoa_df, method="bray", permutations = 10000)

#Pairwise PERMANOVA
set.seed(1312)
pairwise_permanova_bact_pcoa <- permanova_pairwise(exp2_bact_otu_rel, grp = permanova_pcoa_df$BPD_Sex, permutations = 10000, method = "bray", padj = "fdr")

#PERMDISP
bact_pcoa_dist <- vegdist(exp2_bact_otu_rel, method = "bray")
disp_bact_pcoa <- betadisper(bact_pcoa_dist, permanova_pcoa_df$BPD_Sex)
set.seed(1312)
permdisp_bact_pcoa <- permutest(disp_bact_pcoa, permutations = 10000)

print(permanova_bact_pcoa)
print(pairwise_permanova_bact_pcoa)
print(permdisp_bact_pcoa)

############################FIGURE 3F: PPRD AMAB SPIEC-EASI NETWORK############################

###DATA PREP###
exp2_combined_30_gen <- filter_taxa(exp2_combined_30_gen, function(x) sum(x >= 1) > (0.30*length(x)), TRUE)
exp2_combined_30_gen_noNA = subset_taxa(exp2_combined_30_gen, Genus!="NA")

exp2_combined_16s = subset_taxa(exp2_combined_30_gen, Kingdom=="Bacteria")

exp2_combined_16s_bpd = subset_samples(exp2_combined_16s, BPD=="BPD")
exp2_combined_16s_bpd_amab = subset_samples(exp2_combined_16s_bpd, BPD_Sex=="BPD_AMAB")
exp2_combined_16s_bpd_afab = subset_samples(exp2_combined_16s_bpd, BPD_Sex=="BPD_AFAB")

exp2_combined_16s_pprd = subset_samples(exp2_combined_16s, BPD=="No_BPD")
exp2_combined_16s_pprd_amab = subset_samples(exp2_combined_16s_pprd, BPD_Sex=="No_BPD_AMAB")
exp2_combined_16s_pprd_afab = subset_samples(exp2_combined_16s_pprd, BPD_Sex=="No_BPD_AFAB")


exp2_combined_its = subset_taxa(exp2_combined_30_gen, Kingdom=="Fungi")

exp2_combined_its_bpd = subset_samples(exp2_combined_its, BPD=="BPD")
exp2_combined_its_bpd_amab = subset_samples(exp2_combined_its_bpd, BPD_Sex=="BPD_AMAB")
exp2_combined_its_bpd_afab = subset_samples(exp2_combined_its_bpd, BPD_Sex=="BPD_AFAB")

exp2_combined_its_pprd = subset_samples(exp2_combined_its, BPD=="No_BPD")
exp2_combined_its_pprd_amab = subset_samples(exp2_combined_its_pprd, BPD_Sex=="No_BPD_AMAB")
exp2_combined_its_pprd_afab = subset_samples(exp2_combined_its_pprd, BPD_Sex=="No_BPD_AFAB")

###Run SPIEC-EASI###
set.seed(1312)
se.exp2.gen.pprd.amab <- spiec.easi(list(exp2_combined_16s_pprd_amab, exp2_combined_its_pprd_amab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))
spiec.graph.gen.pprd.amab=adj2igraph(getRefit(se.exp2.gen.pprd.amab), vertex.attr=list(name=taxa_names(exp2_combined_30_gen_noNA)))

#Plot network without labels
prev_gen_spiec_plot_pprd_amab_nolegend <- plot_network(spiec.graph.gen.pprd.amab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_amab_nolegend + scale_color_manual(values = spiec.colors)
#Save

#Plot network with labels
prev_gen_spiec_plot_pprd_amab_legend <- plot_network(spiec.graph.gen.pprd.amab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

#Get relevant stats
nodes_pprd_amab <- gorder(spiec.graph.gen.pprd.amab)
edges_pprd_amab <- gsize(spiec.graph.gen.pprd.amab)


############################Figure 3G: PPRD AFAB SPIEC-EASI NETWORK############################

#Run SPIEC-EASI
se.exp2.gen.pprd.afab <- spiec.easi(list(exp2_combined_16s_pprd_afab, exp2_combined_its_pprd_afab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

spiec.graph.gen.pprd.afab=adj2igraph(getRefit(se.exp2.gen.pprd.afab), vertex.attr=list(name=taxa_names(exp2_combined_30_gen_noNA)))

#Plot network without labels
prev_gen_spiec_plot_pprd_afab_nolegend <- plot_network(spiec.graph.gen.pprd.afab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_afab_nolegend + scale_color_manual(values = spiec.colors)
#Save

#Plot network with labels
prev_gen_spiec_plot_pprd_afab_legend <- plot_network(spiec.graph.gen.pprd.afab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
ggsave("Fig 3G without NAs.tiff", width = 10, height = 10, device='tiff', dpi=600)
#Save

nodes_pprd_afab <- gorder(spiec.graph.gen.pprd.afab)
edges_pprd_afab <- gsize(spiec.graph.gen.pprd.afab)


############################Figure 3H: BPD AMAB SPIEC-EASI NETWORK############################

set.seed(1312)
se.exp2.gen.bpd.amab <- spiec.easi(list(exp2_combined_16s_bpd_amab, exp2_combined_its_bpd_amab), method='mb', nlambda=99,
                                   lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))
spiec.graph.gen.bpd.amab=adj2igraph(getRefit(se.exp2.gen.bpd.amab), vertex.attr=list(name=taxa_names(exp2_combined_30_gen_noNA)))

#Plot network without labels
prev_gen_spiec_plot_bpd_amab_nolegend <- plot_network(spiec.graph.gen.bpd.amab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_amab_nolegend + scale_color_manual(values = spiec.colors)
#Save

#Plot network with labels
prev_gen_spiec_plot_bpd_amab_legend <- plot_network(spiec.graph.gen.bpd.amab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

nodes_bpd_amab <- gorder(spiec.graph.gen.bpd.amab)
edges_bpd_amab <- gsize(spiec.graph.gen.bpd.amab)


############################Figure 3I: BPD AFAB SPIEC-EASI NETWORK############################

se.exp2.gen.bpd.afab <- spiec.easi(list(exp2_combined_16s_bpd_afab, exp2_combined_its_bpd_afab), method='mb', nlambda=99,
                                   lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))
spiec.graph.gen.bpd.afab=adj2igraph(getRefit(se.exp2.gen.bpd.afab), vertex.attr=list(name=taxa_names(exp2_combined_30_gen_noNA)))


#Plot network without labels
prev_gen_spiec_plot_bpd_afab_nolegend <- plot_network(spiec.graph.gen.bpd.afab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_afab_nolegend + scale_color_manual(values = spiec.colors)
#Save

#Plot network with labels
prev_gen_spiec_plot_bpd_afab_legend <- plot_network(spiec.graph.gen.bpd.afab, exp2_combined_30_gen_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

nodes_bpd_afab <- gorder(spiec.graph.gen.bpd.afab)
edges_bpd_afab <- gsize(spiec.graph.gen.bpd.afab)
