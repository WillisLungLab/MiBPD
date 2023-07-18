library(tidyverse)
library(phyloseq)
library(SpiecEasi)
library(igraph)

############################FIGURE E7A: MULTIKINGDOM BPD-SAAB BARPLOT############################

###Merge all samples in each group and aggregate to Genus level###
exp2_combined_merged_bpdsex = merge_samples(exp2_combined_prev, "BPD_Sex")
exp2_combined_gen_merged_bpdsex <- tax_glom(exp2_combined_merged_bpdsex, taxrank = 'Genus')

###Get top 19 genera and merge all other genera into the "Other" category###
top20_combined_merged_bpdsex_list <- names(sort(taxa_sums(exp2_combined_gen_merged_bpdsex), decreasing=TRUE)[1:19])
top20_merged_bpdsex_combined_rel <- transform_sample_counts(exp2_combined_gen_merged_bpdsex, function(x) x / sum(x) )
top20_merged_bpdsex_combined_df <- psmelt(top20_merged_bpdsex_combined_rel)
top20_merged_bpdsex_combined_df[!(top20_merged_bpdsex_combined_df$OTU %in% top20_combined_merged_bpdsex_list),]$Genus <- 'Other'

###Plot barplot###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(20)
barplot_gen_bpd_combined <- ggplot(top20_merged_bpdsex_combined_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_combined + scale_fill_manual(values = barplot_colors)


############################FIGURE E7B: Multikingdom BPD-SAAB DESEQ2 SIGNIFICANT BOXPLOT############################

###Run DESeq2###
diff = phyloseq_to_deseq2(exp2_combined_gen_prev, ~ BPD_Sex)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")

###Pull out all pairwise comparisons###
deseq2resultsBPD_Sex1 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_F","BPD_M"))
deseq2resultsBPD_Sex2 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_F","PPRD_F"))
deseq2resultsBPD_Sex3 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_F","PPRD_M"))
deseq2resultsBPD_Sex4 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_M","PPRD_F"))
deseq2resultsBPD_Sex5 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_M","PPRD_M"))
deseq2resultsBPD_Sex6 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","PPRD_F","PPRD_M"))

###This function creates a data frame with the signficantly differentially abundant genera###
diffabund_deseq2_gen <- function(deseq2results){
  deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
  deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp2_combined_gen_prev)[rownames(deseq2results), ], "matrix"))
  sigtab = deseq2results[(deseq2results$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp2_combined_gen_prev)[rownames(sigtab), ], "matrix"))
  sigtabgenus = subset(sigtab, !is.na(Genus))
  print(sigtabgenus)
}

###Get significantly differentiall abundant genera for each pairwise comparison
deseq2BPD_Sex1 <- diffabund_deseq2_gen(deseq2resultsBPD_Sex1)
deseq2BPD_Sex2 <- diffabund_deseq2_gen(deseq2resultsBPD_Sex2)
deseq2BPD_Sex3 <- diffabund_deseq2_gen(deseq2resultsBPD_Sex3)
deseq2BPD_Sex4 <- diffabund_deseq2_gen(deseq2resultsBPD_Sex4)
deseq2BPD_Sex5 <- diffabund_deseq2_gen(deseq2resultsBPD_Sex5)
deseq2BPD_Sex6 <- diffabund_deseq2_gen(deseq2resultsBPD_Sex6)

###Make dataframe of selected significant taxa###
exp2_combined_rel_gen_prev <- transform_sample_counts(exp2_combined_gen_prev, function(x) x / sum(x) )
exp2_combined_signif_genus <- subset_taxa(exp2_combined_rel_gen_prev, rownames(otu_table(exp2_combined_rel_gen_prev)) %in% c("ASV3","ASV7","ASV25","ASV45","ASV72","ASV175","ASV18049","ASV18069","ASV18071"))
exp2_combined_signif_genus_df <- psmelt(exp2_combined_signif_genus)

###Plot###
diffabund_boxplot <- ggplot(data = exp2_combined_signif_genus_df, aes(x = BPD_Sex, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = BPD), height = 0, width = .2) +
  labs(x = "", y = "Abundance" ) +
  facet_wrap(~ Genus, scales = "free", labeller = labeller(OTU=names) ) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))

diffabund_boxplot + scale_color_manual(values = c("#007016","#94D57F"))
#Save

############################FIGURE E7C: Multikingdom BPD-SAAB DESEQ2 LOG2 FOLD CHANGE BARPLOT############################

###This function plots the significant taxa for each comparison###
deseq2_colors <- c("#007016","#94D57F")
deseq2plot <- function(sigtabgenus){
  diffabund <- ggplot(sigtabgenus,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
    geom_col() +
    coord_flip() +
    scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
  diffabund2 <- diffabund + theme_bw()
  diffabund2
}

###BPD_AMAB vs PPRD_AMAB###
deseq2plot(deseq2BPD_Sex5)

###PPRD_AFAB vs PPRD_AMAB###
deseq2plot(deseq2BPD_Sex6)

###BPD_AFAB vs PPRD_AMAB###
deseq2plot(deseq2BPD_Sex3)

###BPD_AFAB vs BPD_AMAB###
deseq2plot(deseq2BPD_Sex1)

#Save all images

############################FIGURE E7D: Multikingdom BPD-SAAB SPIEC-EASI EDGE PLOTS############################

###Filter out ASVs in <30% of samples since SPIEC-EASI can take a long time or return errors with too many ASVs###
exp2_combined_50_gen <- filter_taxa(exp2_combined_gen, function(x) sum(x >= 1) > (0.30*length(x)), TRUE)
exp2_combined_gen_prev_noNA = subset_taxa(exp2_combined_50_gen, Genus!="NA")

exp2_prev_gen_otu <- as.data.frame(exp2_combined_gen_prev_noNA@otu_table)
exp2_prev_gen_tax <- as.data.frame(exp2_combined_gen_prev_noNA@tax_table)
exp2_prev_gen_meta <- as.data.frame(exp2_combined_gen_prev_noNA@sam_data)

###Split up 16s and ITS###
exp2_prev_gen_bact = subset_taxa(exp2_combined_gen_prev_noNA, Kingdom == "Bacteria")
exp2_prev_gen_fung = subset_taxa(exp2_combined_gen_prev_noNA, Kingdom == "Fungi")

###Split 16s by BPD status and SAAB###
exp2_prev_gen_bact_bpd_afab = subset_samples(exp2_prev_gen_bact, BPD_Sex == "BPD_AFAB")
exp2_prev_gen_bact_bpd_amab = subset_samples(exp2_prev_gen_bact, BPD_Sex == "BPD_AMAB")
exp2_prev_gen_bact_pprd_afab = subset_samples(exp2_prev_gen_bact, BPD_Sex == "PPRD_AFAB")
exp2_prev_gen_bact_pprd_amab = subset_samples(exp2_prev_gen_bact, BPD_Sex == "PPRD_AMAB")

###Split ITS by BPD status and SAAB###
exp2_prev_gen_fung_bpd_afab = subset_samples(exp2_prev_gen_fung, BPD_Sex == "BPD_AFAB")
exp2_prev_gen_fung_bpd_amab = subset_samples(exp2_prev_gen_fung, BPD_Sex == "BPD_AMAB")
exp2_prev_gen_fung_pprd_afab = subset_samples(exp2_prev_gen_fung, BPD_Sex == "PPRD_AFAB")
exp2_prev_gen_fung_pprd_amab = subset_samples(exp2_prev_gen_fung, BPD_Sex == "PPRD_AMAB")

spiec.colors <- c("#aaaaaa","#74c476")

###SpiecEasi###
set.seed(1312)
se.exp2.gen.bpd.afab <- spiec.easi(list(exp2_prev_bact_bpd_afab, exp2_prev_fung_bpd_afab), method='mb', nlambda=99,
                                   lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

se.exp2.gen.bpd.amab <- spiec.easi(list(exp2_prev_bact_bpd_amab, exp2_prev_fung_bpd_amab), method='mb', nlambda=99,
                                   lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

se.exp2.gen.pprd.afab <- spiec.easi(list(exp2_prev_bact_pprd_afab, exp2_prev_fung_pprd_afab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

se.exp2.gen.pprd.amab <- spiec.easi(list(exp2_prev_bact_pprd_amab, exp2_prev_fung_pprd_amab), method='mb', nlambda=99,
                                    lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))


###BPD_AFAB Stats###
nodes_bpd_afab <- gorder(spiec.graph.gen.bpd.afab)
edges_bpd_afab <- gsize(spiec.graph.gen.bpd.afab)
betweenness_bpd_afab <- as.list(betweenness(spiec.graph.gen.bpd.afab, normalized = TRUE))
degree_bpd_afab <- as.data.frame(degree(spiec.graph.gen.bpd.afab))

###BPD_AMAB Stats###
nodes_bpd_amab <- gorder(spiec.graph.gen.bpd.amab)
edges_bpd_amab <- gsize(spiec.graph.gen.bpd.amab)
betweenness_bpd_amab <- as.list(betweenness(spiec.graph.gen.bpd.amab, normalized = TRUE))
degree_bpd_amab <- as.data.frame(degree(spiec.graph.gen.bpd.amab))

###PPRD_AFAB Stats###
nodes_pprd_afab <- gorder(spiec.graph.gen.pprd.afab)
edges_pprd_afab <- gsize(spiec.graph.gen.pprd.afab)
betweenness_pprd_afab <- as.list(betweenness(spiec.graph.gen.pprd.afab, normalized = TRUE))
degree_pprd_afab <- as.data.frame(degree(spiec.graph.gen.pprd.afab))

###PPRD_AMAB Stats###
nodes_pprd_amab <- gorder(spiec.graph.gen.pprd.amab)
edges_pprd_amab <- gsize(spiec.graph.gen.pprd.amab)
betweenness_pprd_amab <- as.list(betweenness(spiec.graph.gen.pprd.amab, normalized = TRUE))
degree_pprd_amab <- as.data.frame(degree(spiec.graph.gen.pprd.amab))

###Edges Plot###
edges_df_bpdsex <- data.frame(edges_bpd_afab,edges_bpd_amab,edges_pprd_afab,edges_pprd_amab)
colnames(edges_df_bpdsex) <- c("BPD_AFAB","BPD_AMAB","PPRD_AFAB","PPRD_AMAB")
edges_long <- edges_df_bpdsex %>%
  pivot_longer(c("BPD_AFAB","BPD_AMAB","PPRD_AFAB","PPRD_AMAB"), names_to = "BPD_Sex", values_to = "Edges")

ggplot(edges_long,aes(x= BPD_Sex, y = Edges))+
  geom_bar(stat="identity", aes(fill=BPD_Sex))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save

############################FIGURE E7E: Multikingdom BPD-SAAB SPIEC-EASI NODE PLOTS############################

nodes_df_bpdsex <- data.frame(nodes_bpd_afab,nodes_bpd_amab,nodes_pprd_afab,nodes_pprd_amab)
colnames(nodes_df_bpdsex) <- c("BPD_AFAB","BPD_AMAB","PPRD_AFAB","PPRD_AMAB")
nodes_long <- nodes_df_bpdsex %>%
  pivot_longer(c("BPD_AFAB","BPD_AMAB","PPRD_AFAB","PPRD_AMAB"), names_to = "BPD_Sex", values_to = "Nodes")

ggplot(nodes_long,aes(x= BPD_Sex, y = Nodes))+
  geom_bar(stat="identity", aes(fill=BPD_Sex))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save


############################FIGURE E7F: Multikingdom PPRD-SAAB SPIEC-EASI NETWORKS WITH LABELS############################

###PPRD_AFAB###
prev_gen_spiec_plot_pprd_afab_legend <- plot_network(spiec.graph.gen.pprd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

###PPRD_AMAB###
prev_gen_spiec_plot_pprd_amab_legend <- plot_network(spiec.graph.gen.pprd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

############################FIGURE E7G: Multikingdom BPD-SAAB SPIEC-EASI NETWORKS WITH LABELS############################

###BPD_AFAB###
prev_gen_spiec_plot_bpd_afab_legend <- plot_network(spiec.graph.gen.bpd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

###BPD_AMAB###
prev_gen_spiec_plot_bpd_amab_legend <- plot_network(spiec.graph.gen.bpd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

