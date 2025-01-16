library(tidyverse)
library(phyloseq)
library(SpiecEasi)
library(igraph)
library(DESeq2)


############################FIGURE S6A: ITS BPD-SAAB DESEQ2 SIGNIFICANT BOXPLOT############################

###Run DESeq2###
diff = phyloseq_to_deseq2(exp2_fung_gen_prev, ~ BPD_Sex)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")

###Pull out all pairwise comparisons###
deseq2resultsBPD_Sex1 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_AFAB","BPD_AMAB"))
deseq2resultsBPD_Sex2 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_AFAB","No_BPD_AFAB"))
deseq2resultsBPD_Sex3 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_AFAB","No_BPD_AMAB"))
deseq2resultsBPD_Sex4 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_AMAB","No_BPD_AFAB"))
deseq2resultsBPD_Sex5 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","BPD_AMAB","No_BPD_AMAB"))
deseq2resultsBPD_Sex6 = results(diff, pAdjustMethod = "fdr",contrast = c("BPD_Sex","No_BPD_AFAB","No_BPD_AMAB"))

###This function creates a data frame with the signficantly differentially abundant genera###
diffabund_deseq2_gen <- function(deseq2results){
  deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
  deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp2_fung_gen_prev)[rownames(deseq2results), ], "matrix"))
  sigtab = deseq2results[(deseq2results$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp2_fung_gen_prev)[rownames(sigtab), ], "matrix"))
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
exp2_fung_rel_gen_prev <- transform_sample_counts(exp2_fung_gen_prev, function(x) x / sum(x) )
exp2_fung_signif_genus <- subset_taxa(exp2_fung_rel_gen_prev, rownames(otu_table(exp2_fung_rel_gen_prev)) %in% c("ASV18049","ASV18054","ASV18071","ASV18080","ASV18090"))
exp2_fung_signif_otu_genus <- as.data.frame(t(exp2_fung_signif_genus@otu_table))
exp2_fung_signif_meta_genus <- as.data.frame(exp2_fung_signif_genus@sam_data)
colnames(exp2_fung_signif_otu_genus) <- c("ASV18049: Candida","ASV18054: Plectosphaerella","ASV18071: Gibellulopsis","ASV18080: Pseudeurotium","ASV18090: Mortierella")

###Add BPD data###
exp2_fung_diffabund_signif_genus <- exp2_fung_signif_otu_genus %>%
  dplyr::mutate(BPD = exp2_fung_signif_meta_genus$BPD,Sex = exp2_fung_signif_meta_genus$Sex, BPD_Sex = exp2_fung_signif_meta_genus$BPD_Sex)
#Save

exp2_fung_signif_genus_df <- psmelt(exp2_fung_signif_genus)

###Plot###
diffabund_boxplot <- ggplot(data = exp2_fung_signif_genus_df, aes(x = BPD_Sex, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = BPD), height = 0, width = .2) +
  labs(x = "", y = "Abundance" ) +
  facet_wrap(~ Genus, scales = "free", labeller = labeller(OTU=names) ) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))

diffabund_boxplot + scale_color_manual(values = c("#007016","#94D57F"))
#Save


############################FIGURE S6B: ITS BPD-SAAB DESEQ2 LOG2 FOLD CHANGE BARPLOT############################

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

###BPD_AFAB vs No_BPD_AMAB###
deseq2plot(deseq2BPD_Sex3)
#Save 

###BPD_AMAB vs No_BPD_AFAB###
deseq2plot(deseq2BPD_Sex4)
#Save

###BPD_AMAB vs No_BPD_AMAB###
deseq2plot(deseq2BPD_Sex5)
#Save

###No_BPD_AFAB vs No_BPD_AMAB###
deseq2plot(deseq2BPD_Sex6)
#Save


############################FIGURE S6C: 16S BPD-SAAB SPIEC-EASI PLOTS############################

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

###No_BPD_AFAB Stats###
nodes_No_BPD_afab <- gorder(spiec.graph.gen.pprd.afab)
edges_No_BPD_afab <- gsize(spiec.graph.gen.pprd.afab)
betweenness_No_BPD_afab <- as.list(betweenness(spiec.graph.gen.pprd.afab, normalized = TRUE))
degree_No_BPD_afab <- as.data.frame(degree(spiec.graph.gen.pprd.afab))

###No_BPD_AMAB Stats###
nodes_No_BPD_amab <- gorder(spiec.graph.gen.pprd.amab)
edges_No_BPD_amab <- gsize(spiec.graph.gen.pprd.amab)
betweenness_No_BPD_amab <- as.list(betweenness(spiec.graph.gen.pprd.amab, normalized = TRUE))
degree_No_BPD_amab <- as.data.frame(degree(spiec.graph.gen.pprd.amab))


###ALL TAXA EDGE PLOT###
edges_df_bpdsex <- data.frame(edges_bpd_afab,edges_bpd_amab,edges_pprd_afab,edges_pprd_amab)
colnames(edges_df_bpdsex) <- c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB")
edges_long <- edges_df_bpdsex %>%
  pivot_longer(c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB"), names_to = "BPD_Sex", values_to = "Edges")

ggplot(edges_long,aes(x= BPD_Sex, y = Edges))+
  geom_bar(stat="identity", aes(fill=BPD_Sex))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save


###ALL TAXA NODE PLOT###

nodes_df_bpdsex <- data.frame(nodes_bpd_afab,nodes_bpd_amab,nodes_pprd_afab,nodes_pprd_amab)
colnames(nodes_df_bpdsex) <- c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB")
nodes_long <- nodes_df_bpdsex %>%
  pivot_longer(c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB"), names_to = "BPD_Sex", values_to = "Nodes")

ggplot(nodes_long,aes(x= BPD_Sex, y = Nodes))+
  geom_bar(stat="identity", aes(fill=BPD_Sex))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save


############################FIGURE S6D: FUNGI SUBNETWORK BPD-SAAB SPIEC-EASI PLOTS############################

#Function for extracting subnetworks written by Dr. Tipton
extract_surnodes <- function(selntwrk,node){
  refit <- getOptNet(selntwrk)
  assoc <- which(refit[,node[1]]==1)
  assoc <- c(node, assoc)
  if(length(node) > 1){
    for(nod in 2:length(node)){
      tempassoc <- which(refit[,node[nod]]==1)
      assoc <- unique(c(assoc, tempassoc))
    }
  }
  refit2 <- refit[assoc,assoc]
  colnames(refit2) = rownames(refit2) <- assoc
  return(refit2)
  rm(assoc)
}

###Get Subnetwork of all Fungi and all adjacent nodes for BPD network###

###BPD AFAB###
#Extract all nodes connected to all fungal ASVs
fungi_bpd_afab_subnetwork <- extract_surnodes(se.exp2.gen.bpd.afab, c(51:69))
fungi.bpd.afab.graphable <- as.matrix(fungi_bpd_afab_subnetwork)

#Label ASVs, get metrics like main network
fungi.bpd.afab.subnet <- graph.adjacency(fungi.bpd.afab.graphable, mode="undirected", weighted=NULL)
nodes_fungi_bpd_afab <- gorder(fungi.bpd.afab.subnet)
edges_fungi_bpd_afab <- gsize(fungi.bpd.afab.subnet)
degree_fungi_bpd_afab <- as.data.frame(degree(fungi.bpd.afab.subnet))



###BPD AMAB###
fungi_bpd_amab_subnetwork <- extract_surnodes(se.exp2.gen.bpd.amab, c(51:69))
fungi.bpd.amab.graphable <- as.matrix(fungi_bpd_amab_subnetwork)
row.names(fungi.bpd.amab.graphable)

#Label ASVs, get metrics like main network
fungi.bpd.amab.subnet <- graph.adjacency(fungi.bpd.amab.graphable, mode="undirected", weighted=NULL)
nodes_fungi_bpd_amab <- gorder(fungi.bpd.amab.subnet)
edges_fungi_bpd_amab <- gsize(fungi.bpd.amab.subnet)
degree_fungi_bpd_amab <- as.data.frame(degree(fungi.bpd.amab.subnet))


###NO BPD AFAB###
#Extract all nodes connected to all fungal ASVs
fungi_pprd_afab_subnetwork <- extract_surnodes(se.exp2.gen.pprd.afab, c(51:69))
fungi.pprd.afab.graphable <- as.matrix(fungi_pprd_afab_subnetwork)
row.names(fungi.pprd.afab.graphable)

#Label ASVs, get metrics like main network
fungi.pprd.afab.subnet <- graph.adjacency(fungi.pprd.afab.graphable, mode="undirected", weighted=NULL)
nodes_fungi_pprd_afab <- gorder(fungi.pprd.afab.subnet)
edges_fungi_pprd_afab <- gsize(fungi.pprd.afab.subnet)
degree_fungi_pprd_afab <- as.data.frame(degree(fungi.pprd.afab.subnet))



###NO BPD AMAB###
fungi_pprd_amab_subnetwork <- extract_surnodes(se.exp2.gen.pprd.amab, c(51:69))
fungi.pprd.amab.graphable <- as.matrix(fungi_pprd_amab_subnetwork)
row.names(fungi.pprd.amab.graphable)

#Label ASVs, get metrics like main network
fungi.pprd.amab.subnet <- graph.adjacency(fungi.pprd.amab.graphable, mode="undirected", weighted=NULL)
nodes_fungi_pprd_amab <- gorder(fungi.pprd.amab.subnet)
edges_fungi_pprd_amab <- gsize(fungi.pprd.amab.subnet)
degree_fungi_pprd_amab <- as.data.frame(degree(fungi.pprd.amab.subnet))

###FUNGI SUBNET EDGE PLOT###
edges_fungi_df_bpdsex <- data.frame(edges_fungi_bpd_afab,edges_fungi_bpd_amab,edges_fungi_pprd_afab,edges_fungi_pprd_amab)
colnames(edges_fungi_df_bpdsex) <- c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB")
edges_fungi_long <- edges_fungi_df_bpdsex %>%
  pivot_longer(c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB"), names_to = "BPD_Sex", values_to = "Edges")

ggplot(edges_fungi_long,aes(x= BPD_Sex, y = Edges))+
  geom_bar(stat="identity", aes(fill=BPD_Sex))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save


###FUNGI SUBNET NODE PLOT###
nodes_fungi_df_bpdsex <- data.frame(nodes_fungi_bpd_afab,nodes_fungi_bpd_amab,nodes_fungi_pprd_afab,nodes_fungi_pprd_amab)
colnames(nodes_fungi_df_bpdsex) <- c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB")
nodes_fungi_long <- nodes_fungi_df_bpdsex %>%
  pivot_longer(c("BPD_AFAB","BPD_AMAB","No_BPD_AFAB","No_BPD_AMAB"), names_to = "BPD_Sex", values_to = "Nodes")

ggplot(nodes_fungi_long,aes(x= BPD_Sex, y = Nodes))+
  geom_bar(stat="identity", aes(fill=BPD_Sex))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save


############################FIGURE S6E: 16S NOBPD-SAAB SPIEC-EASI NETWORKS WITH LABELS############################

###No_BPD_AFAB###
prev_gen_spiec_plot_No_BPD_afab_legend <- plot_network(spiec.graph.gen.pprd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_No_BPD_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

###No_BPD_AMAB###
prev_gen_spiec_plot_No_BPD_amab_legend <- plot_network(spiec.graph.gen.pprd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_No_BPD_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save


############################FIGURE S6F: 16S BPD-SAAB SPIEC-EASI NETWORKS WITH LABELS############################

###BPD_AFAB###
prev_gen_spiec_plot_bpd_afab_legend <- plot_network(spiec.graph.gen.bpd.afab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_afab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

###BPD_AMAB###
prev_gen_spiec_plot_bpd_amab_legend <- plot_network(spiec.graph.gen.bpd.amab, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_amab_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save
