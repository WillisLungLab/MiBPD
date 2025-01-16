library(phyloseq)
library(SpiecEasi)
library(plyr)
library(tidyverse)
library(readxl)
library(vegan)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(factoextra)


############################FIGURE S2A: ITS DESEQ2 BOXPLOT############################

###Run DESeq2###
diff = phyloseq_to_deseq2(exp2_fung_gen_prev, ~ BPD)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")
deseq2results = results(diff, pAdjustMethod = "fdr")
deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp2_fung_gen_prev)[rownames(deseq2results), ], "matrix"))
sigtab = deseq2results[(deseq2results$padj < 0.05), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp2_fung_gen_prev)[rownames(sigtab), ], "matrix"))
sigtabgenus = subset(sigtab, !is.na(Genus))
print(sigtabgenus)
#Save significance table as .csv

###Make table of relative abundance of significantly differentially abundant genera###
exp2_fung_rel_gen_prev <- transform_sample_counts(exp2_fung_gen_prev, function(x) x / sum(x) )

exp2_fung_signif_genus <- subset_taxa(exp2_fung_rel_gen_prev, rownames(otu_table(exp2_fung_rel_gen_prev)) %in% rownames(sigtabgenus))
exp2_fung_signif_otu_genus <- as.data.frame(t(exp2_fung_signif_genus@otu_table))
exp2_fung_signif_meta_genus <- as.data.frame(exp2_fung_signif_genus@sam_data)
colnames(exp2_fung_signif_otu_genus) <- c("Candida","Mortierella")

###Add BPD data###
exp2_fung_diffabund_signif_genus <- exp2_fung_signif_otu_genus %>%
  dplyr::mutate(BPD = exp2_fung_signif_meta_genus$BPD)
#Save .csv and make boxplot in GraphPad

exp2_pankaj_diffabund_df <- psmelt(exp2_fung_signif_genus)

diffabund_boxplot <- ggplot(data = exp2_pankaj_diffabund_df, aes(x = BPD, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = BPD), height = 0, width = .2) +
  labs(x = "", y = "Abundance" ) +
  facet_wrap(~ Genus, scales = "free") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))

diffabund_boxplot
#Save


############################FIGURE S2B: ITS DESEQ2 SIGNIFICANCE BARPLOT############################

###Creates barplot of log2 fold change of significant genera###
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
#Save image


############################FIGURE S2C: 16s DESEQ2 BOXPLOT############################

###Run DESeq2###
diff = phyloseq_to_deseq2(exp2_bact_gen_prev, ~ BPD)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")
deseq2results = results(diff, pAdjustMethod = "fdr")
deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp2_bact_gen_prev)[rownames(deseq2results), ], "matrix"))
sigtab = deseq2results[(deseq2results$padj < 0.05), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp2_bact_gen_prev)[rownames(sigtab), ], "matrix"))
sigtabgenus = subset(sigtab, !is.na(Genus))
print(sigtabgenus)
#Save significance table as .csv

###Make table of relative abundance of significantly differentially abundant genera###
exp2_bact_rel_gen_prev <- transform_sample_counts(exp2_bact_gen_prev, function(x) x / sum(x) )

exp2_bact_signif_genus <- subset_taxa(exp2_bact_rel_gen_prev, rownames(otu_table(exp2_bact_rel_gen_prev)) %in% rownames(sigtabgenus))
exp2_bact_signif_otu_genus <- as.data.frame(t(exp2_bact_signif_genus@otu_table))
exp2_bact_signif_meta_genus <- as.data.frame(exp2_bact_signif_genus@sam_data)
colnames(exp2_bact_signif_otu_genus) <- c("Veillonella","Unidentified Staphylococcaceae","Cutibacterium")

###Add BPD data###
exp2_bact_diffabund_signif_genus <- exp2_bact_signif_otu_genus %>%
  dplyr::mutate(BPD = exp2_bact_signif_meta_genus$BPD)

#Save .csv and make boxplot in GraphPad

exp2_pankaj_diffabund_df <- psmelt(exp2_bact_signif_genus)

diffabund_boxplot <- ggplot(data = exp2_pankaj_diffabund_df, aes(x = BPD, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = BPD), height = 0, width = .2) +
  labs(x = "", y = "Abundance" ) +
  facet_wrap(~ Genus, scales = "free") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))

diffabund_boxplot
#Save


############################FIGURE S2D: 16S DESEQ2 SIGNIFICANCE BARPLOT############################

###Creates barplot of log2 fold change of significant genera###
deseq2_colors <- c("#94D57F")
deseq2plot <- function(sigtabgenus){
  diffabund <- ggplot(sigtabgenus,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
    geom_col() +
    coord_flip() +
    scale_fill_manual(values=deseq2_colors, labels=c("positive"))
  diffabund2 <- diffabund + theme_bw()
  diffabund2
}
sigtabgenus_real <- sigtabgenus[c(1,2),]
deseq2plot(sigtabgenus_real)
#Save image


############################FIGURE S2D: 16S SPIEC-EASI EDGE PLOT############################

###Use SPIEC-EASI data from Fig. 1H###

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

#Extract all nodes connected to all fungal ASVs
fungi_bpd_subnetwork <- extract_surnodes(se.exp2.gen.bpd, c(51:69))
fungi.bpd.graphable <- as.matrix(fungi_bpd_subnetwork)

#Get metrics like main network
fungi.bpd.subnet <- graph.adjacency(fungi.bpd.graphable, mode="undirected", weighted=NULL)
nodes_fungi_bpd <- gorder(fungi.bpd.subnet)
edges_fungi_bpd <- gsize(fungi.bpd.subnet)


###Get Subnetwork of all Fungi and all adjacent nodes for NoBPD network###

fungi_pprd_subnetwork <- extract_surnodes(se.exp2.gen.pprd, c(51:69))
fungi.pprd.graphable <- as.matrix(fungi_pprd_subnetwork)

#Get metrics
fungi.pprd.subnet <- graph.adjacency(fungi.pprd.graphable, mode="undirected", weighted=NULL)
nodes_fungi_pprd <- gorder(fungi.pprd.subnet)
edges_fungi_pprd <- gsize(fungi.pprd.subnet)



###EDGE PLOT WITH MAIN NETWORK AND FUNGAL SUBNETWORKS###
#Make data frame with edge count for fungi and connected bacteria
edges_df <- data.frame(edges_bpd,edges_pprd)
colnames(edges_df) <- c("BPD_All","No_BPD_All")
edges_fungi_df <- data.frame(edges_fungi_bpd,edges_fungi_pprd)
colnames(edges_fungi_df) <- c("BPD_Fungi","No_BPD_Fungi")

#Combine with edge counts for overall plot
edges_df_all <- cbind(edges_df,edges_fungi_df)
edges_long_all <- edges_df_all %>%
  pivot_longer(c("BPD_All","No_BPD_All","BPD_Fungi","No_BPD_Fungi"), names_to = "BPD", values_to = "Edges")
col_order <- c("BPD_All","No_BPD_All","BPD_Fungi","No_BPD_Fungi")
#Plot both
ggplot(edges_long_all,aes(x= factor(BPD, col_order), y = Edges))+
  geom_bar(stat="identity", aes(fill=BPD))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  labs(x = "BPD",
       y = "Edges") + 
  theme_bw()
#Save 


############################FIGURE S2E: 16S SPIEC-EASI NODE PLOT############################

nodes_df <- data.frame(nodes_bpd,nodes_pprd)
colnames(nodes_df) <- c("BPD_All","No_BPD_All")
nodes_fungi_df <- data.frame(nodes_fungi_bpd,nodes_fungi_pprd)
colnames(nodes_fungi_df) <- c("BPD_Fungi","No_BPD_Fungi")

#Combine with edge counts for overall plot
nodes_df_all <- cbind(nodes_df,nodes_fungi_df)
nodes_long_all <- nodes_df_all %>%
  pivot_longer(c("BPD_All","No_BPD_All","BPD_Fungi","No_BPD_Fungi"), names_to = "BPD", values_to = "nodes")
col_order <- c("BPD_All","No_BPD_All","BPD_Fungi","No_BPD_Fungi")
#Plot both
ggplot(nodes_long_all,aes(x= factor(BPD, col_order), y = nodes))+
  geom_bar(stat="identity", aes(fill=BPD))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  labs(x = "BPD",
       y = "nodes") + 
  theme_bw()
#Save 

ggsave("Fig E1E.tiff", width = 5, height = 5, device='tiff', dpi=600)


############################FIGURE S2G: 16S NOBPD SPIEC-EASI FULL NETWORK WITH TAXA LABELS############################

###Use SPIEC-EASI Network from Fig. 1G###

#NoBPD only network plot with taxa labels
prev_gen_spiec_plot_pprd <- plot_network(spiec.graph.gen.pprd, exp2_bact_gen_prev_noNA, type='taxa',label=NULL,color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save


############################FIGURE S2F: 16S BPD SPIEC-EASI FULL NETWORK WITH TAXA LABELS############################

###Use SPIEC-EASI Network from Fig. 1H###

#BPD only network plot with taxa labels
prev_gen_spiec_plot_bpd_legend <- plot_network(spiec.graph.gen.bpd, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save
