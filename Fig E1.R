library(phyloseq)
library(SpiecEasi)
library(plyr)
library(tidyverse)
library(readxl)
library(vegan)
library(DESeq2)
library(indicspecies)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(factoextra)

############################FIGURE E1A: MULTIKINGDOM TOP 20 BARPLOT############################

###Get a list of 19 most prevalent genera and renames all other genera as "Other"###
exp2_combined_merged = merge_samples(exp2_combined_prev, "BPD")
exp2_combined_gen_merged <- tax_glom(exp2_combined_merged, taxrank = 'Genus')

top20_combined_merged_list <- names(sort(taxa_sums(exp2_combined_gen_merged), decreasing=TRUE)[1:19])
top20_merged_combined_rel <- transform_sample_counts(exp2_combined_gen_merged, function(x) x / sum(x) )
top20_merged_combined_df <- psmelt(top20_merged_combined_rel)
top20_merged_combined_df[!(top20_merged_combined_df$OTU %in% top20_combined_merged_list),]$Genus <- 'Other'


###Barplot of top 20 genera, including the "Other" category###

barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(20)

barplot_gen_bpd_its <- ggplot(top20_merged_combined_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_its + scale_fill_manual(values = barplot_colors)
#Save


############################FIGURE E1B: ITS DESEQ2 BOXPLOT############################

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
colnames(exp2_fung_signif_otu_genus) <- c("Candida","Saccharomyces")

###Add BPD data###
exp2_fung_diffabund_signif_genus <- exp2_fung_signif_otu_genus %>%
  dplyr::mutate(BPD = exp2_fung_signif_meta_genus$BPD)

#Save .csv and make boxplot in GraphPad


############################FIGURE E1C: ITS DESEQ2 SIGNIFICANCE BARPLOT############################

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


############################FIGURE E1D: Multikingdom SPIEC-EASI EDGE PLOT############################

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
fungi_bpd_subnetwork <- extract_surnodes(se.exp2.gen.bpd, c(61:80))
fungi.bpd.graphable <- as.matrix(fungi_bpd_subnetwork)

#Label ASVs, get metrics like main network
colnames(fungi.bpd.graphable) <- row.names(exp2_prev_gen_otu[c(61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,12,49,31,4,56,47,48,15,11,59,40,46,43),])
fungi.bpd.subnet <- graph.adjacency(fungi.bpd.graphable, mode="undirected", weighted=NULL)
nodes_fungi_bpd <- gorder(fungi.bpd.subnet)
edges_fungi_bpd <- gsize(fungi.bpd.subnet)
degree_fungi_bpd <- as.data.frame(degree(fungi.bpd.subnet))

##There is certainly a more elegant and informative way to do this, and more I can do with this information, suggestions welcome##

###Get Subnetwork of all Fungi and all adjacent nodes for PPRD network###

fungi_pprd_subnetwork <- extract_surnodes(se.exp2.gen.pprd, c(61:80))
fungi.pprd.graphable <- as.matrix(fungi_pprd_subnetwork)

#Get metrics
colnames(fungi.pprd.graphable) <- row.names(exp2_prev_gen_otu[c(61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,10,23,47,51,4,20,29,35,36,44,9,38,40,5,11,15,33,3,6,18),])
fungi.pprd.subnet <- graph.adjacency(fungi.pprd.graphable, mode="undirected", weighted=NULL)
nodes_fungi_pprd <- gorder(fungi.pprd.subnet)
edges_fungi_pprd <- gsize(fungi.pprd.subnet)
degree_fungi_pprd <- as.data.frame(degree(fungi.pprd.subnet))


###EDGE PLOT WITH MAIN NETWORK AND FUNGAL SUBNETWORKS###
#Make data frame with edge count for fungi and connected bacteria
edges_fungi_df <- data.frame(edges_fungi_bpd,edges_fungi_pprd)
colnames(edges_fungi_df) <- c("BPD_Fungi","PPRD_Fungi")

#Combine with edge counts for overall plot
edges_df_all <- data.frame(edges_df,edges_fungi_df)
edges_long_all <- edges_df_all %>%
  pivot_longer(c("BPD_All","PPRD_All","BPD_Fungi","PPRD_Fungi"), names_to = "BPD", values_to = "Edges")

#Plot both
ggplot(edges_long_all,aes(x= BPD, y = Edges))+
  geom_bar(stat="identity", aes(fill=BPD))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save 

############################FIGURE E1E: Multikingdom SPIEC-EASI NODE PLOT############################

#Make data frame with node count for fungi and connected bacteria
nodes_fungi_df <- data.frame(nodes_fungi_bpd,nodes_fungi_pprd)
colnames(nodes_fungi_df) <- c("BPD_Fungi","PPRD_Fungi")

#Combine with node count from overall plot
nodes_df_all <- data.frame(nodes_df,nodes_fungi_df)
nodes_long_all <- nodes_df_all %>%
  pivot_longer(c("BPD_All","PPRD_All","BPD_Fungi","PPRD_Fungi"), names_to = "BPD", values_to = "Nodes")

#Plot both together
ggplot(nodes_long_all,aes(x= BPD, y = Nodes))+
  geom_bar(stat="identity", aes(fill=BPD))+
  scale_fill_manual(values = c("#007016","#007016","#94D57F","#94D57F")) +
  theme_bw()
#Save

############################FIGURE E1F: Multikingdom BPD SPIEC-EASI FULL NETWORK WITH TAXA LABELS############################

###Use SPIEC-EASI Network from Fig. 1H###

#BPD only network plot with taxa labels
prev_gen_spiec_plot_bpd_legend <- plot_network(spiec.graph.gen.bpd, exp2_combined_gen_prev_noNA, type='taxa',label=NULL, color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_bpd_legend + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save

############################FIGURE E1G: Multikingdom PPRD SPIEC-EASI FULL NETWORK WITH TAXA LABELS############################

###Use SPIEC-EASI Network from Fig. 1I###

#PPRD only network plot with taxa labels
prev_gen_spiec_plot_pprd <- plot_network(spiec.graph.gen.pprd, exp2_combined_gen_prev_noNA, type='taxa',label=NULL,color="Kingdom", point_size = 8, line_weight = 1)
prev_gen_spiec_plot_pprd + scale_color_manual(values = spiec.colors) + geom_text(mapping = aes(label = Genus), size = 3)
#Save