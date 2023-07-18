############################FIGURE E8A: MULTIKINGDOM DMM PCoA############################

exp2_combined_tse <- makeTreeSEFromPhyloseq(exp2_combined_prev)

set.seed(1312)
tse_dmn <- mia::runDMN(exp2_combined_tse, name = "DMN", k = 1:7)
best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")

prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
# Add column names
colnames(prob) <- c("comp1", "comp2","comp3", "comp4")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec <- colnames(prob)[max.col(prob,ties.method = "first")]

exp2_combined_rel_prev <- transform_sample_counts(exp2_combined_prev, function(x) x / sum(x) )

exp2_combined_otu_rel <- as.data.frame(t(exp2_combined_rel_prev@otu_table))
exp2_combined_tax_rel <- as.data.frame(exp2_combined_rel_prev@tax_table)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_rel_prev@sam_data)
exp2_combined_meta_rel <- cbind(exp2_combined_meta_rel,
                                dmm_component = vec)

###Ordinate using Bray-Curtis distance###
exp2_combined_rel_bray = vegdist(exp2_combined_otu_rel, method='bray')
exp2_combined_rel_pcoa <- ape::pcoa(exp2_combined_rel_bray)
bray_pcoa_df <- exp2_combined_rel_pcoa$vectors[,1:2]

bray_dmm_pcoa_df <- cbind(bray_pcoa_df,
                          dmm_component = vec)

factor_exp2 <- as.factor(bray_dmm_pcoa_df[,3])
type_exp2 <- as.numeric(factor_exp2)
dmm_colors_all <- c("#8FD96C","#64BD4B","#228019","#096620")

###Plot###
dpi=600
tiff("Multikingdom All DMM PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.4, 0.6), font = 2, font.lab = 2, xlab="PC1(5.41% Explained)", ylab="PC2(4.35% Explained)", type="n")
points(exp2_combined_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_all[type_exp2], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2)
ordispider(exp2_combined_rel_pcoa$vectors[,1:2], factor_exp2, label = TRUE)
dev.off()


############################FIGURE E8B: MULTIKINGDOM DMM ALPHA DIVERSITY############################

###Get Chao1 and Shannon diversity###
fr_combined <- prune_taxa(taxa_sums(exp2_combined_rough) > 0, exp2_combined_rough)
richness_est_combined <- estimate_richness(fr_combined, measures = c("Chao1", "Shannon"))

###Add cluster assignments to dataframe with alpha diversity values###
vec_df <- as.data.frame(vec)
row.names(vec_df) <- row.names(exp2_combined_meta_rel)
richness_est_combined_cluster <- richness_est_combined[which(unlist(row.names(richness_est_combined)) %in% row.names(vec_df)),]
richness_est_combined_cluster <- cbind(richness_est_combined_cluster, vec_df)

#Save and create plot in GraphPad


############################FIGURE E8C: MULTIKINGDOM DMM TAXA BARPLOT############################

###Merge samples by cluster###
exp2_combined_merged = merge_samples(exp2_combined_prev, "cluster")
exp2_combined_gen_merged <- tax_glom(exp2_combined_merged, taxrank = 'Genus')

###Get top 19 genera and an Other entry for all others###
top20_combined_prev_list <- names(sort(taxa_sums(exp2_combined_gen_merged), decreasing=TRUE)[1:19])
top20_prev_combined_rel <- transform_sample_counts(exp2_combined_gen_merged, function(x) x / sum(x) )
top20_prev_combined_df <- psmelt(top20_prev_combined_rel)
top20_prev_combined_df[!(top20_prev_combined_df$OTU %in% top20_combined_prev_list),]$Genus <- 'Other'


###Plot barplot###
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(20)

barplot_gen_bpd_combined <- ggplot(top20_prev_combined_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_combined + scale_fill_manual(values = barplot_colors)
#Save


############################FIGURE E8D: MULTIKINGDOM PPRD DMM PCoA############################

exp2_combined_tse_pprd <- makeTreeSEFromPhyloseq(exp2_combined_prev_pprd)

set.seed(1312)
tse_dmn_pprd <- mia::runDMN(exp2_combined_tse_pprd, name = "DMN", k = 1:7)
best_dmn_pprd <- getBestDMNFit(tse_dmn_pprd, type = "laplace")

prob_pprd <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_pprd))
# Add column names
colnames(prob_pprd) <- c("comp1", "comp2","comp3")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec_pprd <- colnames(prob_pprd)[max.col(prob_pprd,ties.method = "first")]

exp2_combined_rel_prev_pprd <- transform_sample_counts(exp2_combined_prev_pprd, function(x) x / sum(x) )

exp2_combined_otu_rel_pprd <- as.data.frame(t(exp2_combined_rel_prev_pprd@otu_table))
exp2_combined_tax_rel_pprd <- as.data.frame(exp2_combined_rel_prev_pprd@tax_table)
exp2_combined_meta_rel_pprd <- as.data.frame(exp2_combined_rel_prev_pprd@sam_data)
exp2_combined_meta_rel_pprd <- cbind(exp2_combined_meta_rel_pprd,
                                     dmm_component = vec_pprd)

exp2_combined_rel_bray_pprd = vegdist(exp2_combined_otu_rel_pprd, method='bray')
exp2_combined_rel_pcoa_pprd <- ape::pcoa(exp2_combined_rel_bray_pprd)
bray_pcoa_df_pprd <- exp2_combined_rel_pcoa_pprd$vectors[,1:2]

bray_dmm_pcoa_df_pprd <- cbind(bray_pcoa_df_pprd,
                               dmm_component = vec_pprd)

factor_exp2_pprd <- as.factor(bray_dmm_pcoa_df_pprd[,3])
type_exp2_pprd <- as.numeric(factor_exp2_pprd)
dmm_colors_pprd <- c("#8FD96C","#409E2F","#096620")

###Plot PCoA###
dpi=600
tiff("Combined DMM Clusters PCoA PPRD Only.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.6), c(-0.5, 0.5), font = 2, font.lab = 2, xlab="PC1(8.01% Explained)", ylab="PC2(6.29% Explained)", type="n")
points(exp2_combined_rel_pcoa_pprd$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_pprd[type_exp2_pprd], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa_pprd$vectors[,1:2], factor_exp2_pprd)
ordispider(exp2_combined_rel_pcoa_pprd$vectors[,1:2], factor_exp2_pprd, label = TRUE)
dev.off()


############################FIGURE E8E: MULTIKINGDOM PPRD DMM ALPHA DIVERSITY############################

vec_pprd_df <- as.data.frame(vec_pprd)
row.names(vec_pprd_df) <- row.names(exp2_combined_meta_rel_pprd)
richness_est_combined_cluster_pprd <- richness_est_combined[which(unlist(row.names(richness_est_combined)) %in% row.names(vec_pprd_df)),]

richness_est_combined_cluster_pprd <- cbind(richness_est_combined_cluster_pprd, vec_pprd_df)
#Save and plot in GraphPad


############################FIGURE E8F: MULTIKINGDOM PPRD DMM TAXA BARPLOT############################

exp2_combined_merged_pprd = merge_samples(exp2_combined_prev_pprd, "cluster")
exp2_combined_gen_merged_pprd <- tax_glom(exp2_combined_merged_pprd, taxrank = 'Genus')

###Get top 19 genera and create an "Other" designation for all other genera###
top20_combined_prev_list_pprd <- names(sort(taxa_sums(exp2_combined_gen_merged_pprd), decreasing=TRUE)[1:19])
top20_prev_combined_rel_pprd <- transform_sample_counts(exp2_combined_gen_merged_pprd, function(x) x / sum(x) )
top20_prev_combined_df_pprd <- psmelt(top20_prev_combined_rel_pprd)
top20_prev_combined_df_pprd[!(top20_prev_combined_df_pprd$OTU %in% top20_combined_prev_list_pprd),]$Genus <- 'Other'

###Plot###
barplot_gen_pprd_combined_pprd <- ggplot(top20_prev_combined_df_pprd, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_pprd_combined_pprd + scale_fill_manual(values = barplot_colors)
#Save


############################FIGURE E8G: MULTIKINGDOM BPD DMM PCoA############################

exp2_combined_tse_bpd <- makeTreeSEFromPhyloseq(exp2_combined_prev_bpd)

set.seed(1312)
tse_dmn_bpd <- mia::runDMN(exp2_combined_tse_bpd, name = "DMN", k = 1:7)
best_dmn_bpd <- getBestDMNFit(tse_dmn_bpd, type = "laplace")

prob_bpd <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_bpd))
# Add column names
colnames(prob_bpd) <- c("comp1", "comp2")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec_bpd <- colnames(prob_bpd)[max.col(prob_bpd,ties.method = "first")]

exp2_combined_rel_prev_bpd <- transform_sample_counts(exp2_combined_prev_bpd, function(x) x / sum(x) )

exp2_combined_otu_rel_bpd <- as.data.frame(t(exp2_combined_rel_prev_bpd@otu_table))
exp2_combined_tax_rel_bpd <- as.data.frame(exp2_combined_rel_prev_bpd@tax_table)
exp2_combined_meta_rel_bpd <- as.data.frame(exp2_combined_rel_prev_bpd@sam_data)
exp2_combined_meta_rel_bpd <- cbind(exp2_combined_meta_rel_bpd,
                                    dmm_component = vec_bpd)

###Ordinate with Bray-Curtis dissimilarity###
exp2_combined_rel_bray_bpd = vegdist(exp2_combined_otu_rel_bpd, method='bray')
exp2_combined_rel_pcoa_bpd <- ape::pcoa(exp2_combined_rel_bray_bpd)
bray_pcoa_df_bpd <- exp2_combined_rel_pcoa_bpd$vectors[,1:2]

bray_dmm_pcoa_df_bpd <- cbind(bray_pcoa_df_bpd,
                              dmm_component = vec_bpd)

factor_exp2_bpd <- as.factor(bray_dmm_pcoa_df_bpd[,3])
type_exp2_bpd <- as.numeric(factor_exp2_bpd)
dmm_colors_all_bpd <- c("#8FD96C","#228019")

###Plot###
dpi=600
tiff("Combined DMM PCoA BPD Only.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.5), c(-0.6, 0.6), font = 2, font.lab = 2, xlab="PC1(11.85% Explained)", ylab="PC2(8.99% Explained)", type="n")
points(exp2_combined_rel_pcoa_bpd$vectors[,1:2], pch = 21, cex = 1.3, bg = dmm_colors_all_bpd[type_exp2_bpd], lwd = 1)
ordiellipse(exp2_combined_rel_pcoa_bpd$vectors[,1:2], factor_exp2_bpd)
ordispider(exp2_combined_rel_pcoa_bpd$vectors[,1:2], factor_exp2_bpd, label = TRUE)
dev.off()


############################FIGURE E8H: MULTIKINGDOM BPD DMM ALPHA DIVERSITY############################

vec_bpd_df <- as.data.frame(vec_bpd)
row.names(vec_bpd_df) <- row.names(exp2_combined_meta_rel_bpd)
richness_est_combined_cluster_bpd <- richness_est_combined[which(unlist(row.names(richness_est_combined)) %in% row.names(vec_bpd_df)),]

richness_est_combined_cluster_bpd <- cbind(richness_est_combined_cluster_bpd, vec_bpd_df)
#Save


############################FIGURE E8I: MULTIKINGDOM BPD DMM TAXA BARPLOT############################

exp2_combined_merged_bpd = merge_samples(exp2_combined_prev_bpd, "cluster")
exp2_combined_gen_merged_bpd <- tax_glom(exp2_combined_merged_bpd, taxrank = 'Genus')

###Get top 20 genera, including "Other"###
top20_combined_prev_list_bpd <- names(sort(taxa_sums(exp2_combined_gen_merged_bpd), decreasing=TRUE)[1:19])
top20_prev_combined_rel_bpd <- transform_sample_counts(exp2_combined_gen_merged_bpd, function(x) x / sum(x) )
top20_prev_combined_df_bpd <- psmelt(top20_prev_combined_rel_bpd)
top20_prev_combined_df_bpd[!(top20_prev_combined_df_bpd$OTU %in% top20_combined_prev_list_bpd),]$Genus <- 'Other'

###Plot barplot###
barplot_gen_bpd_combined_bpd <- ggplot(top20_prev_combined_df_bpd, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_combined_bpd + scale_fill_manual(values = barplot_colors)
#Save