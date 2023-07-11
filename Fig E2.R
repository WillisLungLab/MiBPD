library(tidyverse)
library(readxl)
library(vegan)
library(phyloseq)
library(DESeq2)
library(indicspecies)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(factoextra)

###############DATA PREP###############

#Load ASV table
bact <- bact %>%
  tibble::column_to_rownames("Name")

#Load taxonomy table
#Convert sample names to row names and convert to a matrix to make converting to phyloseq easier
bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("Name")
bact.taxa <- as.matrix(bact.taxa)

#Load metadata and the rest of this is all the same
bact.sort <- bact[,order(colnames(bact))]
md.bact <- roughmetadata_bact[which(unlist(roughmetadata_bact$Name) %in% colnames(bact.sort)),]
summary(colSums(bact.sort[,-1]))

md.bact <- md.bact %>%
  tibble::column_to_rownames("Name")

#Filter out samples with read depth <250, 2 removed
bact.sort2 <- bact.sort[,-c(1,which(colSums(bact.sort[,-1])<250))]
md.bact.sort2 <- md.bact[which(row.names(md.bact) %in% colnames(bact.sort2)),]

#This is the same
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

bact.nsampls <- apply(bact.sort2, 1, gt0)
bact.clean <- bact.sort2[which(bact.nsampls>1),]

#Setup to create the phyloseq object
OTU = otu_table(bact.clean, taxa_are_rows = TRUE)
TAX = tax_table(bact.taxa)
samples = sample_data(md.bact.sort2)

#Convert!
exp2_bact <- phyloseq(OTU, TAX, samples)
exp2_bact_prev <- filter_taxa(exp2_bact , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)

###############Fig. E2A###############

#Merge samples by BPD status
exp2_bact_merged = merge_samples(exp2_bact_prev, "BPD")
exp2_bact_gen_merged <- tax_glom(exp2_bact_merged, taxrank = 'Genus')

#Gets a list of the top 19 bacterial genera and labels everything else as "Other"
top20_bact_merged_list <- names(sort(taxa_sums(exp2_bact_gen_merged), decreasing=TRUE)[1:19])
top20_merged_bact_rel <- transform_sample_counts(exp2_bact_gen_merged, function(x) x / sum(x) )
top20_merged_bact_df <- psmelt(top20_merged_bact_rel)
top20_merged_bact_df[!(top20_merged_bact_df$OTU %in% top20_bact_merged_list),]$Genus <- 'Other'

#Plot barplot
barplot_colors <- colorRampPalette(brewer.pal(12, "Paired"))(20)

barplot_gen_bpd_its <- ggplot(top20_merged_bact_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90)) 

barplot_gen_bpd_its + scale_fill_manual(values = barplot_colors)
#Save

###############Fig. E2B###############

OTU_rough = otu_table(bact.sort, taxa_are_rows = TRUE)
TAX_rough = tax_table(bact.taxa)
samples_rough = sample_data(md.bact)

exp2_bact_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

fr_bact <- prune_taxa(taxa_sums(exp2_bact_rough) > 0, exp2_bact_rough)
richness_est_bact <- estimate_richness(fr_bact, measures = c("Chao1", "Shannon"))

#Wilcox test for significance
wilcox_alpha_bact <- t(sapply(richness_est_bact, function(x) unlist(wilcox.test(x~sample_data(fr_bact)$BPD)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_bact

#Get dataframe with alpha diversity values and metadata values, to plot in GraphPad
richness_est_bact <- richness_est_bact %>%
  mutate(
    BPD = md.bact$BPD,
    Sex = md.bact$Sex,
    BPD_Sex = md.bact$BPD_Sex,
    Oxygen = md.bact$Oxygen_Days
  )
#Save as .csv

###############Fig. E2C###############

exp2_bact_rel_prev <- transform_sample_counts(exp2_bact_prev, function(x) x / sum(x) )

exp2_bact_otu_rel <- as.data.frame(t(exp2_bact_rel_prev@otu_table))
exp2_bact_tax_rel <- as.data.frame(exp2_bact_rel_prev@tax_table)
exp2_bact_meta_rel <- as.data.frame(exp2_bact_rel_prev@sam_data)

#PCoA using Bray-Curtis dissimilarity
exp2_bact_rel_bray = vegdist(exp2_bact_otu_rel, method='bray')
exp2_bact_rel_pcoa <- ape::pcoa(exp2_bact_rel_bray)
exp2_bact_rel_pcoa$values

factor_exp2 <- as.factor(exp2_bact_meta_rel$BPD)
type_exp2 <- as.numeric(factor_exp2)
pca_colors <- c("#007016","#94D57F")

#Spider plot with confidence interval
dpi=600
tiff("Beta Diversity PCoA 16s.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.4), c(-0.4, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp2], lwd = 1)
ordiellipse(exp2_bact_rel_pcoa$vectors[,1:2], kind = "se", conf = 0.95, factor_exp2, col = pca_colors)
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


