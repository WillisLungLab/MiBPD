library(tidyverse)
library(phyloseq)
library(vegan)
library(viridis)
library(ecole)

############################FIGURE S7A: ITS PCOA BY OXYGEN CONCENTRATION############################

exp2_fung_rel_prev <- transform_sample_counts(exp2_fung_prev, function(x) x / sum(x) )

exp2_fung_otu_rel <- as.data.frame(t(exp2_fung_rel_prev@otu_table))
exp2_fung_tax_rel <- as.data.frame(exp2_fung_rel_prev@tax_table)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_rel_prev@sam_data)

exp2_fung_meta_rel2 <- exp2_fung_meta_rel %>% 
  mutate(Oxygen_Weeks_Binned = case_when(  exp2_fung_meta_rel$Oxygen_Weeks <= 5 ~ "0-5",
                                           exp2_fung_meta_rel$Oxygen_Weeks <= 10 ~ "6-10",
                                           exp2_fung_meta_rel$Oxygen_Weeks <= 15 ~ "11-15",
                                           exp2_fung_meta_rel$Oxygen_Weeks <= 20 ~ "16-20",
                                           exp2_fung_meta_rel$Oxygen_Weeks <= 26 ~ "21-26"
  )
  )

###Ordinate using Bray-Curtis dissimilarity###
exp2_fung_rel_bray = vegdist(exp2_fung_otu_rel, method='bray')
exp2_fung_rel_pcoa <- ape::pcoa(exp2_fung_rel_bray)
exp2_fung_rel_pcoa$values

###Set variable of interest and gradient color palette###
factor_exp2_ga <- ordered(exp2_fung_meta_rel2$Oxygen_Weeks_Binned, levels=c('0-5', '6-10', '11-15', '16-20','21-26'))
type_exp2_ga <- as.numeric(factor_exp2_ga)
pcoa_pal_ga <- c("#fde725","#5ec962","#21918c","#3b528b","#440154")


###Plot PCoA###
dpi=600
tiff("Fig E6A Binned.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.5), c(-0.5, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_fung_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pcoa_pal_ga[type_exp2_ga], lwd = 1)
dev.off()

###Create legend###
tiff("Legend_ITS.tif", width=5*dpi, height=5*dpi, res=dpi)
legend_image <- as.raster(matrix(rev(pcoa_pal_ga), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Weeks on Oxygen')
text(x=1.5, y = seq(0,1,l=5), labels = seq(1))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

###PERMANOVA, omit any NA value###
permanova_pcoa_df <- data.frame(exp2_fung_meta_rel2)
set.seed(1312)
permanova_fung_oxygen_weeks <- vegan::adonis2(exp2_fung_otu_rel ~ Oxygen_Weeks_Binned, data = permanova_pcoa_df, method="bray", na.action = na.omit,permutations = 10000)
print(permanova_fung_oxygen_weeks)

###PERMDISP###
fung_pcoa_dist <- vegdist(exp2_fung_otu_rel_noNA, method = "bray")
disp_fung_pcoa <- betadisper(fung_pcoa_dist, permanova_pcoa_df_noNA$Oxygen_Weeks_Binned)
set.seed(1312)
permdisp_fung_pcoa <- permutest(disp_fung_pcoa, permutations = 10000)

print(permanova_fung_oxygen_weeks)
print(permdisp_fung_pcoa)


############################FIGURE S7B: 16S PCOA BY OXYGEN CONCENTRATION############################

exp2_bact_rel_prev <- transform_sample_counts(exp2_bact_prev, function(x) x / sum(x) )

exp2_bact_otu_rel <- as.data.frame(t(exp2_bact_rel_prev@otu_table))
exp2_bact_tax_rel <- as.data.frame(exp2_bact_rel_prev@tax_table)
exp2_bact_meta_rel <- as.data.frame(exp2_bact_rel_prev@sam_data)

exp2_bact_meta_rel2 <- exp2_bact_meta_rel %>% 
  mutate(Oxygen_Weeks_Binned = case_when(  exp2_bact_meta_rel$Oxygen_Weeks <= 5 ~ "0-5",
                                           exp2_bact_meta_rel$Oxygen_Weeks <= 10 ~ "6-10",
                                           exp2_bact_meta_rel$Oxygen_Weeks <= 15 ~ "11-15",
                                           exp2_bact_meta_rel$Oxygen_Weeks <= 20 ~ "16-20",
                                           exp2_bact_meta_rel$Oxygen_Weeks <= 26 ~ "21-26"
  )
  )

###Ordinate using Bray-Curtis dissimilarity###
exp2_bact_rel_bray = vegdist(exp2_bact_otu_rel, method='bray')
exp2_bact_rel_pcoa <- ape::pcoa(exp2_bact_rel_bray)
exp2_bact_rel_pcoa$values

###Set variable of interest and gradient color palette###
factor_exp2_ga <- ordered(exp2_bact_meta_rel2$Oxygen_Weeks_Binned, levels=c('0-5', '6-10', '11-15', '16-20','21-26'))
type_exp2_ga <- as.numeric(factor_exp2_ga)
pcoa_pal_ga <- c("#fde725","#5ec962","#21918c","#3b528b","#440154")

###Plot PCoA, using same legend as ITS so no need to plot that twice###
dpi=600
tiff("Fig E6B Binned.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.4, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp2_bact_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pcoa_pal_ga[type_exp2_ga], lwd = 1)
dev.off()

###PERMANOVA###
permanova_pcoa_df <- data.frame(exp2_bact_meta_rel2)
set.seed(1312)
permanova_bact_oxygen_weeks <- vegan::adonis2(exp2_bact_otu_rel ~ Oxygen_Weeks_Binned, data = permanova_pcoa_df, method="bray", na.action = na.omit,permutations = 10000)
print(permanova_bact_oxygen_weeks)
oxygen_weeks_p <- permanova_bact_oxygen_weeks$`Pr(>F)`

###Pairwise PERMANOVA###
permanova_pcoa_df_noNA <-permanova_pcoa_df[-c(2,4,62,65,68,93),]
exp2_bact_otu_rel_noNA <-exp2_bact_otu_rel[-c(2,4,62,65,68,93),]

pairwise_permanova_bact_oxygen_weeks <- permanova_pairwise(exp2_bact_otu_rel_noNA, grp = permanova_pcoa_df_noNA$Oxygen_Weeks_Binned, permutations = 10000, method = "bray", padj = "fdr")

###PERMDISP###
bact_pcoa_dist <- vegdist(exp2_bact_otu_rel, method = "bray")
disp_bact_pcoa <- betadisper(bact_pcoa_dist, permanova_pcoa_df$Oxygen_Weeks_Binned)
set.seed(1312)
permdisp_bact_pcoa <- permutest(disp_bact_pcoa, permutations = 10000)

print(permanova_bact_oxygen_weeks)
print(permdisp_bact_pcoa)


###Legend###

dpi=600
tiff("Legend_16S.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('0-5', '6-10', '11-15',
                            '16-20', '21-26'), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#fde725","#5ec962","#21918c","#3b528b","#440154"))
mtext("Weeks on Oxygen", at=0.1, cex=1.5)
dev.off()
