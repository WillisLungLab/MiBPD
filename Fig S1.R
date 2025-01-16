library(tidyverse)
library(PERFect)
library(phyloseq)
library(plyr)
library(RColorBrewer)
library(vegan)

### Fig. S1A: ITS Raw Rarefaction Curves ###

#Read in unfiltered ASV tables

tiff("Figure_S1A.tiff", width = 10, height = 5, units = "in", res = 600)
fung_curve <- rarecurve(fung_contam)
dev.off()
fung_contam_depth <- as.data.frame(colSums(fung_contam))
summary(fung_contam_depth)
#Save read count summary as .csv


### Fig. S1B: 16S Raw Rarefaction Curves ###

tiff("Figure_S1B.tiff", width = 10, height = 5, units = "in", res = 600)
bact_curve <- rarecurve(bact_contam)
write.csv(bact-contam_depth, "readcounts_16S_contam.csv")
dev.off()
bact_contam_depth <- as.data.frame(rowSums(bact_contam))
summary(bact_contam_depth)
#Save read count summary as .csv


### Fig. S1C: ITS Decontaminated Rarefaction Curves ###

#Read in filtered ASV tables (code for PERFect filtering in the Source Contamination folder)

tiff("Figure_S1C.tiff", width = 10, height = 5, units = "in", res = 600)
fung_curve <- rarecurve(fung_decontam)
dev.off()

fung_decontam_depth <- as.data.frame(rowSums(fung_decontam))
summary(fung_decontam_depth)
#Save read count summary as .csv


### Fig. S1D: 16S Decontaminated Rarefaction Curves ###

tiff("Figure_S1D.tiff", width = 10, height = 5, units = "in", res = 600)
bact_curve <- rarecurve(bact_decontam)
dev.off()

bact_decontam_depth <- as.data.frame(rowSums(bact_decontam))
summary(bact_decontam_depth)
#Save read count summary as .csv


### Fig. S1E: ITS Filtering Loss Plot ###

#Code for filtering with PERFect is in the Source Contamination folder

#Filtering loss plot
pvals_Plots(PERFect = res_perm2, X = fung_nozero, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
#Save

#Clean up and save results
fung_filt <- as.data.frame(t(res_perm2[["filtX"]]))
new_order = sort(row.names(fung_filt))
fung_filt <- as.data.frame(fung_filt[new_order,])
#Save filtered ASV table as .csv

#Get table of filtered taxa
unremoved <- colnames(fung_nozero[,colnames(fung_nozero) %in% row.names(fung_filt)])
#Read in unfiltered taxa
fung.taxa <- fung.taxa %>%
  tibble::column_to_rownames("Name")
unremoved.taxa <- fung.taxa[row.names(fung.taxa) %in% unremoved,]
#Save filtered taxonomy table as .csv


### Fig. S1F: 16S Filtering Loss Plot ###

#Code for filtering with PERFect is in the Source Contamination folder

pvals_Plots(PERFect = res_perm2, X = otu_for_perfect, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)

bact_filt2 <- as.data.frame(t(res_perm2[["filtX"]]))
new_order2 = sort(row.names(bact_filt2))
bact_filt2 <- as.data.frame(bact_filt2[new_order2,])
#Save filtered ASV table as .csv

#Get table of filtered taxa
unremoved <- row.names(bact_nozero[row.names(bact_nozero) %in% row.names(bact_filt),])
bact.taxa <- read_csv("Human_16s_Taxa.csv")
bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("Name")
unremoved.taxa <- bact.taxa[row.names(bact.taxa) %in% unremoved,]
#Save filtered taxonomy table as .csv

