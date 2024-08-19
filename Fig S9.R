library(mia)
library(miaViz)
library(phyloseq)
library(tidyverse)
library(DirichletMultinomial)

############################FIGURE S9A: 16S DMM LAPLACE PLOT############################

exp2_bact_tse <- makeTreeSEFromPhyloseq(exp2_bact_prev)

set.seed(1312)
tse_dmn <- mia::runDMN(exp2_bact_tse, name = "DMN", k = 1:7)
best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")

plotDMNFit(tse_dmn, type = "laplace")
#Save


############################FIGURE S9B: 16S DMM VARIABLE IMPORTANCE PLOT############################

best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
fitted <- head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))

exp2_bact_tax_rel$Genus <- make.unique(exp2_bact_tax_rel$Genus)

#For each ASV, assign value per cluster. This indicates importance of ASV to cluster composition
for (k in seq(ncol(fitted(best_dmn)))) {
  d <- melt(fitted(best_dmn))
  colnames(d) <- c("OTU", "cluster", "value")
  print(d)
}

exp2_bact_tax_rel_genonly <- data.frame(new_col = c(exp2_bact_tax_rel$Genus, exp2_bact_tax_rel$Genus, exp2_bact_tax_rel$Genus,exp2_bact_tax_rel$Genus))
d <- d %>%
  mutate(Genus = exp2_bact_tax_rel_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  #Only show ASVs with value over 0.7
  filter(abs(value) > quantile(abs(value), 0.70))
print(d)

###Cluster 1###
d1 <- subset(d, d$cluster == 1)
d1 <- d1 %>% top_n(10, value)
p1 <- ggplot(d1, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 1")) +
  theme_bw()
print(p1)
#Save

###Cluster 2###
d2 <- subset(d, d$cluster == 2)
d2 <- d2 %>% top_n(10, value)
p2 <- ggplot(d2, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 2")) +
  theme_bw()
print(p2)
#Save

###Cluster 3###
d3 <- subset(d, d$cluster == 3)
d3 <- d3 %>% top_n(10, value)
p3 <- ggplot(d3, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 3")) +
  theme_bw()
print(p3)
#Save

###Cluster 4###
d4 <- subset(d, d$cluster == 4)
d4 <- d4 %>% top_n(10, value)
p4 <- ggplot(d4, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 4")) +
  theme_bw()
print(p4)
#Save


############################FIGURE S9C: 16S DMM CCAS############################

#This metadata file has not been made available since it contains clinical data
md_bact_newfactors_sort$cluster <- vec
exp2_bact_meta_df <- as.data.frame(md_bact_newfactors_sort[,c(1,3,4,6,7,8,25,28,29,30)])
exp2_bact_meta_df[sapply(exp2_bact_meta_df, is.character)] <- lapply(exp2_bact_meta_df[sapply(exp2_bact_meta_df, is.character)], as.factor)
sapply(exp2_bact_meta_df, base::class)

exp2_bact_meta_df$BPD<- factor(exp2_bact_meta_df$BPD, levels=c('No_BPD', 'BPD'))
exp2_bact_meta_df$Died<- factor(exp2_bact_meta_df$Died, levels=c('Survived', 'Died'))
exp2_bact_meta_df$cluster<- factor(exp2_bact_meta_df$cluster, levels=c('Cluster 4','Cluster 2','Cluster 3','Cluster 1'))
exp2_bact_meta_df_noNA <- exp2_bact_meta_df[complete.cases(exp2_bact_meta_df), ]
exp2_bact_otu_rel_noNA <- exp2_bact_otu_rel[which(row.names(exp2_bact_otu_rel) %in% row.names(exp2_bact_meta_df_noNA)),]

cca_result2 <- vegan::cca(exp2_bact_otu_rel_noNA,exp2_bact_meta_df_noNA)


#summary of the data
summary(cca_result2)
plot(cca_result2)

#extracting data: clinical metadata
veg_1 = as.data.frame(cca_result2$CCA$biplot)
row.names(veg_1) <- c("BPD","Died","AMAB","GA","BW","CRIB.II","CRP","No_Fungi","Low_Evenness","Fungi_Prev","Cluster 2","Cluster 3","Cluster 1")
veg_1["env"] = row.names(veg_1)

#extracting data: samples
veg_2 = as.data.frame(cca_result2$CCA$u)
veg_2["samples"] = row.names(veg_2)

#extracting data: taxa
veg_3 = as.data.frame(cca_result2$CCA$v)
veg_3["genus"] = row.names(veg_3)


#Plot points
plot = ggplot() +
  geom_point(data = veg_3, aes(x = CCA1, y = CCA2), color = "#AFEE91", size = 1.0) +
  geom_point(data = veg_2, aes(x = CCA1, y = CCA2), color = "#0B8700", size = 2.0) +
  geom_point(data = veg_1, aes(x = CCA1, y = CCA2), color = "black")


plot

#Plot arrows and labels
plot +
  theme_bw() +
  geom_segment(
    data = veg_1,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.5, "cm"))
  ) +
  geom_text_repel(
    data = veg_1,
    aes(x = CCA1, y = CCA2, label = veg_1$env),
    nudge_y = -0.05,
    color = "black",
    size = 3
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
#Save


############################FIGURE S9D: 16S NOBPD DMM LAPLACE PLOT############################

exp2_bact_tse_pprd <- makeTreeSEFromPhyloseq(exp2_bact_prev_pprd)

set.seed(1312)
tse_dmn_pprd <- mia::runDMN(exp2_bact_tse_pprd, name = "DMN", k = 1:7)
getDMN(tse_dmn_pprd)

plotDMNFit(tse_dmn_pprd, type = "laplace")
#Save


############################FIGURE S9E: 16S NOBPD DMM VARIABLE IMPORTANCE PLOT############################

#Get DMM for k = 2 as in Fig S8
set.seed(1312)
tse_dmn_pprd <- mia::runDMN(exp2_bact_tse_pprd, name = "DMN", k = 2)
best_dmn_pprd <- getBestDMNFit(tse_dmn_pprd, type = "laplace")

exp2_bact_tax_rel_pprd$Genus <- make.unique(exp2_bact_tax_rel_pprd$Genus)

for (k in seq(ncol(fitted(best_dmn_pprd)))) {
  d_pprd <- melt(fitted(best_dmn_pprd))
  colnames(d_pprd) <- c("OTU", "cluster", "value")
  print(d_pprd)
}

exp2_bact_tax_rel_pprd_genonly <- data.frame(new_col = c(exp2_bact_tax_rel_pprd$Genus, exp2_bact_tax_rel_pprd$Genus))
d_pprd <- d_pprd %>%
  mutate(Genus = exp2_bact_tax_rel_pprd_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.70))
print(d_pprd)

###Cluster 1###
d1_pprd <- subset(d_pprd, d_pprd$cluster == 1)
d1_pprd <- d1_pprd %>% top_n(10, value)
p1_pprd <- ggplot(d1_pprd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 1")) +
  theme_bw()
print(p1_pprd)
#Save

###Cluster 2###
d2_pprd <- subset(d_pprd, d_pprd$cluster == 2)
d2_pprd <- d2_pprd %>% top_n(10, value)
p2_pprd <- ggplot(d2_pprd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 2")) +
  theme_bw()
print(p2_pprd)
#Save

############################FIGURE S9F: 16S NOBPD DMM CCAS############################

###Get metadata###
md_bact_newfactors_sort_pprd <- md_bact_newfactors_sort[md_bact_newfactors_sort$BPD == "No_BPD",]
md_bact_newfactors_sort_pprd$cluster <- vec_pprd
exp2_bact_meta_df_pprd <- as.data.frame(md_bact_newfactors_sort_pprd[,c(1,3,4,6,7,8,25,28,29,30)])
exp2_bact_meta_df_pprd[sapply(exp2_bact_meta_df_pprd, is.character)] <- lapply(exp2_bact_meta_df_pprd[sapply(exp2_bact_meta_df_pprd, is.character)], as.factor)
sapply(exp2_bact_meta_df, base::class)

exp2_bact_meta_df_pprd$BPD<- factor(exp2_bact_meta_df_pprd$BPD, levels=c('No_BPD', 'BPD'))
exp2_bact_meta_df_pprd$Died<- factor(exp2_bact_meta_df_pprd$Died, levels=c('Survived', 'Died'))
exp2_bact_meta_df_pprd$cluster<- factor(exp2_bact_meta_df_pprd$cluster, levels=c('Cluster 2', 'Cluster 1'))
exp2_bact_meta_df_pprd$Dominance<- factor(exp2_bact_meta_df_pprd$Dominance, levels=c('No_Fungi', 'FALSE','TRUE'))

exp2_bact_meta_df_pprd_noNA <- exp2_bact_meta_df_pprd[complete.cases(exp2_bact_meta_df_pprd), ]
exp2_bact_otu_rel_pprd_noNA <- exp2_bact_otu_rel_pprd[which(row.names(exp2_bact_otu_rel_pprd) %in% row.names(exp2_bact_meta_df_pprd_noNA)),]

cca_result_pprd <- vegan::cca(exp2_bact_otu_rel_pprd_noNA, exp2_bact_meta_df_pprd_noNA)

#summary
summary(cca_result_pprd)
plot(cca_result_pprd)

###Extract metadata###
veg_1_pprd = as.data.frame(cca_result_pprd$CCA$biplot)
row.names(veg_1_pprd) <- c("AMAB","GA","BW","CRIB.II","CRP","High_Evenness","Low_Evenness","Fungi_Prev","Cluster 1")
veg_1_pprd["env"] = row.names(veg_1_pprd)

###Extract sample data###
veg_2_pprd = as.data.frame(cca_result_pprd$CCA$u)
veg_2_pprd["samples"] = row.names(veg_2_pprd)

###Extract taxa###
veg_3_pprd = as.data.frame(cca_result_pprd$CCA$v)
veg_3_pprd["genus"] = row.names(veg_3_pprd)

###Plot points###
plot_pprd = ggplot() +
  geom_point(data = veg_3_pprd, aes(x = CCA1, y = CCA2), color = "#AFEE91", size = 1.0) +
  geom_point(data = veg_2_pprd, aes(x = CCA1, y = CCA2), color = "#0B8700", size = 2.0) +
  geom_point(data = veg_1_pprd, aes(x = CCA1, y = CCA2), color = "black")


plot_pprd

###Plot arrows and labels###
plot_pprd +
  theme_bw() +
  geom_segment(
    data = veg_1_pprd,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.5, "cm"))
  ) +
  geom_text_repel(
    data = veg_1_pprd,
    aes(x = CCA1, y = CCA2, label = veg_1_pprd$env),
    nudge_y = -0.05,
    color = "black",
    size = 3
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
#Save


############################FIGURE S9G: 16S BPD DMM LAPLACE PLOT############################

exp2_bact_tse_bpd <- makeTreeSEFromPhyloseq(exp2_bact_prev_bpd)

set.seed(1312)
tse_dmn_bpd <- mia::runDMN(exp2_bact_tse_bpd, name = "DMN", k = 1:7)
getDMN(tse_dmn_bpd)

plotDMNFit(tse_dmn_bpd, type = "laplace")
#Save


############################FIGURE S9H: 16S BPD DMM VARIABLE IMPORTANCE PLOT############################

best_dmn_bpd <- getBestDMNFit(tse_dmn_bpd, type = "laplace")

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn_bpd))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_bpd)))
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn_bpd)))

exp2_bact_tax_rel_bpd$Genus <- make.unique(exp2_bact_tax_rel_bpd$Genus)

for (k in seq(ncol(fitted(best_dmn_bpd)))) {
  d_bpd <- melt(fitted(best_dmn_bpd))
  colnames(d_bpd) <- c("OTU", "cluster", "value")
  print(d_bpd)
}

exp2_bact_tax_rel_bpd_genonly <- data.frame(new_col = c(exp2_bact_tax_rel_bpd$Genus, exp2_bact_tax_rel_bpd$Genus))
d_bpd <- d_bpd %>%
  mutate(Genus = exp2_bact_tax_rel_bpd_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.70))
print(d_bpd)

###Cluster 1###
d1_bpd <- subset(d_bpd, d_bpd$cluster == 1)
d1_bpd <- d1_bpd %>% top_n(10, value)
p1_bpd <- ggplot(d1_bpd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 1")) +
  theme_bw()
print(p1_bpd)
#Save

###Cluster 2###
d2_bpd <- subset(d_bpd, d_bpd$cluster == 2)
d2_bpd <- d2_bpd %>% top_n(10, value)
p2_bpd <- ggplot(d2_bpd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: Cluster 2")) +
  theme_bw()
print(p2_bpd)
Save


############################FIGURE S9I: 16S BPD DMM CCAS############################

md_bact_newfactors_sort_bpd <- md_bact_newfactors_sort[md_bact_newfactors_sort$BPD == "BPD",]
md_bact_newfactors_sort_bpd$cluster <- vec_bpd
exp2_bact_meta_df_bpd <- as.data.frame(md_bact_newfactors_sort_bpd[,c(1,3,4,6,7,8,25,28,29,30)])
exp2_bact_meta_df_bpd[sapply(exp2_bact_meta_df_bpd, is.character)] <- lapply(exp2_bact_meta_df_bpd[sapply(exp2_bact_meta_df_bpd, is.character)], as.factor)
sapply(exp2_bact_meta_df, base::class)

exp2_bact_meta_df_bpd$BPD<- factor(exp2_bact_meta_df_bpd$BPD, levels=c('No_BPD', 'BPD'))
exp2_bact_meta_df_bpd$Died<- factor(exp2_bact_meta_df_bpd$Died, levels=c('Survived', 'Died'))
exp2_bact_meta_df_bpd$cluster<- factor(exp2_bact_meta_df_bpd$cluster, levels=c('Cluster 1', 'Cluster 2'))
exp2_bact_meta_df_bpd$Dominance<- factor(exp2_bact_meta_df_bpd$Dominance, levels=c('No_Fungi', 'FALSE','TRUE'))

exp2_bact_meta_df_bpd_noNA <- exp2_bact_meta_df_bpd[complete.cases(exp2_bact_meta_df_bpd), ]
exp2_bact_otu_rel_bpd_noNA <- exp2_bact_otu_rel_bpd[which(row.names(exp2_bact_otu_rel_bpd) %in% row.names(exp2_bact_meta_df_bpd_noNA)),]

cca_result_bpd <- vegan::cca(exp2_bact_otu_rel_bpd_noNA, exp2_bact_meta_df_bpd_noNA)

summary(cca_result_bpd)
plot(cca_result_bpd)

###Extract metadata###
veg_1_bpd = as.data.frame(cca_result_bpd$CCA$biplot)
row.names(veg_1_bpd) <- c("Died","AMAB","GA","BW","CRIB.II","CRP","High_Evenness","Low_Evenness","Fungi_Prev","Cluster 2")
veg_1_bpd["env"] = row.names(veg_1_bpd)

###Extract sample data###
veg_2_bpd = as.data.frame(cca_result_bpd$CCA$u)
veg_2_bpd["samples"] = row.names(veg_2_bpd)

###Extract taxa###
veg_3_bpd = as.data.frame(cca_result_bpd$CCA$v)
veg_3_bpd["genus"] = row.names(veg_3_bpd)

###Plot points###
plot_bpd = ggplot() +
  geom_point(data = veg_3_bpd, aes(x = CCA1, y = CCA2), color = "#AFEE91", size = 1.0) +
  geom_point(data = veg_2_bpd, aes(x = CCA1, y = CCA2), color = "#0B8700", size = 2.0) +
  geom_point(data = veg_1_bpd, aes(x = CCA1, y = CCA2), color = "black")


plot_bpd

###Plot arrows and labels###
plot_bpd +
  theme_bw() +
  geom_segment(
    data = veg_1_bpd,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.5, "cm"))
  ) +
  geom_text_repel(
    data = veg_1_bpd,
    aes(x = CCA1, y = CCA2, label = veg_1_bpd$env),
    nudge_y = -0.05,
    color = "black",
    size = 3
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
#Save
