library(mia)
library(miaViz)
library(DirichletMultinomial)
library(reshape2)
library(ggrepel)


############################FIGURE S7A: ALL ITS DMM LAPLACE PLOT############################
exp2_fung_tse <- makeTreeSEFromPhyloseq(exp2_fung_prev)

set.seed(1312)
tse_dmn <- mia::runDMN(exp2_fung_tse, name = "DMN", k = 1:7)
getDMN(tse_dmn)

###Since 2 and 5 clusters had similar values, we ran both###
###Fitting model for 5 clusters###
plotDMNFit(tse_dmn, type = "laplace")
#Save
best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")


############################FIGURE S7B: ALL ITS DMM VARIABLE IMPORTANCE############################

#Get fitted values
for (k in seq(ncol(fitted(best_dmn)))) {
  d <- melt(fitted(best_dmn))
  colnames(d) <- c("OTU", "cluster", "value")
  print(d)
}

###Get important variables and arrange them in descending order###
exp2_fung_tax_rel_genonly <- data.frame(new_col = c(exp2_fung_tax_rel$Genus, exp2_fung_tax_rel$Genus, exp2_fung_tax_rel$Genus, exp2_fung_tax_rel$Genus, exp2_fung_tax_rel$Genus))
d <- d %>%
  mutate(Genus = exp2_fung_tax_rel_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.80))
print(d)

###Plot Cluster 1###
d1 <- subset(d, d$cluster == 1)
d1 <- d1 %>% top_n(10, value)
p1 <- ggplot(d1, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 1")) +
  theme_bw()
print(p1)
#Save

###Plot Cluster 2###
d2 <- subset(d, d$cluster == 2)
d2 <- d2 %>% top_n(10, value)
p2 <- ggplot(d2, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 2")) +
  theme_bw()
print(p2)
#Save

###Plot Cluster 3###
d3 <- subset(d, d$cluster == 3)
d3 <- d3 %>% top_n(10, value)
p3 <- ggplot(d3, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 3")) +
  theme_bw()
print(p3)
#Save

###Plot Cluster 4###
d4 <- subset(d, d$cluster ==4)
d4 <- d4 %>% top_n(10, value)
p4 <- ggplot(d4, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 4")) +
  theme_bw()
print(p4)
#Save

###Plot Cluster 5###
d5 <- subset(d, d$cluster == 5)
d5 <- d5 %>% top_n(10, value)
p5 <- ggplot(d5, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 5")) +
  theme_bw()
print(p5)
#Save


############################FIGURE S7C: ALL ITS DMM CCA############################

#Read in metadata file with clinical data (not provided on GitHub)

md_fung_newfactors_sort$cluster <- vec
exp2_fung_meta_df <- as.data.frame(md_fung_newfactors_sort[,c(1,3,4,6,7,8,25,28,29,30)])
exp2_fung_meta_df[sapply(exp2_fung_meta_df, is.character)] <- lapply(exp2_fung_meta_df[sapply(exp2_fung_meta_df, is.character)], as.factor)
sapply(exp2_fung_meta_df, base::class)

exp2_fung_meta_df$BPD<- factor(exp2_fung_meta_df$BPD, levels=c('No_BPD', 'BPD'))
exp2_fung_meta_df$Died<- factor(exp2_fung_meta_df$Died, levels=c('Survived', 'Died'))
exp2_fung_meta_df$cluster<- factor(exp2_fung_meta_df$cluster, levels=c('Cluster 4','Cluster 2','Cluster 3','Cluster 1','Cluster 5'))

cca_result2 <- vegan::cca(exp2_fung_otu_rel,exp2_fung_meta_df)


#summary of the data
summary(cca_result2)
plot(cca_result2)

#extracting data: clinical metadata
veg_1 = as.data.frame(cca_result2$CCA$biplot)
row.names(veg_1) <- c("BPD","Died","AMAB","GA","BW","CRIB.II","CRP","Low_Evenness","Fungi_Prev","Cluster 2","Cluster 3","Cluster 1","Cluster 5")
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


############################FIGURE S7D: NOBPD ONLY ITS DMM LAPLACE PLOT############################

exp2_fung_tse_pprd <- makeTreeSEFromPhyloseq(exp2_fung_prev_pprd)

set.seed(1312)
tse_dmn_pprd <- mia::runDMN(exp2_fung_tse_pprd, name = "DMN", k = 1:7)
getDMN(tse_dmn_pprd)

plotDMNFit(tse_dmn_pprd, type = "laplace")
#Save

best_dmn_pprd <- getBestDMNFit(tse_dmn_pprd, type = "laplace")


############################FIGURE S7E: NOBPD ONLY ITS DMM VARIABLE IMPORTANCE############################

exp2_fung_tax_rel_pprd$Genus <- make.unique(exp2_fung_tax_rel_pprd$Genus)

#Get fitted values
for (k in seq(ncol(fitted(best_dmn_pprd)))) {
  d_pprd <- melt(fitted(best_dmn_pprd))
  colnames(d_pprd) <- c("OTU", "cluster", "value")
  print(d_pprd)
}

#Get top fitted values
exp2_fung_tax_rel_pprd_genonly <- data.frame(new_col = c(exp2_fung_tax_rel_pprd$Genus, exp2_fung_tax_rel_pprd$Genus))
d_pprd <- d_pprd %>%
  mutate(Genus = exp2_fung_tax_rel_pprd_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.70))
print(d_pprd)

#Plot Cluster 1
d1_pprd <- subset(d_pprd, d_pprd$cluster == 1)
d1_pprd <- d1_pprd %>% top_n(10, value)
p1_pprd <- ggplot(d1_pprd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 1")) +
  theme_bw()
print(p1_pprd)
#Save

#Plot Cluster 2
d2_pprd <- subset(d_pprd, d_pprd$cluster == 2)
d2_pprd <- d2_pprd %>% top_n(10, value)
p2_pprd <- ggplot(d2_pprd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 2")) +
  theme_bw()
print(p2_pprd)
#Save


############################FIGURE S7F: NOBPD ONLY ITS DMM CCA############################

md_fung_newfactors_sort_pprd <- md_fung_newfactors_sort[md_fung_newfactors_sort$BPD == "No_BPD",]
md_fung_newfactors_sort_pprd$cluster <- vec_pprd
exp2_fung_meta_df_pprd <- as.data.frame(md_fung_newfactors_sort_pprd[,c(1,3,4,6,7,8,25,28,29,30)])
exp2_fung_meta_df_pprd[sapply(exp2_fung_meta_df_pprd, is.character)] <- lapply(exp2_fung_meta_df_pprd[sapply(exp2_fung_meta_df_pprd, is.character)], as.factor)
sapply(exp2_fung_meta_df, base::class)

exp2_fung_meta_df_pprd$BPD<- factor(exp2_fung_meta_df_pprd$BPD, levels=c('No_BPD', 'BPD'))
exp2_fung_meta_df_pprd$Died<- factor(exp2_fung_meta_df_pprd$Died, levels=c('Survived', 'Died'))
exp2_fung_meta_df_pprd$cluster<- factor(exp2_fung_meta_df_pprd$cluster, levels=c('Cluster 1', 'Cluster 2'))
cca_result_pprd <- cca(exp2_fung_otu_rel_pprd, exp2_fung_meta_df_pprd)

#summary
summary(cca_result_pprd)
plot(cca_result_pprd)

#Extract metadata
veg_1_pprd = as.data.frame(cca_result_pprd$CCA$biplot)
row.names(veg_1_pprd) <- c("AMAB","GA","BW","CRIB.II","CRP","Low_Evenness","Fungi_Abundance","Cluster_1")
veg_1_pprd["env"] = row.names(veg_1_pprd)

#Extract sample data
veg_2_pprd = as.data.frame(cca_result_pprd$CCA$u)
veg_2_pprd["samples"] = row.names(veg_2_pprd)

#Extract taxa
veg_3_pprd = as.data.frame(cca_result_pprd$CCA$v)
veg_3_pprd["genus"] = row.names(veg_3_pprd)

#Plot points
plot_pprd = ggplot() +
  geom_point(data = veg_3_pprd, aes(x = CCA1, y = CCA2), color = "#AFEE91", size = 1.0) +
  geom_point(data = veg_2_pprd, aes(x = CCA1, y = CCA2), color = "#0B8700", size = 2.0) +
  geom_point(data = veg_1_pprd, aes(x = CCA1, y = CCA2), color = "black")


plot_pprd

#Plot arrows and labels
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


############################FIGURE S7F: BPD ONLY ITS DMM LAPLACE PLOT############################

exp2_fung_tse_bpd <- makeTreeSEFromPhyloseq(exp2_fung_prev_bpd)

set.seed(1312)
tse_dmn_bpd <- mia::runDMN(exp2_fung_tse_bpd, name = "DMN", k = 1:7)
getDMN(tse_dmn_bpd)

plotDMNFit(tse_dmn_bpd, type = "laplace")
#Save

best_dmn_bpd <- getBestDMNFit(tse_dmn_bpd, type = "laplace")


############################FIGURE S7G: BPD ONLY ITS DMM VARIABLE IMPORTANCE############################

exp2_fung_tax_rel_bpd$Genus <- make.unique(exp2_fung_tax_rel_bpd$Genus)

#Get fitted values
for (k in seq(ncol(fitted(best_dmn_bpd)))) {
  d_bpd <- melt(fitted(best_dmn_bpd))
  colnames(d_bpd) <- c("OTU", "cluster", "value")
  print(d_bpd)
}

#Get top fitted values
exp2_fung_tax_rel_bpd_genonly <- data.frame(new_col = c(exp2_fung_tax_rel_bpd$Genus, exp2_fung_tax_rel_bpd$Genus))
d_bpd <- d_bpd %>%
  mutate(Genus = exp2_fung_tax_rel_bpd_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.70))
print(d_bpd)

#Plot Cluster 1
d1_bpd <- subset(d_bpd, d_bpd$cluster == 1)
d1_bpd <- d1_bpd %>% top_n(10, value)
p1_bpd <- ggplot(d1_bpd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 1")) +
  theme_bw()
print(p1_bpd)
#Save

#Plot Cluster 2
d2_bpd <- subset(d_bpd, d_bpd$cluster == 2)
d2_bpd <- d2_bpd %>% top_n(10, value)
p2_bpd <- ggplot(d2_bpd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 2")) +
  theme_bw()
print(p2_bpd)
#Save


############################FIGURE S7H: BPD ONLY ITS DMM CCA############################

md_fung_newfactors_sort_bpd <- md_fung_newfactors_sort[md_fung_newfactors_sort$BPD == "BPD",]
md_fung_newfactors_sort_bpd$cluster <- vec_bpd
exp2_fung_meta_df_bpd <- as.data.frame(md_fung_newfactors_sort_bpd[,c(1,3,4,6,7,8,25,28,29,30)])
exp2_fung_meta_df_bpd[sapply(exp2_fung_meta_df_bpd, is.character)] <- lapply(exp2_fung_meta_df_bpd[sapply(exp2_fung_meta_df_bpd, is.character)], as.factor)
sapply(exp2_fung_meta_df, base::class)

exp2_fung_meta_df_bpd$BPD<- factor(exp2_fung_meta_df_bpd$BPD, levels=c('No_BPD', 'BPD'))
exp2_fung_meta_df_bpd$Died<- factor(exp2_fung_meta_df_bpd$Died, levels=c('Survived', 'Died'))
exp2_fung_meta_df_bpd$cluster<- factor(exp2_fung_meta_df_bpd$cluster, levels=c('Cluster 2', 'Cluster 1'))
cca_result_bpd <- cca(exp2_fung_otu_rel_bpd, exp2_fung_meta_df_bpd)

#summary of the data
summary(cca_result_bpd)
plot(cca_result_bpd)

#Extract metadata
veg_1_bpd = as.data.frame(cca_result_bpd$CCA$biplot)
row.names(veg_1_bpd) <- c("Died","AMAB","GA","BW","CRIB.II","CRP","Low_Evenness","Fungi_Abundance","Cluster_1")
veg_1_bpd["env"] = row.names(veg_1_bpd)

#Extract sample data
veg_2_bpd = as.data.frame(cca_result_bpd$CCA$u)
veg_2_bpd["samples"] = row.names(veg_2_bpd)

#Extract taxa
veg_3_bpd = as.data.frame(cca_result_bpd$CCA$v)
veg_3_bpd["genus"] = row.names(veg_3_bpd)

#Plot points
plot_bpd = ggplot() +
  geom_point(data = veg_3_bpd, aes(x = CCA1, y = CCA2), color = "#AFEE91", size = 1.0) +
  geom_point(data = veg_2_bpd, aes(x = CCA1, y = CCA2), color = "#0B8700", size = 2.0) +
  geom_point(data = veg_1_bpd, aes(x = CCA1, y = CCA2), color = "black")


plot_bpd

#Plot arrows and labels
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
