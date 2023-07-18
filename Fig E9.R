library(mia)
library(miaViz)
library(phyloseq)
library(tidyverse)

############################FIGURE E9A: MULTIKINGDOM DMM LAPLACE PLOT############################
exp2_combined_tse <- makeTreeSEFromPhyloseq(exp2_combined_prev)

set.seed(1312)
tse_dmn <- mia::runDMN(exp2_combined_tse, name = "DMN", k = 1:7)
getDMN(tse_dmn)

plotDMNFit(tse_dmn, type = "laplace")
#Save


############################FIGURE E9B: MULTIKINGDOM DMM VARIABLE IMPORTANCE PLOT############################

best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
fitted <- head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))

exp2_combined_tax_rel$Genus <- make.unique(exp2_combined_tax_rel$Genus)

for (k in seq(ncol(fitted(best_dmn)))) {
  d <- melt(fitted(best_dmn))
  colnames(d) <- c("OTU", "cluster", "value")
  print(d)
}

exp2_combined_tax_rel_genonly <- data.frame(new_col = c(exp2_combined_tax_rel$Genus, exp2_combined_tax_rel$Genus, exp2_combined_tax_rel$Genus,exp2_combined_tax_rel$Genus))
d <- d %>%
  mutate(Genus = exp2_combined_tax_rel_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.70))
print(d)

###Cluster 1###
d1 <- subset(d, d$cluster == 1)
d1 <- d1 %>% top_n(10, value)
p1 <- ggplot(d1, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 1")) +
  theme_bw()
print(p1)
#Save

###Cluster 2###
d2 <- subset(d, d$cluster == 2)
d2 <- d2 %>% top_n(10, value)
p2 <- ggplot(d2, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 2")) +
  theme_bw()
print(p2)
#Save

###Cluster 3###
d3 <- subset(d, d$cluster == 3)
d3 <- d3 %>% top_n(10, value)
p3 <- ggplot(d3, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 3")) +
  theme_bw()
print(p3)
#Save

###Cluster 4###
d4 <- subset(d, d$cluster ==4)
d4 <- d4 %>% top_n(10, value)
p4 <- ggplot(d4, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 4")) +
  theme_bw()
print(p4)
#Save


############################FIGURE E9C: MULTIKINGDOM DMM CCAS############################

###Get metadata###
exp2_combined_meta_rel <- as.data.frame(exp2_combined_meta_rel)
cols <- c("BPD","BPD_Severity","Died","Sex","BPD_Sex","Feeding","Multiples","IUGR","UTI","Pre.E","CS","NEC...2","Pneumonia","Intubation","C..Sepsis","C..Sepsis.1","Dominance","dmm_component")
exp2_combined_meta_rel[cols] <- lapply(exp2_combined_meta_rel[cols], factor)

###Reorder factors and perform CCA###
sapply(exp2_combined_meta_rel, class)
exp2_combined_meta_rel <- as.data.frame(exp2_combined_meta_rel[-c(64),c(1,3,4,6,7,8,24,25,26,27)])
exp2_combined_meta_rel$BPD<- factor(exp2_combined_meta_rel$BPD, levels=c('PPRD', 'BPD'))
exp2_combined_meta_rel$Died<- factor(exp2_combined_meta_rel$Died, levels=c('Survived', 'Died'))
exp2_combined_meta_rel$Dominance<- factor(exp2_combined_meta_rel$Dominance, levels=c('No_Fungi', 'Low_Evenness','High_Evenness'))
exp2_combined_meta_rel$dmm_component<- factor(exp2_combined_meta_rel$dmm_component, levels=c('comp2', 'comp1', 'comp3', 'comp4'))
exp2_combined_otu_rel <- exp2_combined_otu_rel[-c(64),]
cca_result <- cca(exp2_combined_otu_rel, exp2_combined_meta_rel)

###Check results###
summary(cca_result)
plot(cca_result)

###Extract metadata###
veg_1 = as.data.frame(cca_result$CCA$biplot)
row.names(veg_1) <- c("BPD","Died","AMAB","GA","BW","CRIB.II","Low_Evenness","High_Evenness","Fungi_Abundance","CRP","comp2","comp3","comp4" )
veg_1["env"] = row.names(veg_1)

###Extract samples###
veg_2 = as.data.frame(cca_result$CCA$u)
veg_2["samples"] = row.names(veg_2)

###Extract taxa###
veg_3 = as.data.frame(cca_result$CCA$v)
veg_3["genus"] = row.names(veg_3)

###Plot points###
plot = ggplot() +
  geom_point(data = veg_3, aes(x = CCA1, y = CCA2), color = "#AFEE91", size = 1.0) +
  geom_point(data = veg_2, aes(x = CCA1, y = CCA2), color = "#0B8700", size = 2.0) +
  geom_point(data = veg_1, aes(x = CCA1, y = CCA2), color = "black")


plot

###Plot arrows and labels###
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
    nudge_y = 0.2,
    color = "black",
    size = 3
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
#Save


############################FIGURE E9D: MULTIKINGDOM PPRD DMM LAPLACE PLOT############################

exp2_combined_tse_pprd <- makeTreeSEFromPhyloseq(exp2_combined_prev_pprd)

set.seed(1312)
tse_dmn_pprd <- mia::runDMN(exp2_combined_tse_pprd, name = "DMN", k = 1:7)
getDMN(tse_dmn_pprd)

plotDMNFit(tse_dmn_pprd, type = "laplace")
#Save


############################FIGURE E9E: MULTIKINGDOM PPRD DMM VARIABLE IMPORTANCE PLOT############################

best_dmn_pprd <- getBestDMNFit(tse_dmn_pprd, type = "laplace")

exp2_combined_tax_rel_pprd$Genus <- make.unique(exp2_combined_tax_rel_pprd$Genus)

for (k in seq(ncol(fitted(best_dmn_pprd)))) {
  d_pprd <- melt(fitted(best_dmn_pprd))
  colnames(d_pprd) <- c("OTU", "cluster", "value")
  print(d_pprd)
}

exp2_combined_tax_rel_pprd_genonly <- data.frame(new_col = c(exp2_combined_tax_rel_pprd$Genus, exp2_combined_tax_rel_pprd$Genus, exp2_combined_tax_rel_pprd$Genus))
d_pprd <- d_pprd %>%
  mutate(Genus = exp2_combined_tax_rel_pprd_genonly) %>%
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
  labs(title = paste("Top drivers: community type 1")) +
  theme_bw()
print(p1_pprd)
#Save

###Cluster 2###
d2_pprd <- subset(d_pprd, d_pprd$cluster == 2)
d2_pprd <- d2_pprd %>% top_n(10, value)
p2_pprd <- ggplot(d2_pprd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 2")) +
  theme_bw()
print(p2_pprd)
#Save

###Cluster 3###
d3_pprd <- subset(d_pprd, d_pprd$cluster == 3)
d3_pprd <- d3_pprd %>% top_n(10, value)
p3_pprd <- ggplot(d3_pprd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 3")) +
  theme_bw()
print(p3_pprd)
#Save


############################FIGURE E9F: MULTIKINGDOM PPRD DMM CCAS############################

###Get metadata###
exp2_combined_meta_rel_pprd <- as.data.frame(exp2_combined_meta_rel_pprd)
cols <- c("BPD","BPD_Severity","Died","Sex","BPD_Sex","Feeding","Multiples","IUGR","UTI","Pre.E","CS","NEC...2","Pneumonia","Intubation","C..Sepsis","C..Sepsis.1","Dominance","dmm_component")
exp2_combined_meta_rel_pprd[cols] <- lapply(exp2_combined_meta_rel_pprd[cols], factor)

###Reorder factors and do CCA###
sapply(exp2_combined_meta_rel_pprd, class)
exp2_combined_meta_rel_pprd <- as.data.frame(exp2_combined_meta_rel_pprd[-c(40),c(3,4,6,7,8,24,25,26,27)])
exp2_combined_meta_rel_pprd$Died<- factor(exp2_combined_meta_rel_pprd$Died, levels=c('Survived', 'Died'))
exp2_combined_meta_rel_pprd$Dominance<- factor(exp2_combined_meta_rel_pprd$Dominance, levels=c('No_Fungi', 'Low_Evenness','High_Evenness'))
exp2_combined_meta_rel_pprd$dmm_component<- factor(exp2_combined_meta_rel_pprd$dmm_component, levels=c('comp2','comp1', 'comp3'))
exp2_combined_otu_rel_pprd <- exp2_combined_otu_rel_pprd[-c(40),]
cca_result_pprd <- cca(exp2_combined_otu_rel_pprd, exp2_combined_meta_rel_pprd)

###Summarize results###
summary(cca_result_pprd)
plot(cca_result_pprd)

###Extract metadata###
veg_1_pprd = as.data.frame(cca_result_pprd$CCA$biplot)
row.names(veg_1_pprd) <- c("AMAB","GA","BW","CRIB.II","Low_Evenness","High_Evenness","Fungi_Abundance","CRP","comp1","comp3")
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


############################FIGURE E9G: MULTIKINGDOM BPD DMM LAPLACE PLOT############################

exp2_combined_tse_bpd <- makeTreeSEFromPhyloseq(exp2_combined_prev_bpd)

set.seed(1312)
tse_dmn_bpd <- mia::runDMN(exp2_combined_tse_bpd, name = "DMN", k = 1:7)
getDMN(tse_dmn_bpd)

plotDMNFit(tse_dmn_bpd, type = "laplace")
#Save


############################FIGURE E9H: MULTIKINGDOM BPD DMM VARIABLE IMPORTANCE PLOT############################

best_dmn_bpd <- getBestDMNFit(tse_dmn_bpd, type = "laplace")

DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn_bpd))
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn_bpd)))
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn_bpd)))

exp2_combined_tax_rel_bpd$Genus <- make.unique(exp2_combined_tax_rel_bpd$Genus)

for (k in seq(ncol(fitted(best_dmn_bpd)))) {
  d_bpd <- melt(fitted(best_dmn_bpd))
  colnames(d_bpd) <- c("OTU", "cluster", "value")
  print(d_bpd)
}

exp2_combined_tax_rel_bpd_genonly <- data.frame(new_col = c(exp2_combined_tax_rel_bpd$Genus, exp2_combined_tax_rel_bpd$Genus))
d_bpd <- d_bpd %>%
  mutate(Genus = exp2_combined_tax_rel_bpd_genonly) %>%
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
  labs(title = paste("Top drivers: community type 1")) +
  theme_bw()
print(p1_bpd)
#Save

###Cluster 2###
d2_bpd <- subset(d_bpd, d_bpd$cluster == 2)
d2_bpd <- d2_bpd %>% top_n(10, value)
p2_bpd <- ggplot(d2_bpd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 2")) +
  theme_bw()
print(p2_bpd)
#Save


############################FIGURE E9I: MULTIKINGDOM BPD DMM CCAS############################

###Get metadata###
exp2_combined_meta_rel_bpd <- as.data.frame(exp2_combined_meta_rel_bpd)
cols <- c("BPD","BPD_Severity","Died","Sex","BPD_Sex","Feeding","Multiples","IUGR","UTI","Pre.E","CS","NEC...2","Pneumonia","Intubation","C..Sepsis","C..Sepsis.1","Dominance","dmm_component")
exp2_combined_meta_rel_bpd[cols] <- lapply(exp2_combined_meta_rel_bpd[cols], factor)

###Reorder factors and perform CCA###
sapply(exp2_combined_meta_rel_bpd, class)
exp2_combined_meta_rel_bpd <- as.data.frame(exp2_combined_meta_rel_bpd[,c(1,3,4,6,7,8,24,25,26,27)])
exp2_combined_meta_rel_bpd$BPD<- factor(exp2_combined_meta_rel_bpd$BPD, levels=c('PPRD', 'BPD'))
exp2_combined_meta_rel_bpd$Died<- factor(exp2_combined_meta_rel_bpd$Died, levels=c('Survived', 'Died'))
exp2_combined_meta_rel_bpd$Dominance<- factor(exp2_combined_meta_rel_bpd$Dominance, levels=c('No_Fungi', 'Low_Evenness','High_Evenness'))
exp2_combined_meta_rel_bpd$dmm_component<- factor(exp2_combined_meta_rel_bpd$dmm_component, levels=c('comp1','comp2'))
exp2_combined_otu_rel_bpd <- exp2_combined_otu_rel_bpd
cca_result_bpd <- cca(exp2_combined_otu_rel_bpd, exp2_combined_meta_rel_bpd)

summary(cca_result_bpd)
plot(cca_result_bpd)

###Extract metadata###
veg_1_bpd = as.data.frame(cca_result_bpd$CCA$biplot)
row.names(veg_1_bpd) <- c("Died","AMAB","GA","BW","CRIB.II","Low_Evenness","High_Evenness","Fungi_Abundance","CRP","comp2")
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
