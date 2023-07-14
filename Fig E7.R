library(mia)
library(miaViz)
library(DirichletMultinomial)
library(reshape2)

############################FIGURE E7A: ALL ITS DMM LAPLACE PLOT############################
exp2_fung_tse <- makeTreeSEFromPhyloseq(exp2_fung_prev)

set.seed(1312)
tse_dmn <- mia::runDMN(exp2_fung_tse, name = "DMN", k = 1:7)
getDMN(tse_dmn)

###Since 2 and 5 clusters had similar values, we ran both###
###Fitting model for 5 clusters###
plotDMNFit(tse_dmn, type = "laplace")
best_dmn <- getBestDMNFit(tse_dmn, type = "laplace")

###Fitting model for 2 clusters###
tse_dmn2 <- mia::runDMN(exp2_fung_tse, name = "DMN", k = 2)
getDMN(tse_dmn2)
best_dmn2 <- getBestDMNFit(tse_dmn2, type = "laplace")

############################FIGURE E7B: ALL ITS DMM 2 CLUSTERS VARIABLE IMPORTANCE############################

###Identifies taxa present in each cluster###
for (k in seq(ncol(fitted(best_dmn2)))) {
  d2 <- melt(fitted(best_dmn2))
  colnames(d2) <- c("OTU", "cluster", "value")
  print(d2)
}

###Gets the most important variables (>0.80) and arranges them in descending order###
exp2_fung_tax_rel_genonly2 <- data.frame(new_col = c(exp2_fung_tax_rel$Genus, exp2_fung_tax_rel$Genus))
d2 <- d2 %>%
  mutate(Genus = exp2_fung_tax_rel_genonly) %>%
  arrange(value) %>%
  mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
  # Only show the most important drivers
  filter(abs(value) > quantile(abs(value), 0.80))
print(d2)

###Plots top variables for Cluster 1###
d2_1 <- subset(d2, d2$cluster == 1)
d2_1 <- d2_1 %>% top_n(10, value)
p2_1 <- ggplot(d2_1, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 1")) +
  theme_bw()
print(p2_1)
#Save

###Plots top variables for Cluster 2###
d2_2 <- subset(d2, d2$cluster == 2)
d2_2 <- d2_2 %>% top_n(10, value)
p2_2 <- ggplot(d2_2, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 2")) +
  theme_bw()
print(p2_2)
#Save


############################FIGURE E7C: ALL ITS DMM 2 CLUSTERS CCA############################

###Create a dataframe of all metadata factors to plot###
exp2_fung_meta_rel2 <- as.data.frame(exp2_fung_meta_rel2)
cols <- c("BPD","BPD_Severity","Died","Sex","BPD_Sex","Feeding","Multiples","IUGR","UTI","Pre.E","CS","NEC...2","Pneumonia","Intubation","C..Sepsis","C..Sepsis.1","Dominance","dmm_component")
exp2_fung_meta_rel2[cols] <- lapply(exp2_fung_meta_rel2[cols], factor)

#Reorder factors as needed so arrows show the main variables of interest
sapply(exp2_fung_meta_rel2, class)
exp2_fung_meta_rel2 <- as.data.frame(exp2_fung_meta_rel2[,c(1,3,4,6,7,8,24,25,26,28)])
exp2_fung_meta_rel2$BPD<- factor(exp2_fung_meta_rel2$BPD, levels=c('PPRD', 'BPD'))
exp2_fung_meta_rel2$Died<- factor(exp2_fung_meta_rel2$Died, levels=c('Survived', 'Died'))
exp2_fung_meta_rel2$dmm_component<- factor(exp2_fung_meta_rel2$dmm_component, levels=c('comp1','comp2'))
cca_result2 <- cca(exp2_fung_otu_rel, exp2_fung_meta_rel2)

#summary of the data
summary(cca_result2)
plot(cca_result2)

#extract metadata for arrows
veg_12 = as.data.frame(cca_result2$CCA$biplot)
row.names(veg_12) <- c("BPD","Died","AMAB","GA","BW","CRIB.II","Low_Evenness","Fungi_Abundance","CRP","comp2")
veg_12["env"] = row.names(veg_12)

#extract samples for dark green points
veg_22 = as.data.frame(cca_result2$CCA$u)
veg_22["samples"] = row.names(veg_22)

#extract taxa for light green points
veg_32 = as.data.frame(cca_result2$CCA$v)
veg_32["genus"] = row.names(veg_32)

#Plots points only
plot2 = ggplot() +
  geom_point(data = veg_32, aes(x = CCA1, y = CCA2), color = "#AFEE91", size = 1.0) +
  geom_point(data = veg_22, aes(x = CCA1, y = CCA2), color = "#0B8700", size = 2.0) +
  geom_point(data = veg_12, aes(x = CCA1, y = CCA2), color = "black")


plot2

#Adds arrows and labels, ensures labels overlap as little as possible
plot2 +
  theme_bw() +
  geom_segment(
    data = veg_12,
    aes(
      x = 0,
      y = 0,
      xend = CCA1,
      yend = CCA2
    ),
    arrow = arrow(length = unit(0.5, "cm"))
  ) +
  geom_text_repel(
    data = veg_12,
    aes(x = CCA1, y = CCA2, label = veg_12$env),
    nudge_y = -0.05,
    color = "black",
    size = 3
  ) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
#Save


############################FIGURE E7D: ALL ITS DMM 5 CLUSTERS VARIABLE IMPORTANCE############################

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


############################FIGURE E7E: ALL ITS DMM 5 CLUSTERS CCA############################

#Get clinical data
exp2_fung_meta_rel <- as.data.frame(exp2_fung_meta_rel)
cols <- c("BPD","BPD_Severity","Died","Sex","BPD_Sex","Feeding","Multiples","IUGR","UTI","Pre.E","CS","NEC...2","Pneumonia","Intubation","C..Sepsis","C..Sepsis.1","Dominance","cluster")
exp2_fung_meta_rel[cols] <- lapply(exp2_fung_meta_rel[cols], factor)

#Reorder factors
sapply(exp2_fung_meta_rel, class)
exp2_fung_meta_rel <- as.data.frame(exp2_fung_meta_rel[,c(1,3,4,6,7,8,24,25,26,27)])
exp2_fung_meta_rel$BPD<- factor(exp2_fung_meta_rel$BPD, levels=c('PPRD', 'BPD'))
exp2_fung_meta_rel$Died<- factor(exp2_fung_meta_rel$Died, levels=c('Survived', 'Died'))
exp2_fung_meta_rel$cluster<- factor(exp2_fung_meta_rel$cluster, levels=c('Cluster 4', 'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 5'))
cca_result <- cca(exp2_fung_otu_rel, exp2_fung_meta_rel)

#summary of the data
summary(cca_result)
plot(cca_result)

#extracting data: clinical metadata
veg_1 = as.data.frame(cca_result$CCA$biplot)
row.names(veg_1) <- c("BPD","Died","AMAB","GA","BW","CRIB.II","Low_Evenness","Fungi_Abundance","CRP","comp1","comp2","comp3","comp5" )
veg_1["env"] = row.names(veg_1)

#extracting data: samples
veg_2 = as.data.frame(cca_result$CCA$u)
veg_2["samples"] = row.names(veg_2)

#extracting data: taxa
veg_3 = as.data.frame(cca_result$CCA$v)
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


############################FIGURE E7F: BPD ONLY ITS DMM LAPLACE PLOT############################

exp2_fung_tse_bpd <- makeTreeSEFromPhyloseq(exp2_fung_prev_bpd)

#Get best k value
set.seed(1312)
tse_dmn_bpd <- mia::runDMN(exp2_fung_tse_bpd, name = "DMN", k = 1:7)
getDMN(tse_dmn_bpd)

#Laplace plot
plotDMNFit(tse_dmn_bpd, type = "laplace")
#DMM with k=3
best_dmn_bpd <- getBestDMNFit(tse_dmn_bpd, type = "laplace")


############################FIGURE E7G: BPD ONLY ITS DMM VARIABLE IMPORTANCE############################

exp2_fung_tax_rel_bpd$Genus <- make.unique(exp2_fung_tax_rel_bpd$Genus)

#Get fitted values
for (k in seq(ncol(fitted(best_dmn_bpd)))) {
  d_bpd <- melt(fitted(best_dmn_bpd))
  colnames(d_bpd) <- c("OTU", "cluster", "value")
  print(d_bpd)
}

#Get top fitted values
exp2_fung_tax_rel_bpd_genonly <- data.frame(new_col = c(exp2_fung_tax_rel_bpd$Genus, exp2_fung_tax_rel_bpd$Genus, exp2_fung_tax_rel_bpd$Genus))
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

#Plot Cluster 3
d3_bpd <- subset(d_bpd, d_bpd$cluster == 3)
d3_bpd <- d3_bpd %>% top_n(10, value)
p3_bpd <- ggplot(d3_bpd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 3")) +
  theme_bw()
print(p3_bpd)
#Save


############################FIGURE E7H: BPD ONLY ITS DMM CCA############################

#Get clinical data
exp2_fung_meta_rel_bpd <- as.data.frame(exp2_fung_meta_rel_bpd)
cols <- c("BPD","BPD_Severity","Died","Sex","BPD_Sex","Feeding","Multiples","IUGR","UTI","Pre.E","CS","NEC...2","Pneumonia","Intubation","C..Sepsis","C..Sepsis.1","Dominance","dmm_component")
exp2_fung_meta_rel_bpd[cols] <- lapply(exp2_fung_meta_rel_bpd[cols], factor)

#Reorder factors
sapply(exp2_fung_meta_rel_bpd, class)
exp2_fung_meta_rel_bpd <- as.data.frame(exp2_fung_meta_rel_bpd[,c(1,3,4,6,7,8,24,25,26,27)])
exp2_fung_meta_rel_bpd$BPD<- factor(exp2_fung_meta_rel_bpd$BPD, levels=c('PPRD', 'BPD'))
exp2_fung_meta_rel_bpd$Died<- factor(exp2_fung_meta_rel_bpd$Died, levels=c('Survived', 'Died'))
exp2_fung_meta_rel_bpd$dmm_component<- factor(exp2_fung_meta_rel_bpd$dmm_component, levels=c('Cluster 3','Cluster 1', 'Cluster 2'))
cca_result_bpd <- cca(exp2_fung_otu_rel_bpd, exp2_fung_meta_rel_bpd)

#summary of the data
summary(cca_result_bpd)
plot(cca_result_bpd)

#Extract metadata
veg_1_bpd = as.data.frame(cca_result_bpd$CCA$biplot)
row.names(veg_1_bpd) <- c("Died","AMAB","GA","BW","CRIB.II","Low_Evenness","Fungi_Abundance","CRP","comp1","comp2")
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


############################FIGURE E7I: PPRD ONLY ITS DMM LAPLACE PLOT############################

exp2_fung_tse_pprd <- makeTreeSEFromPhyloseq(exp2_fung_prev_pprd)

#Get best k value
set.seed(1312)
tse_dmn_pprd <- mia::runDMN(exp2_fung_tse_pprd, name = "DMN", k = 1:7)
getDMN(tse_dmn_pprd)

#Laplace plot
plotDMNFit(tse_dmn_pprd, type = "laplace")
#Best model, k=3
best_dmn_pprd <- getBestDMNFit(tse_dmn_pprd, type = "laplace")


############################FIGURE E7J: PPRD ONLY ITS DMM VARIABLE IMPORTANCE############################

exp2_fung_tax_rel_pprd$Genus <- make.unique(exp2_fung_tax_rel_pprd$Genus)

#Get fitted values
for (k in seq(ncol(fitted(best_dmn_pprd)))) {
  d_pprd <- melt(fitted(best_dmn_pprd))
  colnames(d_pprd) <- c("OTU", "cluster", "value")
  print(d_pprd)
}

#Get top fitted values
exp2_fung_tax_rel_pprd_genonly <- data.frame(new_col = c(exp2_fung_tax_rel_pprd$Genus, exp2_fung_tax_rel_pprd$Genus, exp2_fung_tax_rel_pprd$Genus))
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

#Plot Cluster 3
d3_pprd <- subset(d_pprd, d_pprd$cluster == 3)
d3_pprd <- d3_pprd %>% top_n(10, value)
p3_pprd <- ggplot(d3_pprd, aes(x = reorder(Genus$new_col, value), y = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Top drivers: community type 3")) +
  theme_bw()
print(p3_pprd)
#Save


############################FIGURE E7K: PPRD ONLY ITS DMM CCA############################

#Metadata
exp2_fung_meta_rel_pprd <- as.data.frame(exp2_fung_meta_rel_pprd)
cols <- c("BPD","BPD_Severity","Died","Sex","BPD_Sex","Feeding","Multiples","IUGR","UTI","Pre.E","CS","NEC...2","Pneumonia","Intubation","C..Sepsis","C..Sepsis.1","Dominance","dmm_component")
exp2_fung_meta_rel_pprd[cols] <- lapply(exp2_fung_meta_rel_pprd[cols], factor)

#Reorder metadata factors
sapply(exp2_fung_meta_rel_pprd, class)
exp2_fung_meta_rel_pprd <- as.data.frame(exp2_fung_meta_rel_pprd[,c(3,4,6,7,8,24,25,26,27)])
exp2_fung_meta_rel_pprd$Died<- factor(exp2_fung_meta_rel_pprd$Died, levels=c('Survived', 'Died'))
exp2_fung_meta_rel_pprd$dmm_component<- factor(exp2_fung_meta_rel_pprd$dmm_component, levels=c('Cluster 3','Cluster 1', 'Cluster 2'))
cca_result_pprd <- cca(exp2_fung_otu_rel_pprd, exp2_fung_meta_rel_pprd)

#summary
summary(cca_result_pprd)
plot(cca_result_pprd)

#Extract metadata
veg_1_pprd = as.data.frame(cca_result_pprd$CCA$biplot)
row.names(veg_1_pprd) <- c("AMAB","GA","BW","CRIB.II","Low_Evenness","Fungi_Abundance","CRP","comp1","comp2")
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