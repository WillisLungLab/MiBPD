library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(randomForest)
library(xgboost)
library(e1071)
library(caret)
library(glmnet)
library(MASS)
library(pROC)
library(stringr)
library(MLeval)

########ITS Random Forest########
exp2_fung_gen_prev_noNA = subset_taxa(exp2_fung_gen_prev, Genus!="NA")
exp2_fung_gen_prev_rel <- transform_sample_counts(exp2_fung_gen_prev_noNA, function(x) x / sum(x) )

###Prep OTU table###
predictors <- as.data.frame(exp2_fung_gen_prev_rel@otu_table)
taxlabels <- as.data.frame(exp2_fung_gen_prev_rel@tax_table)
predictors <- t(predictors)
colnames(predictors) <- taxlabels$Genus
#Get relevant metadata
response <- as.factor(sample_data(exp2_fung_gen_prev_rel)$BPD)

###Make dataframe combining OTUs and Metadata###
rf.data <- data.frame(response, predictors)

###Set random seed and partition data into training and testing###
set.seed(1312)
fung_idx = createDataPartition(rf.data$response, p = 0.75, list = FALSE)
fung_train = rf.data[fung_idx, ]
fung_test = rf.data[-fung_idx, ]

###Set cross-validation, in this case Leave One Out Cross-Validation###
fung_cv <- trainControl(method='LOOCV', classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

#Weight model so both groups are weighted equally
model_weights <- ifelse(fung_train$response == "BPD",
                        (1/table(fung_train$response)[1]) * 0.5,
                        (1/table(fung_train$response)[2]) * 0.5
)

###Set parameter tuning grid, testing a few parameter combinations###
tgrid <- expand.grid(
  .mtry = c(2,16,32,63),
  .splitrule = c("gini","extratrees"),
  .min.node.size = c(1,3,5,10,18)
)

###Run random forest###
set.seed(1312)
rf_fit <- train(as.factor(response) ~ ., 
                data = fung_train, 
                method = "ranger",
                tuneGrid = tgrid,
                num.trees = 1000,
                trControl = fung_cv,
                metric = "ROC",
                verbose = FALSE,
                verbosity = 0,
                importance = "permutation"
)

rf_fit

###New grid, using the best hyperparameter combination###
best_grid <- expand.grid(
  .mtry = 63,
  .splitrule = "gini",
  .min.node.size = 18
)

###Final random forest model###
set.seed(1312)
rf_fit2 <- train(as.factor(response) ~ ., 
                 data = fung_train, 
                 method = "ranger",
                 tuneGrid = best_grid,
                 num.trees = 1000,
                 trControl = fung_cv,
                 metric = "ROC",
                 verbose = FALSE,
                 verbosity = 0,
                 importance = "permutation"
)

rf_fit2


############################FIGURE S4A: ITS RANDOM FOREST VARIABLE IMPORTANCE############################
features <- varImp(rf_fit2)
features <- features$importance
features <- rownames_to_column(features)
features <- features[order(-features$Overall),]
features_top25 <- as.data.frame(features[1:25,])

ggplot2::ggplot(features_top25, aes(x=reorder(rowname, Overall), y=Overall)) +
  geom_point( color="#FB0207", size=4, alpha=1.0)+
  geom_segment( aes(x=rowname, xend=rowname, y=0, yend=Overall), 
                color='#EF94A2') +
  xlab('Variable')+
  ylab('Overall Importance')+
  theme_light()+
  coord_flip()
#Save


############################FIGURE S4B: ITS RANDOM FOREST ROC############################

res_rf <- evalm(rf_fit2)
#Save

###Check against test data, don't worry too much about overfitting, this is exploratory so the confusion matrix is mostly a sanity check###
testclass <- predict(rf_fit2, newdata = fung_test)
cfMatrix <- confusionMatrix(data = testclass, fung_test$response)



########16s Random Forest########
exp2_bact_gen_prev_noNA = subset_taxa(exp2_bact_gen_prev, Genus!="NA")
exp2_bact_gen_prev_rel <- transform_sample_counts(exp2_bact_gen_prev_noNA, function(x) x / sum(x) )

###Prep OTU table###
predictors <- as.data.frame(exp2_bact_gen_prev_rel@otu_table)
taxlabels <- as.data.frame(exp2_bact_gen_prev_rel@tax_table)
predictors <- t(predictors)
colnames(predictors) <- taxlabels$Genus
#Get predictors from clinical data
#Get relevant metadata
response <- as.factor(sample_data(exp2_bact_gen_prev_rel)$BPD)
#Surprise tool that will help us later
tax.labels <- as.data.frame(exp2_bact_gen_prev_rel@tax_table)

###Make dataframe combining OTUs and Metadata###
rf.data <- data.frame(response, predictors)

###Set random seed and partition data into training and testing###
set.seed(1312)
bact_idx = createDataPartition(rf.data$response, p = 0.75, list = FALSE)
bact_train = rf.data[bact_idx, ]
bact_test = rf.data[-bact_idx, ]

###Set cross-validation, in this case Leave One Out Cross-Validation###
bact_cv <- trainControl(method='LOOCV', classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)


###Set parameter tuning grid, testing a few parameter combinations###
tgrid <- expand.grid(
  .mtry = c(2,16,32,63),
  .splitrule = c("gini","extratrees"),
  .min.node.size = c(1,3,5,10,18)
)

###Run random forest###
set.seed(1312)
rf_fit <- train(as.factor(response) ~ ., 
                data = bact_train, 
                method = "ranger",
                tuneGrid = tgrid,
                num.trees = 1000,
                trControl = bact_cv,
                metric = "ROC",
                verbose = FALSE,
                verbosity = 0,
                importance = "permutation"
)

rf_fit

###New grid, using the best hyperparameter combination###
best_grid <- expand.grid(
  .mtry = 16,
  .splitrule = "gini",
  .min.node.size = 1
)

###Final random forest model###
set.seed(1312)
rf_fit2 <- train(as.factor(response) ~ ., 
                 data = bact_train, 
                 method = "ranger",
                 tuneGrid = best_grid,
                 num.trees = 1000,
                 trControl = bact_cv,
                 metric = "ROC",
                 verbose = FALSE,
                 verbosity = 0,
                 importance = "permutation"
)

rf_fit2

############################FIGURE S4C: 16S RANDOM FOREST VARIABLE IMPORTANCE############################
features <- varImp(rf_fit2)
features <- features$importance
features <- rownames_to_column(features)
features <- features[order(-features$Overall),]
features_top50 <- as.data.frame(features[1:50,])
features_top25 <- as.data.frame(features[1:25,])

ggplot2::ggplot(features_top25, aes(x=reorder(rowname, Overall), y=Overall)) +
  geom_point( color="#FB0207", size=4, alpha=1.0)+
  geom_segment( aes(x=rowname, xend=rowname, y=0, yend=Overall), 
                color='#EF94A2') +
  xlab('Variable')+
  ylab('Overall Importance')+
  theme_light() +
  coord_flip() 
#Save


############################FIGURE S4D: 16S RANDOM FOREST VARIABLE IMPORTANCE############################

res_rf <- evalm(rf_fit2)
#Save

###Check against test data###
testclass <- predict(rf_fit2, newdata = bact_test)
cfMatrix <- confusionMatrix(data = testclass, bact_test$response)
