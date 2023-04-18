library(tidyverse)
library(rpart)
library(rpart.plot)
library(caret)
library(zoo)
library(ggbiplot)
library(reshape2)
library(mlr)
library(randomForest)
library(readxl)
library(ggrepel)

dataset_orig <- read_excel("ExractE3fromAll_Dups_Fixed_SourcesMerged_Apr11_ErG.xlsx", sheet = "ExractE3fromAll_Dups_Fixed_Sour")

names(dataset_orig) <- make.names(names(dataset_orig))

dataset_orig2 <- dataset_orig[,-c(2:3)]

dataset_orig2 <- dataset_orig2 %>% select(id, E3.Ligase, source, Class, everything())

dataset_pca <- dataset_orig2[,-c(1:3)]

toremove <- nearZeroVar(dataset_pca[,-c(1)])

toremove <- toremove + 1

dataset_pca <- dataset_pca[,-toremove]

############################################################
#                                                          #
#                           PCA                            #
#                                                          #
############################################################

pca_ErG <- prcomp(dataset_pca[,-c(1)], scale. = T, center = T)

pca_ErG_df <- as.data.frame(pca_ErG$x[,c(1,2)], stringsAsFactors = F)
pca_ErG_df$source <- dataset_orig$source
pca_ErG_df$Class <- dataset_orig$Class
pca_ErG_df$ID <- dataset_orig$id

p <- ggbiplot(pca_ErG,  groups = pca_ErG_df$Class, circle = F, ellipse = T, var.axes = F)#, labels = pca_ErG_df$ID)
p
p + geom_point(aes(color = as.factor(pca_ErG_df$Class)), size = 3) +
  theme_light(13)+
  theme(legend.position = "top")+
  guides(color=guide_legend(title="E3 Ligase"))#+


ggplot(pca_ErG_df, aes(PC1, PC2))+
  geom_point(aes(color = as.factor(pca_ErG_df$Class)), size = 3)+
  theme(legend.position = "top")+
  guides(color=guide_legend(title="E3 Ligase"))#+

#pca_ErG_df2 <- pca_ErG_df %>% mutate(ID = ifelse(PC1 > -6 | PC2 > 0 | PC2 < -2 , "", ID)) 

# ggplot(pca_ErG_df2, aes(PC1, PC2))+
#   geom_point(aes(color = as.factor(pca_ErG_df2$Class)), size = 3)+
#   theme(legend.position = "top")+
#   guides(color=guide_legend(title="E3 Ligase"))+
#   geom_jitter()+
#   geom_text_repel(box.padding = 0.5, max.overlaps = 500, label = pca_ErG_df2$ID)


############################################################
#                                                          #
#                          t-SNE                           #
#                                                          #
############################################################

library(Rtsne)

fortsne <- dataset_orig[!duplicated(dataset_orig[-c(1:6)]), ]

`%ni%` <- negate(`%in%`)

names(fortsne)

tsne_out <- Rtsne(as.matrix(fortsne[,-c(1:6)]))


# Conversion of matrix to dataframe
tsne_plot <- data.frame(Dim1 = tsne_out$Y[,1],
                        Dim2 = tsne_out$Y[,2],
                        ID = fortsne$id,
                        source = fortsne$source,
                        Class = fortsne$Class, 
                        stringsAsFactors = F)

# Plotting the plot using ggplot() function
ggplot2::ggplot(tsne_plot) + 
  geom_point(aes(x = Dim1, y = Dim2, color = Class), size = 3) +
  #geom_jitter()+
  ggtitle("t-SNE Plot for the E3 ligase binders collection") +
  theme_light(13)+
  theme(legend.position = "top")+
  guides(color=guide_legend(title="E3 Ligase"))

############################################################
#                                                          #
#                        tree part                         #
#                                                          #
############################################################

first_tree <- rpart(dataset_pca$Class ~ .,method = "class",data = dataset_pca[,-c(1)],
                    control =rpart.control(minsplit = 5,minbucket=1, cp=0.01))

rpart.plot(first_tree, extra = 104, # show fitted class, probs, percentages
           box.palette = "GnBu", # color scheme
           branch.lty = 3, # dotted branch lines
           shadow.col = "gray", # shadows under the node boxes
           nn = TRUE,
           main = "Tree for selectivity classification")


pred = predict(first_tree, type="class")
confusionMatrix(as.factor(pred), as.factor(dataset_orig$Class))

selection_tree <- rpart:::labels.rpart(first_tree)[-1]
selection_tree

selection_tree <- sub("\\>|\\>.*", "", selection_tree)
selection_tree

selection_tree <- unique(selection_tree)
selection_tree <- selection_tree[-5]

selection_RF <- c("Hf_Ac_d15", "Hf_Ac_d14", "Hf_Ac_d2", "Hf_D_d12", "Hf_Ac_d10", 
                  "Hf_D_d14", "Hf_D_d11", "D_Ac_d14", "D_Ac_d13")

length(unique(selection_tree))

newsel <- c("Class", "id", selection_tree)

toplot <- dataset_orig %>% dplyr::select(all_of(c(newsel))) %>% melt()

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

library(ggsignif)

ggplot(toplot, aes(x=Class, y=value, fill=Class)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
  geom_signif(comparisons =list(c("Active", "Inactive")), test='t.test', map_signif_level = T, vjust = 1.5)+
  facet_wrap(~variable, scales="free")+
  theme(legend.position = "none") +
  ggtitle("ErG bits distribution of most influential ErG bits for Ligase classification")

ggplot(toplot, aes(x=Class, y=value, fill=Class)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") +
  geom_signif(comparisons =list(c("CRBN", "IAP"), c("CRBN","Other"), c("CRBN", "VHL"), c("CRBN", "XIAP"),
                                c("IAP", "Other"), c("IAP", "VHL"), c("IAP", "XIAP"),
                                c("Other", "VHL"), c("Other", "XIAP"), c("VHL", "XIAP")), test='t.test', margin_top = 0.04,
              step_increase = 0.15, map_signif_level = T)+
  facet_wrap(~variable, scales="free")+
  theme_light(11) +
  theme(legend.position = "none") +
  ggtitle("ErG bits distribution of most influential ErG bits for Ligase classification")


############################################################
#                                                          #
#                      Random Forrest                      #
#                                                          #
############################################################

rf_df <- as.data.frame(dataset_pca, stringsAsFactors = F)

RF_task <- makeClassifTask(data=rf_df, target = "Class")

RF_learner<-makeLearner("classif.randomForest")

paraspace<-makeParamSet(makeIntegerParam("ntree", lower = 20, upper = 100), 
                        makeIntegerParam("mtry", lower = 5, upper = 50),
                        makeIntegerParam("nodesize", lower = 5, upper = 100), 
                        makeIntegerParam("maxnodes", lower = 10, upper = 50))

randsearch<-makeTuneControlRandom(maxit=50)

Kfold<-makeResampleDesc(method = "RepCV", folds=2, reps=5, stratify = T)

Turning_para<-tuneParams(RF_learner,RF_task,resampling = Kfold, 
                         par.set = paraspace, control = randsearch)

save(Turning_para, file = "RF_Hypertune_apr11_ErG.RData")

RF_Set <- setHyperPars(RF_learner, par.vals = Turning_para$x)

train_RF <- train(RF_Set, RF_task)

get_RF <- getLearnerModel(train_RF)

plot(get_RF, main = "RF optimization with ErG description")

labels <- colnames(get_RF$err.rate)

legend("topright", labels,
       col = 1:length(labels),
       lty = 1:length(labels))

outer <- makeResampleDesc("RepCV", reps = 2, stratify = T, folds=5)

paraspace2<-makeParamSet(makeIntegerParam("ntree", lower = 50, upper = 200), 
                         makeIntegerParam("mtry", lower = 5, upper = 10),
                         makeIntegerParam("nodesize", lower = 3, upper = 60), 
                         makeIntegerParam("maxnodes", lower = 45, upper = 100))

randsearch2<-makeTuneControlRandom(maxit=100)

forestWrapper <- makeTuneWrapper("classif.randomForest", resampling = outer,
                                 par.set = paraspace2,
                                 control = randsearch2)


cvWithTuning <- resample(forestWrapper, RF_task, resampling = outer)

save(cvWithTuning, file = "RF_selectivity_2_5_Apr11_v1_ErG.RData")

calculateConfusionMatrix(cvWithTuning$pred,relative = T)

cvWithTuning$pred$data$truth

cvWithTuning$pred$instance$group

############################################################
#                                                          #
#                  extract best RF model                   #
#                                                          #
############################################################

confusionMatrix(as.factor(cvWithTuning$pred$data$truth), as.factor(cvWithTuning$pred$data$response))

bestmodel <- randomForest(as.factor(Class) ~ ., data = rf_df,
                          ntree = 100, mtry = 9, 
                          nodesize =6, maxnodes = 82)

summary(bestmodel)

confusionMatrix(as.factor(bestmodel$predicted), as.factor(dataset_orig$Class))


# features importance plot 
# https://stats.stackexchange.com/questions/153663/how-to-get-the-most-important-variables-in-random-forests-in-r

feat_bestmodel <- importance(bestmodel) %>% 
  data.frame() %>% mutate(feature = row.names(.))

# make dataframe from importance() output

feat_imp_df <- importance(bestmodel) %>% 
  data.frame() %>% 
  mutate(feature = row.names(.)) %>% arrange(desc(MeanDecreaseGini)) %>% head(10)

# plot dataframe
ggplot(feat_imp_df, aes(x = reorder(feature, MeanDecreaseGini), 
                        y = MeanDecreaseGini)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_classic() +
  labs(
    x     = "Feature",
    y     = "Importance",
    title = "Top10 Descriptor Importance for RandomForest model"
  ) 

save(bestmodel, file = "bestmodel_ErG_RF_Apr11_Erg.RData")

############################################################
#                                                          #
#                         XGBoost                          #
#                                                          #
############################################################
library(xgboost)

#split into training (80%) and testing set (20%)
parts = createDataPartition(dataset_pca$Class, p = .8, list = F)
train = dataset_pca[parts,]
test = dataset_pca[-parts,]

#define predictor and response variables in training set
train_x = data.matrix(train[,-1])
train_y = as.numeric(unlist(as.factor(train$Class)))-1
test_x = data.matrix(test[,-1])
test_y = as.numeric(unlist(as.factor(test$Class)))-1

#define final training and testing sets
xgb_train = xgb.DMatrix(data = train_x, label = train_y)
xgb_test = xgb.DMatrix(data = test_x, label = test_y)

#define watchlist
watchlist = list(train=xgb_train, test=xgb_test)


# Create empty lists
lowest_error_list = list()
parameters_list = list()

# Create 10,000 rows with random hyperparameters
set.seed(1206)
for (iter in 1:100){
  param <- list(booster = "gbtree", objective = "multi:softmax", num_class = 5, 
                gamma=0,
                max_depth = sample(4:8, 1),
                eta = runif(1, .01, 1),
                subsample = runif(1, .7, 1),
                colsample_bytree = runif(1, .7, 1),
                min_child_weight = sample(1:7, 1)
  )
  parameters <- as.data.frame(param)
  parameters_list[[iter]] <- parameters
}

# Create object that contains all randomly created hyperparameters
parameters_df = do.call(rbind, parameters_list)

# Use randomly created parameters to create 10,000 XGBoost-models
set.seed(1206)
for (row in 1:nrow(parameters_df)){
  mdcv <- xgb.train(data=xgb_train,
                    booster = "gbtree", objective = "multi:softmax",
                    max_depth = parameters_df$max_depth[row],
                    eta = parameters_df$eta[row],
                    subsample = parameters_df$subsample[row],
                    colsample_bytree = parameters_df$colsample_bytree[row],
                    min_child_weight = parameters_df$min_child_weight[row],
                    nrounds= 200,
                    num_class = 6,
                    metric = "mlogloss",
                    early_stopping_rounds= 30,
                    print_every_n = 25,
                    watchlist = watchlist
  )
  lowest_error <- as.data.frame(min(mdcv$evaluation_log$test_mlogloss))
  lowest_error_list[[row]] <- lowest_error
}

# Create object that contains all accuracy's
lowest_error_df = do.call(rbind, lowest_error_list)

# Bind columns of accuracy values and random hyperparameter values
randomsearch = cbind(lowest_error_df, parameters_df)

# Quickly display highest accuracy
min(randomsearch$`min(mdcv$evaluation_log$test_mlogloss)`)

write_csv(randomsearch, "randomsearch_XGB_Apr11_200.csv")

# Prepare table
randomsearch <- as.data.frame(randomsearch) %>%
  mutate(val_acc = `min(mdcv$evaluation_log$test_mlogloss)`) %>% 
  arrange(val_acc)

randomsearch <- randomsearch[,-1]

# Tuned-XGBoost model
set.seed(1206)
params <- list(booster = "gbtree", objective = "multi:softmax",
               max_depth = randomsearch[1,]$max_depth,
               eta = randomsearch[1,]$eta,
               num_class=6,
               subsample = randomsearch[1,]$subsample,
               colsample_bytree = randomsearch[1,]$colsample_bytree,
               min_child_weight = randomsearch[1,]$min_child_weight)
xgb_tuned <- xgb.train(params = params,
                       data = xgb_train,
                       nrounds =800,
                       print_every_n = 10,
                       #eval_metric = "auc",
                       eval_metric = "mlogloss",
                       early_stopping_rounds = 30,
                       watchlist = list(train= xgb_train, test = xgb_test))

bestxgb <- xgboost(objective = "multi:softmax",
                   max_depth = 13, 
                   eta = 0.1430, 
                   subsample = 0.92, 
                   colsample_bytree = 0.83,
                   min_child_weight = 2,
                   data = xgb_train,
                   nrounds = 100,
                   num_class=6,
                   print_every_n = 10,
                   #eval_metric = "auc",
                   eval_metric = "mlogloss",
                   early_stopping_rounds = 50)#,
#watchlist = list(train= xgb_train, test = xgb_test))

# Make prediction on test
predicted <- predict(xgb_tuned, xgb_test)
predicted2 <- as.data.frame(predicted, stringsAsFactors=FALSE)

View(predicted2)
table(predicted2$predicted)

save(bestxgb, file = "bestxgb_E3selectivity_apr11.rda")

# Check accuracy with the confusion matrix
confusionMatrix(as.factor(test$Class), 
                factor(predicted2$predicted,
                       labels=c("CIAP1","CRBN","IAP", "Other", "VHL", "XIAP")),
                positive = "VHL", 
                dnn = c("Prediction", "Actual Data"))

importance_matrix = xgb.importance(colnames(xgb_train), model = xgb_tuned)

xgb.plot.importance(importance_matrix[1:10,], main = "Top 10 most important ErG bit", cex = 2)

selection_XGB <- c("Hf_D_d2","Hf_Ac_d2","Ac_Ac_d4", "D_Ac_d5","Hf_D_d4","D_Ac_d3", "Hf_Ac_d9", "Ar_Ac_d9", "Hf_Ac_d5", "Hf_Ac_d3" )

length(unique(selection_XGB))

newsel_XGB <- c("Class", "id", "source", selection_XGB)

toplot <- dataset_orig %>% dplyr::select(all_of(newsel_XGB)) %>% melt()

toplot$source

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

library(ggsignif)
library(ggforce)

ggplot(toplot, aes(x=source, y=value, fill=source)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
  geom_signif(comparisons =list(c("Literature", "Patent")), test='t.test', map_signif_level = T, vjust = 1.5)+
  facet_wrap(~variable, scales="free", ncol = 4)+
  theme(legend.position = "none") +
  ggtitle("ErG bits distribution of most influential ErG bits for Ligase classification")

ggplot(toplot, aes(x=Class, y=value, fill=Class)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") +
  geom_signif(comparisons =list(c("CIAP1", "CRBN"), c("CIAP1","Other"), c("CIAP1", "VHL"), c("CIAP1", "XIAP"),
                                c("CRBN", "IAP"), c("CRBN","Other"), c("CRBN", "VHL"), c("CRBN", "XIAP"),
                                c("IAP", "Other"), c("IAP", "VHL"), c("IAP", "XIAP"),
                                c("Other", "VHL"), c("Other", "XIAP"), c("VHL", "XIAP")), test='t.test', margin_top = 0.04,
              step_increase = 0.15, map_signif_level = T)+
  facet_wrap_paginate(~variable, scales="free", ncol = 2, nrow = 2, page = 3)+
  theme_light(14) +
  theme(legend.position = "none") +
  ggtitle("ErG bits distribution of most influential ErG bits for Ligase classification")

############################################################
#                                                          #
#                       Predictions                        #
#                                                          #
############################################################

Asinex <- read_excel("../data/2022_04_asinex_1257_molgluer_ErG.xlsx")

used <- names(train)[-1]

Asinex <- Asinex %>% select(all_of(c(used)))

xgb_Asinex = xgb.DMatrix(data = as.matrix(Asinex))

predicted_asinex <- predict(xgb_tuned, xgb_Asinex)

predicted_asinex_class <- factor(predicted_asinex, labels=c("CRBN","IAP", "Other", "VHL", "XIAP")) 

round(table(predicted_asinex_class)/1257,2)

# Broad library

Broad_ErG <- read_delim("../data/Broad_CIS_3D_ErG.txt", delim = ";")

Broad_ErG <- Broad_ErG %>% select(all_of(c(used)))

xgb_Broad = xgb.DMatrix(data = as.matrix(Broad_ErG))

predicted_Broad <- predict(xgb_tuned, xgb_Broad)

predicted_Broad_class <- factor(predicted_Broad, labels=c("CIAP1","CRBN","IAP", "Other", "VHL", "XIAP")) 

round(table(predicted_Broad_class)/9195,2)
