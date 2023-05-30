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

dataset_orig <- read_excel("e3_selectivity_FP.xlsx", sheet = "MACCS")

names(dataset_orig) <- make.names(names(dataset_orig))

############################################################
#                                                          #
#                          t-SNE                           #
#                                                          #
############################################################

library(Rtsne)

fortsne <- dataset_orig[!duplicated(dataset_orig[-c(1:6)]), ]

#fortsne_dupl <- dataset_orig[duplicated(dataset_orig[-c(1:4,6)]), ]

#dupl <- unlist(fortsne_dupl$id)

`%ni%` <- negate(`%in%`)

#fortsne <- fortsne %>% filter(id %ni% dupl)

names(fortsne)

tsne_out <- Rtsne(as.matrix(fortsne[,-c(1:6)]))

# Conversion of matrix to dataframe
tsne_plot <- data.frame(Dim1 = tsne_out$Y[,1],
                        Dim2 = tsne_out$Y[,2],
                        ID = fortsne$id,
                        source = fortsne$source,
                        Class = fortsne$Class, 
                        stringsAsFactors = F)

#tsne_plot2 <- tsne_plot %>% mutate(ID2 = ifelse(ID %in% tobechecked, ID, ""))


# Plotting the plot using ggplot() function
ggplot2::ggplot(tsne_plot) + 
  geom_point(aes(x = Dim1, y = Dim2, color = Class), size = 3) +
  #geom_jitter()+
  ggtitle("t-SNE Plot for the E3 ligase binders collection") +
  theme_light(13)+
  theme(legend.position = "top")+
  guides(color=guide_legend(title="E3 Ligase"))

ggplot2::ggplot(tsne_plot2) + geom_point(aes(x=Dim1,y=Dim2, color=Class), size = 3) +
  ggtitle("t-SNE Plot for the E3 ligase binders collection") +
  geom_jitter()+
  geom_text_repel(aes(x=Dim1,y=Dim2),max.overlaps = 50,  label = tsne_plot2$ID2)+
  theme(legend.position = "top")+
  guides(color=guide_legend(title="E3 Ligase"))#+
#geom_text_repel(aes(x=tsnel_plot$Dim1,y=tsnel_plot$Dim2), max.overlaps = 500, label = pca_ErG_df2$ID2, color = tsne_plot$Class)


############################################################
#                                                          #
#               Multidimensional Scaling MDS               #
#                                                          #
############################################################

library(ggpubr)

mds <- dataset_pca[,-c(1:4)] %>% dist() %>% cmdscale() %>% as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = dataset_pca$E3.Ligase,
          color = as.numeric(as.factor(dataset_pca$source)),
          size = 3,
          repel = TRUE)


############################################################
#                                                          #
#                    K-means clustering                    #
#                                                          #
############################################################

#optimal number of cluster

k.max <- 15
data <- tsne_plot[,c(1:2)]
sum(is.na(data))
wss <- sapply(1:k.max, function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
#wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# K=4

#another plot variation (https://arxiv.org/pdf/2212.12189.pdf)
mydata <- dataset_pca[,-c(1:4)]
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
wsse <- function(centers){
  sum(kmeans(mydata, centers=centers)$withinss)
}
#wsse <- c()
VDR <- c()

for (i in 2:12) {
  VDR[i] <- (wsse(1)-wsse(i)/i-1)/(wsse(i)/10-i)
}
plot(2:12, VDR[2:12], type = "line") # first knick = 3, second knick = 5


clust <- kmeans(mds, 5)$cluster %>% as.factor()
mds <- mds %>% mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = dataset_orig$E3.Ligase,
          color = "groups",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          title = "MDS for 685 E3_Ligase binders",
          ellipse.type = "convex",
          repel = TRUE)


############################################################
#                                                          #
#                           PLOT                           #
#                                                          #
############################################################

selection_XGB <- c("Ac_Ac_d4", "D_D_d3","Hf_Ar_d9","Hf_D_d7", "Hf_Ac_d4", "D_Ac_d12", "Ar_Ar_d2","Hf_D_d3","Hf_Ac_d2","Hf_Ac_d6" )

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

Asinex <- read_excel("../../2022_04_asinex_1257_molgluer_ErG.xlsx")

used <- names(train)[-1]

Asinex <- Asinex %>% select(all_of(c(used)))

xgb_Asinex = xgb.DMatrix(data = as.matrix(Asinex))

predicted_asinex <- predict(xgb_tuned, xgb_Asinex)

predicted_asinex_class <- factor(predicted_asinex, labels=c("CRBN","IAP", "Other", "VHL", "XIAP")) 

round(table(predicted_asinex_class)/1257,2)

# Broad library

Broad_ErG <- read_delim("../../Broad_CIS_3D_ErG.txt", delim = ";")

Broad_ErG <- Broad_ErG %>% select(all_of(c(used)))

xgb_Broad = xgb.DMatrix(data = as.matrix(Broad_ErG))

predicted_Broad <- predict(xgb_tuned, xgb_Broad)

predicted_Broad_class <- factor(predicted_Broad, labels=c("CIAP1","CRBN","IAP", "Other", "VHL", "XIAP")) 

round(table(predicted_Broad_class)/9195,2)

