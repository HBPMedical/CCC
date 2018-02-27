## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE, tidy=TRUE----------------------------
data(c3_sample2)
data(c3_sample2_categories)
head(c3_sample2_categories) 
table(c3_sample2_categories[,"varCategory"])
head(c3_sample2)[, 1:5]

## ---- message=FALSE, warning=FALSE, tidy=TRUE----------------------------
library(CCC)  
CM_AD_matrix <- get_xy_from_DATA_C2(c3_sample2, c3_sample2_categories)
  x <- CM_AD_matrix[[1]]  # CM matrix
  y <- CM_AD_matrix[[2]]  # DD vector
  
  C2_results <- C2(x[,1:50] , y, feature_selection_method="RF", num_clusters_method="Manhattan", clustering_method="Manhattan", plot.num.clus=TRUE, plot.clustering=TRUE, k=6)

C2_results

## ---- message=FALSE, warning=FALSE, tidy=TRUE----------------------------
# Creating PB matrix
PBx <- get_PBx_from_DATA_C3(c3_sample2, c3_sample2_categories)
newy <- unlist(C2_results[[3]])
C3_results <- C3(PBx, newy, feature_selection_method="RF", classification_method="RF") 
C3_results
table(newy, C3_results[[2]])

