## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,eval = F, # change to eval = T to run the code and get the output, it takes a while.
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE, tidy=TRUE----------------------------
#  DATA <- read.csv("../data/DATA.csv")
#  META_DATA <- read.csv("../data/META_DATA.csv")[1:200, -1]
#  head(META_DATA)
#  table(META_DATA[,"varCategory"])
#  head(DATA)[, 1:5]

## ---- message=FALSE, warning=FALSE, tidy=TRUE----------------------------
#  library(CCC)
#  CM_AD_matrix <- get_xy_from_DATA_C2(DATA, META_DATA)
#    x <- CM_AD_matrix[[1]]  # CM matrix
#    y <- CM_AD_matrix[[2]]  # DD vector
#  
#    C2_results <- C2(x, y, feature_selection_method="RF", num_clusters_method="Manhattan", clustering_method="Manhattan", plot.num.clus=TRUE, plot.clustering=TRUE, k=6)
#  
#  C2_results

## ---- message=FALSE, warning=FALSE, tidy=TRUE----------------------------
#  # Creating PB matrix
#  PBx <- get_PBx_from_DATA_C3(DATA, META_DATA)
#  newy <- unlist(C2_results[[3]])
#  C3_results <- C3(PBx, newy, feature_selection_method="RF", classification_method="RF")
#  C3_results
#  table(newy, C3_results[[2]])
#  error.rate <- sum(newy != C3_results[[2]]) / length(newy)
#  error.rate

