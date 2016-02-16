#### TASKS:
# 1. create a one FULL example of running 2C - V
# 1.2 create a one FULL example of running 3C 
# 2. document each of the functions - V (still needs to be done for C2)
# 3. create a larger function to choose between the functions (C1, C2, C3)  - V

# 4. decide which functions to add/remove
# 5. decide on proper format to ask from our HBP friends (i.e.: how they should give us their files, and what output we will provide them with)
# 6. send to people (+bundle as a package)


# Consistent names in functions Consistent output for each function
# Speed - make it run faser. (for choosing num of clusters)

# add to CCC parameter for potential biomarkers



# Add marginal screening F test function

# Decide what to do with multinomial regression
# https://cran.r-project.org/web/packages/mnlogit/vignettes/mnlogit.pdf
# http://www.r-bloggers.com/how-to-multinomial-regression-models-in-r/
# https://cran.r-project.org/web/packages/mlogitBMA/mlogitBMA.pdf
# file:///C:/Users/user/Downloads/v34i12.pdf


RUNC3 <- FALSE
RUN <- FALSE
# RUN <- TRUE



if (RUN) {
  
  ############################################################# Generating the clinical measurements variables and the DX #
  
  # We'll generate 1200 observations (200 for each real diagnosis, 0
  # denotes Normal)
  n0 <- 200
  n1.1 <- 200
  n1.2 <- 200
  n1 <- n1.1 + n1.2
  n2 <- 200
  n3 <- 200
  n4 <- 200
  n <- n0 + n1.1 + n1.2 + n2 + n3 + n4
  real_DX <- c(rep(0, n0), rep(1.1, n1.1), rep(1.2, n1.2), rep(2, n2), 
               rep(3, n3), rep(23, n4))
  Doctors_DX <- c(rep(0, n0), rep(1, n1), rep(2, n2), rep(3, n3), rep(2, 
                                                                      n4/2), rep(3, n4/2))
  
  # Total of 36 CMs, we'll have 12 continuous, 12 discrete and 12 binary
  # variables. Out of 12 continuous 3 are meaningfull, Out of 12 discrete
  # 2 are meaningfull, Out of 12 binary 2 are meaningfull
  CM <- matrix(NA, nrow = n, ncol = 36)
  parameters <- matrix(c(10, 14, 17, 40, 60, 40, 15, 18, 23, 50, 25, 
                         50, 35, 7, 14, 10, 30, 18, rep(0, 9 * 6)  #none influencing continuous vars       
                         , 
                         1, 20, 27, 13, 5, 5, 29, 9, 5, 18, 24, 24, 4, 4, 4, 4, 4, 4  #none influencing discrete var
                         , rep(0, 9 * 6)  #none influencing continuous vars
                         , 0.1, 0.2, 0.05, 0.95, 
                         0.975, 0.975, 0.089, 0.9, 0.93, 0.05, 0.1, 0.1, 0.5, 0.5, 0.5, 
                         0.5, 0.5, 0.5  #none influencing binary var
                         , 
                         rep(0, 9 * 6)), nrow = 6)  #none influencing binary vars
  
  
  
  
  # Forth Scenario: Noise from influencing vars + Noise from none
  # influencing
  
  min_sd <- 0.1
  max_sd <- 3
  s.d <- matrix(runif(12 * 6, min = min_sd, max = max_sd), nrow = 6)
  noise_dis <- 2.55
  
  # Uncorrelated continous variables
  mu_uncorr <- c(4, 20, 55, 2, 19, 37)
  s.d_uncorr <- c(1, 5, 13, 0.5, 12, 22)
  # Uncorrelated discrete variables
  dis_uncor <- c(40, 12, 3, 77, 10, 14)
  # probabilities for uncorrelated binary variables
  bin_uncor <- c(0.2, 0.4, 0.1, 0.7, 0.5, 0.55)
  # rho for correlated variables
  rho <- c(0.2, 0.4, 0.7)
  
  # create correlated variables function
  getBiCop <- function(n, rho, mar.fun = rnorm, x = NULL, ...) {
    if (!is.null(x)) {
      X1 <- x
    } else {
      X1 <- mar.fun(n, ...)
    }
    if (!is.null(x) & length(x) != n) 
      warning("Variable x does not have the same length as n!")
    
    C <- matrix(rho, nrow = 2, ncol = 2)
    diag(C) <- 1
    
    C <- chol(C)
    
    X2 <- mar.fun(n)
    X <- cbind(X1, X2)
    
    # induce correlation (does not change X1)
    df <- X %*% C
    
    return(df)
  }
  
  
  # for each variable we assign the relevant values according to the real
  # diagnosis
  for (j in 1:3) {
    CM[, j] <- c(rnorm(n0, parameters[1, j], s.d[1, j]), rnorm(n1.1, 
                                                               parameters[2, j], s.d[2, j]), rnorm(n1.2, parameters[3, j], 
                                                                                                   s.d[3, j]), rnorm(n2, parameters[4, j], s.d[4, j]), rnorm(n3, 
                                                                                                                                                             parameters[5, j], s.d[5, j]), rnorm(n4, parameters[6, j], s.d[6, 
                                                                                                                                                                                                                           j]))
    set.seed(80)
    CM[, j + 3] <- getBiCop(n, rho[j], x = CM[, 1])[, 2]  #3 correlated  
    CM[, j + 6] <- c(rnorm(n, mu_uncorr[j], s.d_uncorr[j]))
    CM[, j + 9] <- c(rnorm(n, mu_uncorr[j + 3], s.d_uncorr[j + 3]))
    
    CM[, j + 12] <- c(parameters[1, j + 12] + round(runif(n0, min = 0, 
                                                          max = noise_dis), 0), parameters[2, j + 12] + round(runif(n1.1, 
                                                                                                                    min = 0, max = noise_dis), 0), parameters[3, j + 12] + round(runif(n1.2, 
                                                                                                                                                                                       min = 0, max = noise_dis), 0), parameters[4, j + 12] + round(runif(n2, 
                                                                                                                                                                                                                                                          min = 0, max = noise_dis), 0), parameters[5, j + 12] + round(runif(n3, 
                                                                                                                                                                                                                                                                                                                             min = 0, max = noise_dis), 0), parameters[6, j + 12] + round(runif(n4, 
                                                                                                                                                                                                                                                                                                                                                                                                min = 0, max = noise_dis), 0))
    CM[, j + 15] <- round(getBiCop(n, rho[j], x = CM[, j + 12])[, 2], 
                          0)  #3 correlated to discrete
    CM[, j + 18] <- c(dis_uncor[j] + round(runif(n, min = 0, max = noise_dis), 
                                           0))  #6 non correlated
    CM[, j + 21] <- c(dis_uncor[j + 3] + round(runif(n, min = 0, max = noise_dis), 
                                               0))
    
    CM[, j + 24] <- c(rbinom(n0, 1, parameters[1, j + 24]), rbinom(n1.1, 
                                                                   1, parameters[2, j + 24]), rbinom(n1.2, 1, parameters[3, j + 
                                                                                                                           24]), rbinom(n2, 1, parameters[4, j + 24]), rbinom(n3, 1, parameters[5, 
                                                                                                                                                                                                j + 24]), rbinom(n4, 1, parameters[6, j + 24]))
    CM[, j + 27] <- rbinom(n, 1, ifelse(CM[, j + 24] == 0, 0.2, 0.8))  #3 correlated            
    CM[, j + 30] <- c(rbinom(n, 1, bin_uncor[j]))  #6 uncorrelated
    CM[, j + 33] <- c(rbinom(n, 1, bin_uncor[j + 3]))
  }
  
  # Looking at the DB
  library(MASS)
  CM_true <- data.frame(CM[, 1:3], CM[, 13:14], CM[, 25:26])
  parcoord(CM_true, col = 1, var.label = TRUE, main = "Influencing vars")
  parcoord(CM[, 1:12], col = 1, var.label = TRUE, main = "Normal vars")
  parcoord(CM[, 13:24], col = 1, var.label = TRUE, main = "Discrete vars")
  parcoord(CM[, 25:36], col = 1, var.label = TRUE, main = "binary vars")
  CM_DX <- data.frame(CM, real_DX, Doctors_DX)
  colnames(CM_DX)[37:38] <- c("real_DX", "Doctors_DX")
  
  # Perform transformations- will be needed in advanced settings
  CM_DX_trans <- CM_DX
  
  ## Simulated data ## View(CM_DX_trans)
  x <- CM_DX_trans[, 1:36]
  y <- CM_DX_trans[, 38]
  
  
}

# creating the PB observed matrix

if (FALSE) {
  
  # first we create a matrix for the diseases (800 patients), only
  # meaningfull PBs Each disease is defined by 3 PBs (each patient has
  # exactly one of them)
  PB_sick <- matrix(NA, nrow = sum(n_vec_real[2:6]), ncol = real_PB)
  begin <- 0
  offset <- 0
  for (i in 1:(length(n_vec_real) - 1)) {
    begin <- 1 + offset
    offset <- offset + num_sick[1]
    # normally distributed PB, will be meaningfull for 35% of people with
    # true DX
    PB_sick[begin:offset, i] <- rnorm(length(begin:offset), mean = marker[i], 
                                      s.d_PB)
    PB_sick[is.na(PB_sick[, i]), i] <- rnorm(length(PB_sick[is.na(PB_sick[, 
                                                                          i]), i]), mean = no_PB[i], sd = s.d_PB)
    
    begin <- begin + offset
    offset <- offset + num_sick[2]
    PB_sick[begin:offset, i + (length(n_vec_real) - 1)] <- rnorm(length(begin:offset), 
                                                                 mean = marker[i + (length(n_vec_real) - 1)], s.d_PB)
    PB_sick[is.na(PB_sick[, i + (length(n_vec_real) - 1)]), i + (length(n_vec_real) - 
                                                                   1)] <- rnorm(length(PB_sick[is.na(PB_sick[, i + (length(n_vec_real) - 
                                                                                                                      1)]), i + (length(n_vec_real) - 1)]), mean = no_PB[i + (length(n_vec_real) - 
                                                                                                                                                                                1)], sd = s.d_PB)
    
    begin <- begin + offset
    offset <- offset + num_sick[3]
    # Binary PB, will be meaningfull for 30% of people with true DX
    PB_sick[begin:offset, i + (length(n_vec_real) - 1) * 2] <- marker[i + 
                                                                        (length(n_vec_real) - 1) * 2]
    PB_sick[is.na(PB_sick[, i + (length(n_vec_real) - 1) * 2]), i + 
              (length(n_vec_real) - 1) * 2] <- no_PB[i + (length(n_vec_real) - 
                                                            1) * 2]
    
    offset <- n_vec_real[i + 1]
  }
  
  # Now we create the healthy patients matrix (real DX=0). Only
  # meaningfull PBs For each PB a small portion of the healthy patients
  # will have the PB
  offset <- 0
  begin <- 0
  PB0 <- matrix(NA, nrow = n_vec_real[1], ncol = real_PB)
  for (i in 1:5) {
    begin <- offset + 1
    offset <- offset + num_health[1]
    PB0[begin:offset, i] <- rnorm(length(begin:offset), mean = marker[i], 
                                  sd = s.d_PB)
    PB0[is.na(PB0[, i]), i] <- rnorm(length(PB0[is.na(PB0[, i]), i]), 
                                     mean = no_PB[i], sd = s.d_PB)
    
    begin <- begin + num_health[1]
    offset <- offset + num_health[2]
    PB0[begin:offset, i + (length(n_vec_real) - 1)] <- rnorm(length(begin:offset), 
                                                             marker[i + (length(n_vec_real) - 1)], sd = s.d_PB)
    PB0[is.na(PB0[, i + (length(n_vec_real) - 1)]), i + (length(n_vec_real) - 
                                                           1)] <- rnorm(length(PB0[is.na(PB0[, i + (length(n_vec_real) - 
                                                                                                      1)]), i + (length(n_vec_real) - 1)]), no_PB[i + (length(n_vec_real) - 
                                                                                                                                                         1)], sd = s.d_PB)
    
    begin <- begin + num_health[2]
    offset <- offset + num_health[3]
    PB0[begin:offset, i + (length(n_vec_real) - 1) * 2] <- marker[i + 
                                                                    (length(n_vec_real) - 1) * 2]
    PB0[is.na(PB0[, i + (length(n_vec_real) - 1) * 2]), i + (length(n_vec_real) - 
                                                               1) * 2] <- no_PB[i + (length(n_vec_real) - 1) * 2]
    
  }
  PB <- matrix(data = rbind(PB0, PB_sick), ncol = real_PB)
  
  # Now we'll add the none influancing variables using the given
  # correlation. We want a total of 100 PBs
  PB <- matrix(data = cbind(PB, matrix(data = NA, ncol = (100 - col(PB)), 
                                       nrow = n)), ncol = 100, nrow = n)
  for (i in (real_PB + 1):dim(PB)[2]) {
    PB[, i] <- getBiCop(n, rho_PB, X1 = PB[, (i - length(real_PB))])[, 
                                                                     2]
  }
  
  # turning influancing variables 6-10 to categorial and influancing
  # variables 11-15 to binary
  begin <- 0
  offset <- 0
  for (i in 1:length(n_vec_real)) {
    begin <- offset + 1
    offset <- offset + n_vec_real[i]
    PB[begin:offset, 6] <- as.numeric(as.character(cut(PB[begin:offset, 
                                                          6], breaks = c(min(PB[begin:offset, 6]) - 1, (mean(PB[begin:offset, 
                                                                                                                6]) - sd(PB[begin:offset, 6])), (mean(PB[begin:offset, 6]) + 
                                                                                                                                                   sd(CM[begin:offset, 6])), Inf), labels = c(round(mean(PB[begin:offset, 
                                                                                                                                                                                                            6]), 0) - 1, round(mean(PB[begin:offset, 6]), 0), round(mean(PB[begin:offset, 
                                                                                                                                                                                                                                                                            6]), 0) + 1))))
    PB[begin:offset, 7] <- as.numeric(as.character(cut(PB[begin:offset, 
                                                          7], breaks = c(min(PB[begin:offset, 7]) - 1, (mean(PB[begin:offset, 
                                                                                                                7]) - sd(PB[begin:offset, 7])), (mean(PB[begin:offset, 7]) + 
                                                                                                                                                   sd(CM[begin:offset, 7])), Inf), labels = c(round(mean(PB[begin:offset, 
                                                                                                                                                                                                            7]), 0) - 1, round(mean(PB[begin:offset, 7]), 0), round(mean(PB[begin:offset, 
                                                                                                                                                                                                                                                                            7]), 0) + 1))))
    PB[begin:offset, 8] <- as.numeric(as.character(cut(PB[begin:offset, 
                                                          8], breaks = c(min(PB[begin:offset, 8]) - 1, (mean(PB[begin:offset, 
                                                                                                                8]) - sd(PB[begin:offset, 8])), (mean(PB[begin:offset, 8]) + 
                                                                                                                                                   sd(CM[begin:offset, 8])), Inf), labels = c(round(mean(PB[begin:offset, 
                                                                                                                                                                                                            8]), 0) - 1, round(mean(PB[begin:offset, 8]), 0), round(mean(PB[begin:offset, 
                                                                                                                                                                                                                                                                            8]), 0) + 1))))
    PB[begin:offset, 9] <- as.numeric(as.character(cut(PB[begin:offset, 
                                                          9], breaks = c(min(PB[begin:offset, 9]) - 1, (mean(PB[begin:offset, 
                                                                                                                9]) - sd(PB[begin:offset, 9])), (mean(PB[begin:offset, 9]) + 
                                                                                                                                                   sd(CM[begin:offset, 9])), Inf), labels = c(round(mean(PB[begin:offset, 
                                                                                                                                                                                                            9]), 0) - 1, round(mean(PB[begin:offset, 9]), 0), round(mean(PB[begin:offset, 
                                                                                                                                                                                                                                                                            9]), 0) + 1))))
    PB[begin:offset, 10] <- as.numeric(as.character(cut(PB[begin:offset, 
                                                           10], breaks = c(min(PB[begin:offset, 10]) - 1, (mean(PB[begin:offset, 
                                                                                                                   10]) - sd(PB[begin:offset, 10])), (mean(PB[begin:offset, 10]) + 
                                                                                                                                                        sd(CM[begin:offset, 10])), Inf), labels = c(round(mean(PB[begin:offset, 
                                                                                                                                                                                                                  10]), 0) - 1, round(mean(PB[begin:offset, 10]), 0), round(mean(PB[begin:offset, 
                                                                                                                                                                                                                                                                                    10]), 0) + 1))))
    PB[begin:offset, 11] <- ifelse(PB[begin:offset, 11] <= quantile(PB[begin:offset, 
                                                                       11], probs = 0.95), round(mean(PB[begin:offset, 11])), 1 - 
                                     round(PB[begin:offset, 11]))
    PB[begin:offset, 12] <- ifelse(PB[begin:offset, 12] <= quantile(PB[begin:offset, 
                                                                       12], probs = 0.95), round(mean(PB[begin:offset, 12])), 1 - 
                                     round(PB[begin:offset, 12]))
    PB[begin:offset, 13] <- ifelse(PB[begin:offset, 13] <= quantile(PB[begin:offset, 
                                                                       13], probs = 0.95), round(mean(PB[begin:offset, 13])), 1 - 
                                     round(PB[begin:offset, 13]))
    PB[begin:offset, 14] <- ifelse(PB[begin:offset, 14] <= quantile(PB[begin:offset, 
                                                                       14], probs = 0.95), round(mean(PB[begin:offset, 14])), 1 - 
                                     round(PB[begin:offset, 14]))
    PB[begin:offset, 15] <- ifelse(PB[begin:offset, 15] <= quantile(PB[begin:offset, 
                                                                       15], probs = 0.95), round(mean(PB[begin:offset, 15])), 1 - 
                                     round(PB[begin:offset, 15]))
    
  }
  
  # Transform some of the correlated variables to Discrete
  begin <- 0
  offset <- 0
  # going disease by disease (according to the number of observations in
  # each disease)
  for (j in 1:length(n_vec_real)) {
    begin <- offset + 1
    offset <- offset + n_vec_real[j]
    for (i in (2 * real_PB):(3 * real_PB)) # changing columns 30 to 45
    {
      PB[begin:offset, i] <- as.numeric(as.character(cut(PB[begin:offset, 
                                                            i], breaks = c(min(PB[begin:offset, i]) - 1, (mean(PB[begin:offset, 
                                                                                                                  i]) - sd(PB[begin:offset, i])), (mean(PB[begin:offset, 
                                                                                                                                                           i]) + sd(PB[begin:offset, i])), Inf), labels = c(round(mean(PB[begin:offset, 
                                                                                                                                                                                                                          i]), 0) - 1, round(mean(PB[begin:offset, i]), 0), round(mean(PB[begin:offset, 
                                                                                                                                                                                                                                                                                          i]), 0) + 1))))
    }
  }
  colnames(PB) <- paste(rep("PB", dim(PB)[2]), c(1:dim(PB)[2]))
  
  
  full_DB <- data.frame(real_DX_f, assigned_DX_f, CM, PB)
  return(full_DB)
}


visualization <- function(db, real) {
  influence_features <- db[, 1:real]
  none_influence <- db[, (real + 1):dim(db)[2]]
  parcoord(influence_features, col = c(1:6), var.label = TRUE, main = "Influencing vars")
  parcoord(none_influence[, 1:floor(dim(none_influence)[2]/2)], col = c(1:6), 
           var.label = TRUE, main = "None influencing vars, part I")
  parcoord(none_influence[, ceiling(dim(none_influence)[2]/2):dim(none_influence)[2]], 
           col = c(1:6), var.label = TRUE, main = "None influencing vars, part II")
}







if (FALSE) {
  
  # checks for each function - no need to run.
  
  
  
  feature_selection(x, y, feature_selection_method = "AIC_MSFDR")
  
  Feature_Selection_RF(x, y, Importance = T, a = 10, b = 14)  # Checking functionality
  
  tmp <- Feature_Selection_RF(x, y, Importance = T, a = 10, b = 14)  # Checking functionality
  important_var_RF_large <- tmp[[3]]
  
  Feature_Selection_BIC(x, y, nvmax = 5, nbest = 12, plot = TRUE)
  
  Feature_Selection_AIC_MSFDR(x, y, q = 0.5, summary.MSFDR = FALSE, print.the.steps = T)
  
  Feature_Selection_AIC(x, y, print.summary.AIC = TRUE)
  
  vars_AICMSFDR <- c(1:4, 7, 10, 13:19, 22, 25:36)
  Reduced_X <- x[, vars_AICMSFDR]
  
  # scaling it
  Scaled_Reduced_CM_trans <- scale(Reduced_CM_trans)
  
  clusters_medoids(Scaled_Reduced_CM_trans)
  
  clusters_manhattan(Scaled_Reduced_CM_trans, 5)
  
  
  clusGap_best <- clusGap(Scaled_Reduced_CM_trans, FUN = hclust_k_euc, 
                          K.max = 10, B = 60, verbose = FALSE)
  plot(clusGap_best, main = "Gap statistic, hclust Euclidean")
  
  clusGap_best <- clusGap(Scaled_Reduced_CM_trans, FUN = hclust_k_man, 
                          K.max = 10, B = 60, verbose = FALSE)
  plot(clusGap_best, main = "Gap statistic, hclust Manhattan")
  
  Gap_euclidean(Reduced_CM_trans, 5)
  Gap_manhattan(Reduced_CM_trans, 5)
  
  k.gap <- 4
  #### k.gap<-3 k.gap<-7 k.gap<-5
  k.true <- 6
  factor_DX_real <- as.factor(CM_DX_trans[, 37])
  
  cluster_euclidean(Scaled_Reduced_CM_trans, k.gap)
  
  cluster_manhattan(Scaled_Reduced_CM_trans)
  
  cluster_hierarchical_euc(Scaled_Reduced_CM_trans)
  
  cluster_hierarchical_man(Scaled_Reduced_CM_trans)
  
  
  
  
}





if (RUN) {
  
  # TODO: give me one FULL example of running 3C
  
  # xy <- get_xy_from_DATA_C2(DATA, META_DATA) x <- xy$x y <- xy$y
  
  
  CM_imp_var <- feature_selection(x, y, method = "AIC_MSFDR", summary.MSFDR = TRUE)
  CM_imp_var
  CM_final_vars <- CM_imp_var[[1]][2]  # Extracting a list of important CM variables
  CM_final_vars
  newX <- x[, unlist(CM_final_vars)]  # Creating a reduced new X matrix accordingly
  num_clust <- number_of_clusters(newX, method = "Manhattan", B = 90, 
                                  k.max = 10, plot = TRUE)
  num_clust
  k.gap <- 6  # From looking at the number of clusters function's plot
  final_cluster <- clustering(newX, y, k.gap, plot = TRUE)
  final_cluster
  
  dim(x)
  # k=7,
  C2(x, y, feature_selection_method = "AIC_MSFDR", num_clusters_method = "Manhattan", 
     clustering_method = " hclust_man")
  
}

if (RUNC3) {
  PB_imp_var <- feature_selection(x = CM_DX_comp[, -194], y = CM_DX_comp[, 
                                                                         194], method = "RF")
  PB_final_vars <- PB_imp_var[[1]][1]
  PB_final_vars
  newPB <- data.frame(y, PB[, unlist(PB_final_vars)])
  classification.PB <- C3(newPB, method = "RF", downSampling = TRUE)
}



create_data_metadata <- FALSE

# install.packages('ADNIMERGE', )
# http://stackoverflow.com/questions/1474081/how-do-i-install-an-r-package-from-source


if (create_data_metadata) {
  if (FALSE) {
    install.packages("inst\\packages\\ADNIMERGE_0.0.1.tar.gz", repos = NULL, 
                     type = "source")
    install.packages("reshape")
    install.packages("mice")
  }
  
  library(ADNIMERGE)
  library(plyr)
  library(reshape)
  library(psych)
  # combine the assigned DX vector with the CM (after pre-processing)
  # matrix
  DX <- data.frame(adnimerge[adnimerge$VISCODE == "bl", "DX.bl"], adnimerge[adnimerge$VISCODE == 
                                                                              "bl", "RID"])
  names(DX)[1] <- "DX"
  names(DX)[2] <- "RID"
  CM <- read.csv("C:/Users/tal/Dropbox/ADNI 3-C/Pre-processing/OUTPUT trans_Db.csv")
  # CM <- read.csv('inst/data/OUTPUT trans_Db.csv')
  
  CM <- CM[, 2:dim(CM)[[2]]]
  CM_DX <- merge(CM, DX, by = "RID", all.x = TRUE)
  subset_rows <- complete.cases(CM_DX[, "DX"])
  CM_DX_comp <- CM_DX[subset_rows, ]
  
  # Variables 'PTACOGBEG' and 'PTAADDX' has non missing values only for
  # diagnosed patients so, we'll exclude them from the feature selection
  # step:
  CM_DX_comp <- subset(CM_DX_comp, select = -c(PTACOGBEG.None, PTAADDX.power_2))
  # We'll check for highly correlated variables and exclude one of each
  # pair
  CM_comp <- subset(CM_DX_comp, select = -c(RID, DX))
  cor_mat <- as.matrix(cor(CM_comp))
  cor_mat_melt <- arrange(melt(cor_mat), -abs(value))
  t <- subset(cor_mat_melt, value > 0.9)
  t <- subset(t, X1 != X2)
  # pairs.panels(CM_comp[, t[, 'X1']])
  CM_DX_comp <- subset(CM_DX_comp, select = -c(MMSE.logit_frac, delayed_recall_mmse.logit_frac, 
                                               ADAS11.logit, EcogSPMem.ordered, CDRSB.ordered))
  RID <- CM_DX_comp[, "RID"]
  
  library(mice)
  
  # x <- CM_DX_comp[, -c(1, 189)] y <- CM_DX_comp$DX
  
  
  # # PB matrix
  PB_imputed <- read.csv("C:/Users/Tal/Dropbox/ADNI 3-C/07-04-2015/results/assigned_clusters_with_PB_after_imputation.csv")
  real.diagnosis <- read.csv("C:\\Users\\Tal\\Dropbox\\ADNI 3-C\\Tal\\25.10.15\\adni_DX_bl_v1.csv")
  # PB_imputed <-
  # read.csv('inst/data/assigned_clusters_with_PB_after_imputation.csv')
  # real.diagnosis <- read.csv('inst/data/adni_DX_bl_v1.csv')
  
  real.diagnosis <- real.diagnosis[, -1]
  # attaching
  PB_DX <- merge(real.diagnosis, PB_imputed, by = "RID", all.y = TRUE)
  length(real.diagnosis$RID)
  length(PB_imputed$RID)
  # PB_DX <- PB_DX[,-1] PBx <- PB_DX[, -c(1,2)] PBy <- PB_DX[, 2]
  
  # META_DATA A meta-data file including two columns. One named varName
  # which is the variable/column name from the data file, and another
  # column named varCategory which is the categorization (i.e.: Dx, CM or
  # PB).  DATA A data file with Dx (disease diagnosis), CM (Clinical
  # measurements) PB (Potential biomarkers)
  
  DATA <- merge(CM_DX_comp, PB_DX, by = "RID")
  DATA$DX.y <- NULL
  DATA$DX <- DATA$DX.x
  DATA$DX.x <- NULL
  
  # create META_DATA object
  META_DATA <- data.frame(varName = colnames(DATA), varCategory = NA, 
                          stringsAsFactors = FALSE)
  ss <- META_DATA$varName %in% colnames(CM_DX_comp)
  META_DATA[ss, "varCategory"] <- "CM"
  ss <- META_DATA$varName %in% colnames(PB_DX)
  META_DATA[ss, "varCategory"] <- "PB"
  ss <- META_DATA$varName %in% "RID"
  META_DATA[ss, "varCategory"] <- "RID"
  ss <- META_DATA$varName %in% "DX"
  META_DATA[ss, "varCategory"] <- "DX"
  
  head(META_DATA)
  table(META_DATA[, 2])
  
}



if (FALSE) {
  library(CCC)
  tmp <- get_xy_from_DATA_C2(DATA, META_DATA)
  str(tmp)
  x <- tmp[[1]]  # CM matrix
  y <- tmp[[2]]  # DD vector
  PBx <- get_PBx_from_DATA_C3(DATA, META_DATA)  # PB matrix
  feature_selection_methods <- c("RF", "AIC_MSFDR", "AIC", "BIC")
  num_clusters_methods <- c("Euclidean", "Manhattan", "hclust_Euclidean", 
                            "hclust_Manhattan")
  clustering_methods <- c("Euclidean", "Manhattan", "Heuclidean", "Hmanhattan")
  C2_combinations <- expand.grid(feature_selection_methods, num_clusters_methods, 
                                 clustering_methods)  # All possible combinations for C2
  classification_methods <- c("RF", "RF_downsampling", "CART_information", 
                              "CART_gini")
  C3_combinations <- expand.grid(feature_selection_methods, classification_methods)  # All possible combinations for C3
  results_C2 <- list(NA)
  newy <- c(NA)
  results_C3 <- list(NA)
  results <- list(results_C2 = list(results_C3 = NA))
  for (i in 1:dim(C2_combinations)[1]) {
    if (C2_combinations[i, 1] == "BIC") {
      alter_x <- x[, 1:50]
      results_C2[[i]] <- C2(alter_x, y, C2_combinations[i, 1], C2_combinations[i, 
                                                                               2], C2_combinations[i, 3], plot.num.clust = FALSE, plot.clustering = FALSE, 
                            k = 6)
      newy <- unlist(results_C2[[i]][[3]])
    } else {
      results_C2[[i]] <- C2(x, y, C2_combinations[i, 1], C2_combinations[i, 
                                                                         2], C2_combinations[i, 3], plot.num.clust = FALSE, plot.clustering = FALSE, 
                            k = 6)
      newy <- unlist(results_C2[[i]][[3]])
    }
    for (j in 1:dim(C3_combinations)[1]) {
      print(C2_combinations[i, c(1:3)])
      print(C3_combinations[j, c(1, 2)])
      if (C3_combinations[j, 1] == "BIC") {
        alter_PBx <- PBx[, 1:50]
        results_C3[[i]] <- C3(alter_PBx, newy, C3_combinations[j, 
                                                               1], C3_combinations[j, 2])
      } else {
        
        results_C3[[i]] <- C3(PBx, newy, C3_combinations[j, 1], 
                              C3_combinations[j, 2])
      }
    }
    
    results[[i]] <- results_C2[[i]]
    results[[i]][[j]] <- results_C3[[j]]
  }
  
  
}