
CCC - Categorization, Clustering and Classification
===================================================

The package implements the 3C-strategy\\pipeline for the refinement of disease diagnoses in medical research. The first step in the analysis pipeline is manual **C**ategorization of the feature set to (i) current clinical measures, (ii) potential biomarkers and (iii) assigned diagnosis.

In the beginning of the second step (**C**lustering), a subset of the clinical measures is selected via supervised algorithm (Random Forest, LASSO, etc) where the target variable is the assigned diagnosis. Then, the selected measures are used to determine the number of clusters for the clustering and the clustering algorithm is applied (K-means, K-medoids, and Hierarchical clustering).

In step 3 a new model is trained, using the potential biomarkers as features, to **C**lassify the data according to the cluster labels created in step 2.

Installation
------------

You can install CCC from github with:

``` r
# install.packages("devtools")
devtools::install_github("HBPMedical/CCC")
```

Example
-------

This is a basic example:

``` r
## basic example code...
```
