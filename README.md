# Bayesian dimensionality reduction with PCA using penalized semi-integrated likelihood

## Overview

In this repository you can find a code that can be used to recreate figures from the paper **Bayesian dimensionality reduction with PCA using penalized semi-integrated likelihood**.

## Important informations

1. Running this code requires both R and MATLAB. In files **figure_2.R**, **figure_3.R**, 
**data_generation_robustness_fixed_effects.R** one needs to specify path to MATLAB executable file.
2. Please note that it might take a lot of time to run simulations (several days). We recomend to separate performing simulations and plotting results. To speeed up calculations one might reduce number of repetitions (*numb.repetitions*), reduce the set of number of variables included in simulations (*vars*) or reduce the set of tested signal to noise ratios (*SNRs*).
3. Packages required for this simulations:
   * [varclust](https://github.com/psobczyk/varclust)
   * FactoMineR
   * pryr
   * softImpute
   * ggplot2 (for plots)

To install all those dependencies just run the code below in R console
```
library(devtools)
install_github("psobczyk/varclust")
install.packages(c("pryr", "softImpute", "FactoMineR", "ggplot2"))
```

TO DO:

* remove pesel dependency for figure 1
* split creating figure into: a) data generation b) results plotting
* update path to temp files for figures 2,3,7
