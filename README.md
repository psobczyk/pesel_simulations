# Bayesian dimensionality reduction with PCA using penalized semi-integrated likelihood

## Overview

In this repository you can find a code used to recreate figures from the paper **Bayesian dimensionality reduction with PCA using penalized semi-integrated likelihood**.

## Important informations

1. Running this code requires both R and MATLAB. In files **figure_2.R**, **figure_3.R**, **figure_4_5_6.R** you need to specify path to MATLAB executable file.
2. Please note that it might take a lot of time to run simulations (several days). We recomend to separate performing simulations and plotting results. To speeed up calculations one might reduce number of repetitions (*numb.repetitions*), reduce the set of number of variables included in simulations (*vars*) or reduce the set of tested signal to noise rations (*SNRs*).
