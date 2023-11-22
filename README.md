## Introduction
LongDat R package takes longitudinal dataset as input data and analyzes if there is significant change of the features over time (proxy for treatments), while detects and controls for covariates at the same time. LongDat is able to take in several data types as input, including count, proportion, binary, ordinal and continuous data. The output table contains p values, effect sizes and covariates of each feature, making the downstream analysis easy. 


## Install
Install LongDat by typing ```install.packages("LongDat")``` in R. 

If you encounter errors like the one below when installing the package \
```Error: package or namespace load failed for ‘LongDat’ object ‘A’ is not exported by 'namespace:B_package'``` \
please try install the dependency B_package first, and then try to install LongDat again. An example to this kind of problem and solution can be found [here](https://stackoverflow.com/questions/48962946/error-package-or-namespace-load-failed-for-arulesviz-object-cividis-is-not)

## Tutorial
Tutorials for the analysis on continuous time variable (e.g. days) can be found [here](https://CRAN.R-project.org/package=LongDat/vignettes/LongDat_cont_tutorial.html). 

Tutorials for the analysis on discrete time variable (e.g. before/after treatment) can be found [here](https://CRAN.R-project.org/package=LongDat/vignettes/LongDat_disc_tutorial.html). 

Alternatively, you can type ```browseVignettes(“LongDat”)``` in R after installing LongDat to access these tutorials.

## Citation
The paper will be added here once it is published. Before that, please cite:  \
Chia-Yu Chen, Ulrike Löber, Sofia K Forslund, LongDat: an R package for covariate-sensitive longitudinal analysis of high-dimensional data, Bioinformatics Advances, Volume 3, Issue 1, 2023, vbad063, https://doi.org/10.1093/bioadv/vbad063
