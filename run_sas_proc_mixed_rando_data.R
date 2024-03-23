#### RUN SAS PROC MIXED ON 1000 RANDOMISED DATA SETS

### Using Python's saspy package to run PROC MIXED to estimate mutational 
### variance covariance matrices for 1000 randomised datasets.

# R packages and installing Python and SAS (via saspy)
library(reticulate)
py_available(TRUE)
py_install("saspy", pip = TRUE)
py_module_available("saspy")# check that saspy is available
