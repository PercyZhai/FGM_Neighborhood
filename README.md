# FGM_Neighborhood
Functional Graphical Models via Neighborhood Selection Approach

## ROC
This folder contains codes running multiple methods, namely FPCA-gX, FPCA-gY, FGLasso, PSKL, and FPCA-PSKL (for Model D).
Threshold is set to $\epsilon=0$. Penalty parameter $\lambda_n$ varies.
Output results are TPR and FPR for each methods.
Runtime of each method under each run is also recorded.
A script is given for plotting ROC for all methods under each setting.

## ROC Threshold
This folder contains codes running FPCA-gX method.
Threshold varies under each setting. Penalty parameter $\lambda_n$ varies under each threshold parameter.
Output results are TPR and FPR for FPCA-gX method under all thresholds
A script is given for plotting ROC for all thresholds under each setting.

## SCV
This folder contains codes running FPCA-gX method.
Threshold varies under each setting. Penalty parameter $\lambda_n$ varies under each threshold parameter.
Output results are adjacency matrices under AND and OR schemes under the optimal choice of $\lambda_n, t_\epsilon$ pairs that maximizes SCV-BIC.
Runtime of each run is also recorded.
A script is given for computing precision and recall for the optimal adjacency matrices.

## Real Data Analysis
This folder contains codes running FPCA-gX method.
Two datasets: ABIDE for autism and ADHD.
Both raw data and time series extracted from raw data can be found in Data folder.
Time series are extracted by "ABIDE.read.data.R" and "ADHD.read.data.R" respectively.
### ABIDE dataset:
"SCV.autism.R" selects and outputs the optimal adjacency matrices using SCV under FPCA-gX method.
"spa.ctrl.autism.R" outputs an adjacency matrix with 2\% connection density under FPCA-gX method.
### ADHD dataset:
"SCV.ADHD.R" selects and outputs the optimal adjacency matrices using SCV under FPCA-gX method.
"spa.ctrl.ADHD.R" outputs an adjacency matrix with 2\% connection density under FPCA-gX method.
"spa.ctrl.ADHD.FGLasso.R" outputs an adjacency matrix estimated by FGLasso method.
