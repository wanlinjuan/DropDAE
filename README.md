
# DropDAE <img src="man/figures/logo.png" align="right" width="120"/>

[![R-CMD-check](https://github.com/wanlinjuan/DropDAE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/wanlinjuan/DropDAE/actions)

## Introduction

**DropDAE** (Dropout Denoising Autoencoder) is a novel deep learning
framework that improves imputation for high-dimensional data like
single-cell RNA sequencing (scRNA-seq).  
It uses: 
- A **denoising autoencoder** structure to reconstruct corrupted input data.
- An optional **triplet loss** based on **consensus clustering** to improve separation between similar and
dissimilar cells.


## Installation

You can install the development version of `DropDAE` from GitHub:

\`\`\`r \# install.packages(“devtools”)
devtools::install_github(“wanlinjuan/DropDAE”)
