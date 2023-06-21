# R package for "Integrative Network Learning for Multi-modaility Biomarker Data"

<img src="https://img.shields.io/badge/Study%20Status-Results%20Available-yellow.svg" alt="Study Status: Results Available"> 

The biomarker networks measured by different modalities of data (e.g., structural magnetic resonance imaging (sMRI), diffusion tensor imaging (DTI)) may share the same true underlying biological model. In this work, we propose a node-wise biomarker graphical model to leverage the shared mechanism between multi-modality data to provide a more reliable estimation of the target modality network and account for the heterogeneity in networks due to differences between subjects and networks of external modality. Latent variables are introduced to represent the shared unobserved biological network and the information from the external modality is incorporated to model the distribution of the underlying biological network. An approximation approach is used to calculate the posterior expectations of latent variables to reduce time.  

- Title: **Integrative Network Learning for Multi-modaility Biomarker Data**
  
- Authors: **Shanghong Xie<sup>a</sup> (sx2168@columbia.edu), Donglin Zeng<sup>b</sup>, and Yuanjia Wang<sup>a</sup>**
- Affiliations: 
  + 1. **Department of Biostatistics, Mailman School of Public Health, Columbia University, New York**
  + 2. **Department of Biostatistics, University of North Carolina, Chapel Hill, North Carolina**
  
- Manuscript: Xie S, Zeng D and Wang Y (2021). [Integrative Network Learning for Multi-modaility Biomarker Data.](https://github.com/shanghongxie/INL) Annals of Applied Statistics 15(1), 64-87. 
  

<p align="center">
<img src="https://github.com/shanghongxie/Integrative-Network/blob/master/Diagram1-1.png" width="900" height="500">
</p>





## Setup Requirements
- R
- Install Rcpp and RcppEigen packages

## Code Instructions

The code for the proposed methodology is included in **Code** folder. Please download all the files in the folder to implement the method.

### To implement the proposed method with approximated posterior expectation of latent variables, source the following R files
  
  sourceCpp('~/INLApproxC.cpp')
  
  source('~/INLApproxRcode.R')
  
  source('~/INLApproxHardThrRcode.R')
  
### To implement the proposed method with exact calculation of posterior expectation of latent variables, source the following R files
  
  sourceCpp('~/INLDirectC.cpp')
  
  source('~/INLDirectRcode.R')
  
  source('~/INLDirectHardThrRcode.R')

### Main functions

- **INLApproxRcode**: estimate network without pruning, the posterior expectations of latent variables are calculated by the approximated approach. It requires sourceCpp('~/INLApproxC.cpp').
  
- **INLApproxHardThrRcode.R**: hard thresholding the estimated network from INLApproxRcode output based on EBIC.
  
- **INLDirectRcode.R**: estimate network without pruning, the posterior expectations of latent variables are calculated by direct calculation. It requires sourceCpp('~/INLDirectC.cpp').
  
- **INLDirectHardThrRcode.R**: hard thresholding the estimated network from INLDirectRcode output based on EBIC.

### Examples

**Viegette.R** provides an example to implement the methods using the simulated data.

