# Circular Outlier Detection Using Mixture of Von Mises Distributions

This repository contains R code related to **circular outlier detection** using a **mixture of Von Mises distributions** and the **Expectation-Maximization (EM) algorithm**.

## Research Context

This code is based on the MSc research work of **Khadijeh Abdi**.

**Thesis:** Outlier Detection in Circular Regression Models using Mixture of Von Mises Distributions  
**Publication:** Circular Outliers Detection Using a Mixture of Von Mises Distributions  
**Journal:** Journal of Applied Mathematics and Modeling (JAMM)  
**DOI:** 10.22055/JAMM.2018.20998.1394

## Main Methods

- Circular data analysis
- Von Mises distribution
- Mixture models
- EM algorithm
- Outlier / anomaly detection
- COVRATIO diagnostic calculation
- Simulation-based model evaluation

## Files

- `circular_outlier_detection_von_mises.R`  
  Main R script containing:
  - COVRATIO calculation
  - EM algorithm for a two-component Von Mises mixture model
  - Example usage sections

## Required R Package

```r
install.packages("CircStats")
library(CircStats)
```

## Important Note

The original program was available as a scanned thesis appendix.  
This version has been cleaned and formatted for GitHub.  
Before using it for publication, reproduction, or formal analysis, the code should be checked against the original thesis and tested carefully.

## Suggested Citation

Abdi, K., Golalizadeh, M., & Baghfalaki, T.  
*Circular Outliers Detection Using a Mixture of Von Mises Distributions*.  
Journal of Applied Mathematics and Modeling.  
DOI: 10.22055/JAMM.2018.20998.1394
