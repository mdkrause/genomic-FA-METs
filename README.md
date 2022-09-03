# Boosting predictive ability of tropical maize hybrids viagenotype-by-environment interaction under multivariateGBLUP models

Krause, et al. (2020); https://doi.org/10.1002/csc2.20253

This manuscript was published in Crop Science in 2020, and it is the research work I have done for my master's degree at the University of São Paulo (2016-2017). After several requests, I decided to share a simplified version of the R code used for the project. This shared script differs from the material & methods presented in the paper as follows:

1. Genomic prediction in multi-environment trial (METs) is now performed with a two stage-wise model, with weights equal to 1/SE^2.
2. In addition to the 5-fold cross-validation (20% of missing hybrids) with 10 replicates, the new script automatically performs analysis for 5%, 10%, 15%, 20%, 25%, 30%, 35%, 40%, and 50% of missing hybrids within environments, with 100 replicates each. Parallel computing was employed. 
3. The first model assumes a compound symmetry covariance matrix, and the second a factor analysis (FA). Both include the genomic relationship (GBLUP model).
4. In the manuscript, we presented a way to include checks and non-single-hybrid cross materials as fixed effects. There is no magic; one needs to include a dummy variable in the model (see column X1 in the shared phenotypic data). In the shared script, I removed those genotypes and am just considering the 147 single-cross hybrids this time.

The following files contain:

1. pheno.rds: the adjusted means (eBLUE values) in a trial-level.
2. Ginv: the inverse of the additive genomic relationship matrix among single-cross hybrids.
3. CV1_CV2.R: codes to perform CV1/CV2 cross-validations schemes. 
4. predFA.R: codes to perform genomic prediction and compute basic results.
5. results.pdf: basic plots with results.

