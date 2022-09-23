# Boosting predictive ability of tropical maize hybrids viagenotype-by-environment interaction under multivariate GBLUP models

Krause, et al. (2020); https://doi.org/10.1002/csc2.20253

This manuscript was published in Crop Science in 2020, and it is the research work I have done for my master's degree at the University of São Paulo (2016-2017). After several requests, I decided to share a simplified version of the R code used for this project. This shared script differs from the material & methods presented in the paper as follows:

1. Genomic prediction in multi-environment trial (METs) is now performed with a two stage-wise model, with weights equal to 1/SE^2.
2. In addition to the 5-fold cross-validation (20% of missing hybrids) with 10 replicates, the new script automatically performs analysis for 10%, 20%, and 40% of missing hybrids within environments, with 50 replicates each. Parallel computing was employed. 
3. The first model assumes a compound symmetry covariance matrix, and the second a factor analysis (FA). Both included the genomic relationship (GBLUP, yielding genomic estimated breeding values within environments).
4. In the manuscript, we presented a modelling way to include checks and non-single-hybrid cross genotypes as fixed effects. There is no magic; one needs to include a dummy variable in the model (see column X1 in the shared phenotypic data). In the shared script, I removed these genotypes and I am just considering the 147 single-cross hybrids.
5. In the manuscript, we considered 12 environments, and now another one was includes (i.e., there are 13 environments).

The following files contain:

1. pheno.rds: the adjusted means (eBLUE values) in a trial-level.
2. Ginv.rds: the inverse of the additive genomic relationship matrix among single-cross hybrids.
3. CV1_CV2.R: codes to perform CV1/CV2 cross-validations schemes. 
4. predFA.R: codes to perform genomic prediction within environments. Note the code was written to run in a batch jobs (i.e., server, with a .job file). Out of curiosity, it took 7.5 days to be completed in 32 Intel(R) Xeon(R) Gold 6152 CPU @2.10GHz processors with ~120 Gb of RAM memory.
5. resultsPred.html: basic R codes & plot with results (similar to Figure 4 in the manuscript).
