
########################################################################################################################
##  Date: September, 2022                                                                                             ##
##                                                                                                                    ##
##  Boosting predictive ability of tropical maize hybrids via genotype-by-environment interaction under multivariate  ##
##  GBLUP models                                                                                                      ## 
##  DOI: 10.1002/csc2.20253                                                                                           ##
##  Krause et al. (2020)                                                                                              ##
##                                                                                                                    ##
##  Author:    MD Krause        <krause.d.matheus@gmail.com>                                                          ##
########################################################################################################################

########################################################################################################################
##                                                  DISCLAIMER                                                        ##
##             I honestly do not recall If I wrote these functions, modified them, or simply found them               ##
##         somewhere else. If you identify yourself as the author, please contact me so I can cite you.               ##
########################################################################################################################


SamplingLines <- function(id,I,seed=NULL,rep=NULL){
        if(is.null(seed)) seed <-1010231
        if(is.null(rep)) rep<-10
        if(is.null(I)) I <-.2
        n <- length(id)
        set1 <- list()
        set.seed(seed)
        for(s in 1:rep) set1[[s]] <- caret::createFolds(1:length(id), k = round(n/(n*I)), list = FALSE, returnTrain = FALSE);set.seed(seed+s)
        return(set1)
}

SamplingEnv <- function(envs,N,seed=NULL,rep=NULL){
        if(is.null(seed)) seed <-1010231
        if(is.null(rep)) rep<-10
        if(is.null(N)) N <-1
        set1 <- list()
        set.seed(seed)
        for(s in 1:rep) set1[[s]] <- sample(envs, N, replace = F);set.seed(seed+s)
        return(set1)
}
