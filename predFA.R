
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

library(dplyr)
library(foreach)
library(doParallel)

I <- c(0.10, 0.20, 0.40)
rep <- 50

pheno <- readRDS("pheno.rds")
Ginv <- readRDS("Ginv.rds")
ids <- rownames(Ginv)
n <- length(ids)

pheno <- pheno %>% filter(geno %in% rownames(Ginv)) %>% droplevels()
envs <- levels(pheno$envs)

source("CV1_CV2.R")

doParallel::registerDoParallel(32) # carefull here

resultsFA <- c()

for(MISS in 1:length(envs)){
        
        for(STR in 1:length(I)){
                
                lines <- SamplingLines(id = ids, I = I[STR], rep = rep)
                env <- SamplingEnv(envs = envs, N = MISS, rep = rep)
                
                results1 <-foreach(REP = 1:rep)%:%
                        foreach(i=1:round(n/(n*I[STR])), 
                                .packages = c('asreml','dplyr'),
                                .combine=rbind)  %dopar% {
					.GlobalEnv$Ginv <- Ginv 
                               m <- asreml(fixed = eblue ~ envs,
                                       random = ~ corh(envs):vm(geno,Ginv),
                                       family = asr_gaussian(dispersion = 1),
                                       weights = weights,
                                       workspace = "2gb",
                                       trace = FALSE, maxit = 100,
                                       data = mutate_at(pheno, vars(eblue),
                                                        list(~ ifelse(geno %in% ids[which(lines[[REP]] == i)] &
                                                                              envs %in% env[[REP]], NA, .))))
                               
                               
                               pred2 <- data.frame(predict.asreml(m, "envs:geno", 
                                                                  pworkspace = "3gb")$pvals)
                                                                  
                               # If one also includes the mean of each environment given how the cross-validation was designed, 
                               # it will (greatly) increase the Pearson correlation between predicted GEBVs and observed eBLUE. 
                               # This won't be informative at all.
                               
                               pred1 <- pred2 %>% group_by(envs) %>% summarise(fixedEnv = mean(predicted.value)) %>% as.data.frame()
                               
                               pred2 <- pred2 %>% filter(geno %in% as.factor(ids[which(lines[[REP]] == i)]) &
                                                         envs %in% as.factor(env[[REP]]))
                               
                               
                               pred <- left_join(pred2, pred1)
                               pred$missEnv <- MISS
                               pred$missGeno <- I[STR]
                               pred$rep <- REP
                               pred$fold <- i
                               return(pred)
			}
                
                results1 <- do.call(rbind,results1)
                results1 <- left_join(results1, pheno[, c(1:3)])
                results1$eblup <- results1$eblue-results1$fixedEnv
                results1$model <- 'CS'
                
                results1 <- results1 %>% group_by(rep, model, missEnv, missGeno) %>% summarise(corr = cor(eblup, predicted.value))
                resultsFA <- rbind(resultsFA, results1) ; rm(results1)
                
                results2 <-foreach(REP = 1:rep)%:%
                        foreach(i=1:round(n/(n*I[STR])), 
                                .packages = c('asreml','dplyr'),
                                .combine=rbind)  %dopar% {
                                        .GlobalEnv$Ginv <- Ginv 
                                        m <- asreml(fixed = eblue ~ envs,
                                                    random = ~ fa(envs,2):vm(geno,Ginv),
                                                    family = asr_gaussian(dispersion = 1),
                                                    weights = weights,
                                                    workspace = "2gb",
                                                    trace = FALSE, maxit = 100,
                                                    data = mutate_at(pheno, vars(eblue),
                                                                     list(~ ifelse(geno %in% ids[which(lines[[REP]] == i)] &
                                                                                           envs %in% env[[REP]], NA, .))))
                                        
                                        pred2 <- data.frame(predict.asreml(m, "envs:geno", 
                                                                           pworkspace = "3gb")$pvals)
                                        
                                        pred1 <- pred2 %>% group_by(envs) %>% summarise(fixedEnv = mean(predicted.value)) %>% as.data.frame()
                                        
                                        pred2 <- pred2 %>% filter(geno %in% as.factor(ids[which(lines[[REP]] == i)]) &
                                                                          envs %in% as.factor(env[[REP]]))
                                        
                                        
                                        pred <- left_join(pred2, pred1)
                                        pred$missEnv <- MISS
                                        pred$missGeno <- I[STR]
                                        pred$rep <- REP
                                        pred$fold <- i
                                        return(pred)
                                }
                
                results2 <- do.call(rbind,results2)
                results2 <- left_join(results2, pheno[, c(1:3)])
                results2$eblup <- results2$eblue-results2$fixedEnv
                results2$model <- 'FA2'
                
                results2 <- results2 %>% group_by(rep, model, missEnv, missGeno) %>% summarise(corr = cor(eblup, predicted.value))
                resultsFA <- rbind(resultsFA, results2)  ; rm(results2)
                
        }
        
}
saveRDS(resultsFA, "resultsFA.rds")
