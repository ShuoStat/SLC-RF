
source("../functions/helpFun.R")
source("./fitFun.R")

library(glmnet)
library(dplyr)
library(doParallel)
library(doSNOW)
library(doRNG)
library(tidyr)

#-------------------------------------------------------------------------------
#- run it
#- Load data
#-------------------------------------------------------------------------------

nsim   <- 10 # number of duplications
nfolds <- 5  # five-fold CV
nCores <- 25 # parallel runing

# used data
datNames <- c("BLCA", "BRCA", "COAD", "HNSC", "KIRC", "LIHC", 
              "LGG",  "LUAD",  "LUSC", "OV", "SKCM", "STAD")

# get basic information and dichotomize y. 

bestCuts <- c()
nVars <- c()
nPoss <- c()
nSams <- c()

for(nam in datNames) {
  
  cat("\n\n", nam, ":", as.character(Sys.time()), "\n\n")
  #- get data
  getData(dir = "../../data/", nam, log2 = F, toBinary  = T, cutTimes = c(1, 2, 3, 5))
  nSams <- c(nSams, length(y))
  nVars <- c(nVars, sum(blocks == 2))
  nPoss <- c(nPoss, sum(y == TRUE))
  bestCuts <- c(bestCuts, bestCut)
  
}
baseInfo <- data.frame(data = datNames, nSams = nSams, nPoss = nPoss, nVars = nVars, bestCuts = bestCuts)

write.csv(baseInfo, file = "../../results/baseInfo.csv")

names(bestCuts) <- datNames

#-------------------------------------------------------------------------------
#- Clinical information, Table 1
#-------------------------------------------------------------------------------


# fold change
# create directory
if (!"realData" %in% dir("../../output/")){
  dir.create("../../output/realData/")
} else {
  unlink("../../output/realData/", recursive = TRUE)
}

for(nam in datNames) {
  
  cat("\n\n", nam, ":", as.character(Sys.time()), "\n\n")
  #- get data
  getData("../../data/", nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  n <- nrow(X)
  getSplits(n, nsim, nfolds, seed = 7513564, foldid.internal = TRUE)
  
  set.seed(1584635)
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(min = 1, max = nfolds * nsim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  foreach(i = 1 : (nfolds * nsim), 
          .packages = c("glmnet", "randomForest", "rngtools"), 
          .options.snow = opts) %dorng% {
            
            samID <- ceiling(i / nfolds)
            #- j, which fold of the i.nsim split being the test
            j <- i - (samID - 1) * nfolds
            #- foldid for i.nsim splits and j fold
            
            sam <- sampleSplits[[samID]]
            
            # include RNAseq data only
            X.t <- X[sam != j, blocks == 2]
            y.t <- y[sam != j]
            X.v <- X[sam == j, blocks == 2]
            y.v <- y[sam == j]
            
            # generate binary outcomes
            grp <- biY(y.t, cutTimes = bestCuts[nam] * 365)
            X.t <- X.t[grp$sel,]
            y.t <- as.numeric(grp$biY)
            
            # get validation data X.v y.v
            grp <- biY(y.v, cutTimes = bestCuts[nam] * 365)
            X.v <- X.v[grp$sel,]
            y.v <- as.numeric(grp$biY)
            
            # models
            pred <- fit.fun(X = X.t, y = y.t, 
                            newdata = X.v, 
                            ntree = 100, 
                            mtry = round(sqrt(ncol(X.t))), 
                            nodesize = 1,
                            max_depth = 100)
            
            # add y.v to the output
            pred <- list(pred = pred, y.v = y.v)
            save(list = "pred", file = paste0("../../output/realData/real_", nam, "_", samID, "_", j, ".RData"))
          }
  
  stopCluster(cl)
  Sys.time()
}

#-------------------------------------------------------------------------------


