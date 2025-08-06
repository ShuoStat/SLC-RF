
# Script for running the simulations

#-------------------------------------------------------------------------------
# environment

library(doParallel)
library(doRNG)
library(doSNOW)

source("./fitFun.R")
source("./simFun.R")
source("../functions/helpFun.R")


# simulation parameters
nCores = 25  # number of cores
rep = 100    # repetitions
n_genes = 1000 # number of genes
n_path = 50  # nunbre of pathways
rho_b = 0.3  # correlation between blocks
rho_w <- c(0.3, 0.5, 0.7) # correlation within blocks
n_obs <- c(100, 400)  # 
nonLinear <- c(TRUE, FALSE)

#-------------------------------------------------------------------------------
# set parallel running

if (!"sim" %in% dir("../../output/")){
  dir.create("../../output/sim/")
} else {
  unlink("../../output/sim/", recursive = TRUE)
}

#-------------------------------------------------------------------------------

# parallel runing
nCores <- pmin(detectCores() - 2, nCores)
cl <- makeCluster(nCores)
registerDoSNOW(cl)

scenario = 1
for(rho_w in c(0.3, 0.5, 0.7)) {
  
  for(n_obs in c(100, 400)) {
    
    for(nonLinear in c(TRUE, FALSE)) {
      
      message(scenario, " in a total of 12 scenarios! \n",
              "rho_w = ", rho_w, "\n", 
              "n_obs = ", n_obs, "\n",
              "nonLinear = ", nonLinear, "\n")
      
      scenario = scenario + 1

      # newX and newy
      newdata <- genX(rho_w, rho_b, n_genes, n_path, n_obs = 10000)
      new     <- genY(newdata, nonLinear = nonLinear)
      newy    <- new$y
      newProb <- new$prob
      
      # processor bars
      pb <- txtProgressBar(min = 1, max = rep, style = 3)
      progress <- function(k) setTxtProgressBar(pb, k)
      opts <- list(progress = progress)
      
      set.seed(7569431)
      
      out <- foreach(times = 1:rep, 
                     .packages = c("MASS"),
                     # seed = rng[(i - 1) * rep + 1:rep]
                     .options.snow = opts) %dorng% {
                       
                       X <- genX(rho_w, rho_b, n_genes, n_path, n_obs = n_obs)
                       y <- genY(X, nonLinear = nonLinear)$y
                       
                       # test data
                       preds <- fit.fun(X, y, newdata = newdata, 
                                        ntree = 100, 
                                        mtry = sqrt(ncol(X)), 
                                        nodesize = 1, 
                                        max_depth = 50)
                       
                       probs <- lapply(preds, `[[`, "mod")
                       unlist(lapply(probs, getAUC, y = newy))
                       
                       class <- lapply(preds, `[[`, "class")
                       unlist(lapply(probs, getAcc, y = newy))
                       
                       mods <- list(preds = preds, y = newy, prob = newProb)
                       nam <- paste0(c(n_obs, rho_w, nonLinear, times), collapse = "_")
                       save(list = "mods", file = paste0("../../output/sim/sim_", nam, ".RData"))
                       
                     }
    }
  }
}

stopCluster(cl)

#-------------------------------------------------------------------------------






























