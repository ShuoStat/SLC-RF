
#-------------------------------------------------------------------------------

biY <- function(y, cutTimes = NULL) {
  
  # times: choice for time cutoffs (years), if NULL, all times will be tried.
  if (is.null(cutTimes))
    cutTimes <- unique(sort(y[ ,1], decreasing = F))
  
  dif = 99999
  bestCut <- 99999
  
  for(i in cutTimes){
    
    time <- y[ ,1]
    status <- y[ ,2]
    
    status <- ifelse(time > i, 1, status)
    time   <- time[status != 0]
    
    #- rm zero
    tmp   <- abs(sum(time > i) - sum(time <= i))
    dif   <- ifelse(dif < tmp, dif, tmp)
    bestCut <- ifelse(dif < tmp, bestCut, i)
  }
  
  s <- y[,1] >= bestCut | y[,2] != 0
  # s <- y[,1] < bestCut & y[,2] != 0
  
  list(biY = y[s, 1] < bestCut,
       sel = s,
       bestCut = bestCut)
}


getData <- function(dir, datName, log2 = T, toBinary = T,
                    cutTimes = NULL) {
  
  # datName, data name
  # log2, whether log2 transformation on TPM (molecular information)
  # cutTimes, a vector of years. If null, all times will be used
  
  datDir  <- paste0(dir, datName, ".RData")
  load(datDir)
  
  y <- survival::Surv(X[ ,"time"], event = X[ ,"status"])
  X <- X[ ,-c(1:3)]
  
  b <- as.numeric(!grepl("clin", colnames(X))) + 1
  
  if (log2) {
    X[ ,b == 2] <- log2(X[ ,b == 2] + 0.001)
  }
  
  if (is.null(cutTimes))
    cutTimes <- sort(y[ ,1], decreasing = F)
  else
    cutTimes <- cutTimes * 365
  
  #- trans to binary response
  if (toBinary){
    
    # cutTime, from years to days
    obj <- biY(y, cutTimes = cutTimes)
    X <- X[obj$sel,]
    
    bestCut <- obj$bestCut
    if (!is.null(cutTimes))
      bestCut = bestCut / 365
    
    surv <- y[obj$sel]
    y <- obj$biY
  }
  
  #- Rename molecules since the original name is too long
  colnames(X) <- c(colnames(X)[b == 1], 
                   paste0("m", seq_len(sum(b == 2))))
  
  #- return
  .GlobalEnv$X <- as.matrix(X)
  .GlobalEnv$y <- y
  .GlobalEnv$blocks <- b
  if (toBinary) {
    .GlobalEnv$bestCut <- bestCut
    .GlobalEnv$surv <- surv 
  }
}

getIPCW <- function(X.t, X.v, surv.t, surv.v, bestCut) {
  
  censT <- Surv(surv.t[ ,1], 1 - surv.t[ ,2])
  censV <- Surv(surv.v[ ,1], 1 - surv.v[ ,2])
  # fit cox model
  mod  <- coxph(censT ~ ., data = as.data.frame(X.t))
  lp.t <- predict(mod, as.data.frame(X.t))
  lp.v <- predict(mod, as.data.frame(X.v))
  
  # high, high risk groups. 1, censoring
  s <- surv.v[,1] >= bestCut | surv.v[,2] == 1
  # if (any(!s)) {
  #   warning("Several early censored cases will have zero weights")
  # }
  # 
  # high risk group
  high <- surv.v[, 1] < bestCut
  
  # get times: high time + bestCut
  # time should be in increasing order
  times <- sort(c(censV[, 1], bestCut), decreasing = F)
  
  # get cumulative hazard
  chz <- cumhz(lp.t, censT, times = times)
  probs <- mapply(function(x) {exp(-x * exp(lp.v))}, x = chz, SIMPLIFY = T)
  
  # get probabilities
  timeUse <- ifelse(high, surv.v[, 1], bestCut)
  ind     <- match(timeUse, times)
  
  p <- c()
  for(i in seq_along(timeUse)) {
    p <- c(p, probs[i, ind[i]])
  }
  
  weights <- ifelse(s, 1 / p, 0)
  setNames(weights, rownames(X.v))
}

getSplits <- function(n, nsim, nfolds, seed, foldid.internal = T){
  
  # n, sample size
  # nsim, number of duplicated runs
  
  set.seed(seed)
  
  sampleSplits <- list()
  foldids <- list()
  
  for(i in seq_len(nsim)){
    samsplit <- sample(rep(1:nfolds, length = n), n)
    sampleSplits[[i]] <- samsplit
    foldids[[i]] <- sapply(1:nfolds, 
                           function(id){ 
                             nlength = sum(samsplit != id)
                             sample(rep(1:10, length = nlength), nlength)}, 
                           simplify = F)
  }
  
  if (nsim == 1) {
    sampleSplits <- unlist(sampleSplits)
    foldids <- foldids[[1]]
  }
  
  if (foldid.internal) {
    
    foldid.internals <- list()
    for(i in seq_len(nsim)){
      tmp <- list()
      for(j in seq_along(foldids[[i]])){
        sams <- foldids[[i]][[j]]
        tmp[[j]] <- sapply(1:nfolds, 
                           function(x){ 
                             nlength = sum(sams != x)
                             sample(rep(1:10, length = nlength), nlength)
                           }, simplify = F)
      }
      foldid.internals[[i]] <- tmp
    }
    .GlobalEnv$foldid.internals <- foldid.internals
  }
  
  #- return
  .GlobalEnv$sampleSplits <- sampleSplits
  .GlobalEnv$foldids <- foldids
  
}


cumhz <- function(lp, y, times){
  
  require(survival)
  if (!is.Surv(y))
    stop("y should be a Surv object")
  
  n <- length(y)
  ord <- order(y[,1], y[,2])
  
  #- time should be ordered
  #- warning: Strata will not be considered
  
  riskRegression::baseHaz_cpp(starttimes = rep(0, n),
                              stoptimes  = y[ord ,1],
                              status = y[ord ,2],
                              eXb = exp(lp)[ord],
                              strata = rep(0, n),
                              nPatients = n,
                              nStrata = 1,
                              emaxtimes = max(y[,1]),
                              predtimes = times,
                              cause = 1,
                              Efron = T)$cumhazard
  
}


#-------------------------------------------------------------------------------

getAUC <- function(haty, y, weights = NULL) {
  
  # remove zero weights
  if (is.null(weights))
    weights <- rep(1, length(y))
  
  indZero <- weights == 0
  roc <- WeightedROC::WeightedROC(as.numeric(haty)[!indZero], 
                                  as.numeric(y)[!indZero], 
                                  weight = weights[!indZero])
  auc <- WeightedROC::WeightedAUC(roc)
  
  # roc <- pROC::roc(as.numeric(y)[!indZero], as.numeric(haty)[!indZero])
  return(as.numeric(auc))
}

getBrier <- function(haty, y, weights = NULL) {
  
  if (is.null(weights))
    weights <- rep(1, length(y))
  
  haty <- drop(haty)
  brier <- sum(weights * (haty - y)^2) / sum(weights)
  return(brier)
}

getDev <- function(haty, y, weights) {
  
  if (is.null(weights))
    weights <- rep(1, length(y))
  
  haty <- drop(haty)
  dev <- -2 * (sum(y * log(haty) * weights) + sum(((1- y)) * log((1 - haty)) * weights))
  return(dev)
}

getAcc <- function(haty, y, weights = NULL){
  
  if (is.null(weights))
    weights <- rep(1, length(y))
  
  haty <- ifelse(haty > 0.5, 1, 0)
  1 - sum(abs(haty - y) * weights) / sum(weights)
  
}


getF1 <- function(haty, y) {
  
  # haty is prediction classes
  recall <- sum(haty + y == 2) / sum(y == 1)
  precision <- sum(y + haty == 2) / sum(haty == 1)
  2 * (recall * precision) / (recall + precision)
  
}


