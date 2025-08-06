# linear combination 
LCRandomForest <- function(X, y, 
                           subsetSize = nrow(X) / 10, 
                           ngen = 1000,
                           ntree = 100, 
                           mtry = round(sqrt(ngen)), 
                           maxnodes = 10, 
                           nodesize = 1) {
    
    
    # generate newX
    n <- nrow(X)
    p <- ncol(X)
    
    # recode mode
    newX <- c()
    submods <- list()
    pIDs <- list()
    
    for(i in seq(ngen)) {
        
        # random faeture
        pID <- sample(seq(p), subsetSize)
        # random sample
        nID <- sample(seq(n), n, replace = TRUE)
        
        mod <- glm.fit(cbind(1, X[nID, pID]), y[nID], family = binomial())
        betas <- coef(mod)
        newX <- cbind(newX, cbind(1, X[,pID]) %*% betas)
        submods[[i]] <- betas
        pIDs[[i]] <- pID
        
    }
    
    require(randomForest)
    rf <- randomForest(x = newX, y = as.factor(y), 
                       ntree = ntree, 
                       mtry = mtry, 
                       maxnodes  = maxnodes, 
                       nodesize  = nodesize)
    
    re <- list(rf = rf,
               submods = submods,
               pIDs = pIDs)
    
    return(re)
}

predict.LCRandomForest <- function(obj, newdata, type = "response"){
    
    pIDs <- obj$pIDs
    newX <- c()
    
    for(i in seq_along(pIDs)) {
        pred <- cbind(1, newdata[,pIDs[[i]], drop = FALSE]) %*% obj$submods[[i]]
        newX <- cbind(newX, pred)
    }
    predict(obj$rf, newX, type = type)
}
