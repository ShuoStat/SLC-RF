
# variance and covariance matrix
# 
# pre-define
# 10 genes in each pathways
# coefficients are fixed 

genX <- function(rho_w, rho_b, n_genes, n_path, n_obs) {
  
  # generate correlation matrix
  m <- matrix(rho_b, n_path * 10,  n_path * 10)
  
  for(i in seq(n_path)) {
    ind <- (i * 10 - 9):(i * 10)
    m[ind, ind] <- rho_w
  }
  
  diag(m) <- 1
  
  # mean and variance of x
  mu <- rep(0, n_path * 10)
  sigma <- 0.65
  
  # generate genes in pathways
  d1 <- mvtnorm::rmvnorm(n_obs, mu, m * sigma)
  
  # generate pathway information
  d2 <- matrix(rnorm(n_obs * (n_genes - n_path * 10), 0, sigma),
               n_obs,  (n_genes - n_path * 10))
  
  # combine data
  genes <- cbind(d1, d2)
  colnames(genes) <- paste0("m", seq(n_genes))
  
  # add noise
  # E <- matrix(rnorm(n_genes * n_obs, nu, tau), ncol = n_genes, nrow = n_obs)
  # multiplicative noise
  # M <- matrix(rnorm(n_genes * n_obs, 0, phi), ncol = n_genes, nrow = n_obs)
  # observations including the noise
  # genes <- genes * exp(M) + E
  
  # thresholding and normalization
  # genes[genes < 10] <- 10
  # genes[genes > 16000] <- 16000
  # genes <- log(genes)
  
  return(genes)
}

# X <- genX(rho_w = 0, rho_b = 0, n_genes = 1000, n_path = 50, n_obs = 1000)
# X <- genX(rho_w = 0.7, rho_b = 0.3, n_genes = 1000, n_path = 50, n_obs = 10000)
# X <- genX(rho_w = 0.5, rho_b = 0.3, n_genes = 1000, n_path = 50, n_obs = 10000)
# X <- genX(rho_w = 0.3, rho_b = 0.3, n_genes = 1000, n_path = 50, n_obs = 10000)
# gen Y

# linear = 2, location = 0
# nonLinear =  1.5, location = 0.75

genY <- function(X, nonLinear = TRUE){

  n_genes <- ncol(X)
  n_obs   <- nrow(X)

  # beta <- c(0.75, 0.75, rep(0, 8),
  #           0.5, rep(0, 9),
  #           0.25, 0.25, 0.25, rep(0, 7),
  #           0.75, 0.25, rep(0, 8),
  #           0.5, 0.5, 0.25, rep(0, 7),
  #           0.5, rep(0, 9),
  #           0.75, rep(0, 9),
  #           rep(0, (n_genes - 70)))
  
  beta <- c(3, -3, 1, -1, 1, rep(0, 5),
            2, -2, 2, rep(0, 7),
            rep(0, 10),
            -3, 2, 1, rep(0, 7),
            -1, -1, -1, -1, -1, rep(0, 5),
            -2, rep(0, 9),
            rep(1, 10),
            1, rep(0, 9),
            rep(0, (n_genes - 80)))
            
  # 
  # if (nonLinear) {
  #   # add non-linear
  #   X[, 31] <- exp(-X[, 31]^2)
  #   X[, 41] <- X[,41]^2
  #   X[, 51] <- abs(X[,51])
  #   X[, 61] <- sin(X[, 61] * pi)
  # }
  
  
  if (nonLinear) {
    # add non-linear
    X[, 31:33] <- exp(-X[, 31:33]^2)
    X[, 41:45] <- X[,41:45]^2
    X[, 51] <- sin(X[, 51] * pi)
    X[, 61:63] <- abs(X[,61:63])
    X[, 64:65] <- sin(X[,64:65] * pi)
    
  }
  
  eta <- X %*% beta
  mean(eta)
  sd(eta)
  
  par_location <- ifelse(nonLinear, -1.35, 0)
  par_scale <- ifelse(nonLinear, 3, 5)
  p <- plogis(eta, par_location, par_scale)
  hist(p)
  y <- rbinom(n_obs, size = 1, prob = p)
  
  # return y and p
  re <- list(y = y, prob = p[,1])
  return(re)
}

#-------------------------------------------------------------------------------












































