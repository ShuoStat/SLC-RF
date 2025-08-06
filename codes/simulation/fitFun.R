
# 
fit.fun <- function(X, y, newdata, 
                    ntree = 100, 
                    mtry = sqrt(ncol(X)), 
                    nodesize = 1, 
                    max_depth = 100) {
  
  # lasso, baseline
  require(glmnet)
  
  # lasso
  a <- Sys.time()
  cvs <- cv.glmnet(X, y, family = "binomial")
  las <- glmnet(X, y, family = "binomial", lambda = cvs$lambda.min)
  pred.las <- as.numeric(predict(las, newdata, type = "response"))
  b <- Sys.time()
  pred.las <- list(mod = pred.las,
                   class = as.numeric(pred.las > 0.5),
                   time = b - a)
  
  
  # xgboosting
  a <- Sys.time()
  require(xgboost)
  params <- list(
    objective = "binary:logistic",  # Binary classification
    eval_metric = "logloss",        # Logarithmic loss (alternatives: "auc", "error")
    eta = 0.3,                      # Learning rate (lower = more robust, slower)
    max_depth = 6, 
    min_child_weight = 1,          # Tree depth (lower = less overfitting)
    nthread = 1                    # Parallel processing
  )
  
  dtrain <- xgb.DMatrix(data = X, label = y)
  xgb <-xgb.train(params,
                  data = dtrain,
                  nrounds = 100, # Number of boosting iterations
                  verbose = 1)
  pred.xgb <- predict(xgb, newdata) 
  b <- Sys.time()
  pred.xgb <- list(mod = pred.xgb,
                   class = as.numeric(pred.xgb > 0.5),
                   time = b - a)
  
  # support vector machine
  
  a <- Sys.time()
  require(e1071)
  svm_model <- svm(x = X,
                   y = y,
                   type = "C-classification",  # For classification
                   kernel = "linear",          # Linear kernel for high-dim data
                   scale = TRUE,               # Scale features (important for SVM)
                   cost = 1,
                   probability = TRUE)                   # Regularization parameter (tune this)
  
  pred.svm <- predict(svm_model, newdata, probability = TRUE)
  pred.svm <- attr(pred.svm, "probabilities")[,2]
  b <- Sys.time()
  pred.svm <- list(mod = pred.svm,
                   class = as.numeric(pred.svm > 0.5),
                   time = b - a)
  
  # random forest
  # rf <- random_forest(X, y, n_trees = ntree, max_depth = max_depth, min_samples_leaf = nodesize, n_features = mtry)
  # pred.rf <- predict_forest(rf, newdata, type = "probabilities")$probabilities
  a <- Sys.time()
  require(randomForest)
  RF <- randomForest(x = X, y = as.factor(y), 
                     ntree = ntree, mtry = mtry, 
                     maxnodes  = max_depth, 
                     nodesize  = nodesize)
  pred.RF <- predict(RF, newdata = newdata, type = "prob")[,2] 
  b <- Sys.time()
  pred.RF <- list(mod = pred.RF,
                  class = as.numeric(pred.RF > 0.5),
                  time = b - a)
  
  # Extreme random forest
  a <- Sys.time()
  require(extraTrees)
  extraRF <- extraTrees(x = X, y = as.factor(y), 
                        ntree = ntree,
                        mtry = mtry,
                        nodesize = nodesize,
                        numRandomCuts = 1,
                        evenCuts = FALSE,
                        numThreads = 1,
                        quantile = F,
                        weights = NULL,
                        subsetSizes = NULL,
                        subsetGroups = NULL,
                        tasks = NULL,
                        numRandomTaskCuts = 1,
                        na.action = "stop")
  
  pred.extraRF <- predict(extraRF, newdata = newdata, probability = TRUE)[ ,2]
  b <- Sys.time()
  pred.extraRF <- list(mod = pred.extraRF,
                       class = as.numeric(pred.extraRF > 0.5),
                       time = b - a)
  
  # Random Ferns Classifier
  a <- Sys.time()
  require(rFerns)
  ferns <- rFerns(x = X,
                 y = as.factor(y),
                 depth = 5,
                 ferns = ntree,
                 importance = "none",
                 saveForest = TRUE,
                 consistentSeed = NULL,
                 threads = 0)
  
  pred.ferns_score <- predict(ferns, x = as.data.frame(newdata), scores = TRUE)[,2]
  pred.ferns_class <- predict(ferns, x = as.data.frame(newdata), scores = FALSE)
  b <- Sys.time()
  pred.ferns <- list(mod = pred.ferns_score,
                     class = as.numeric(pred.ferns_class) - 1,
                     time = b - a)
  
  # # Oblique Random Forest
  a <- Sys.time()
  source("../functions/obliqueRF.R", local = TRUE)
  oblRF <- obliqueRF(x = X, y = y,
                     x.test = NULL,
                     y.test = NULL,
                     mtry = mtry,
                     ntree = ntree,
                     training_method = "ridge",
                     bImportance = FALSE,
                     bProximity = FALSE)
  pred.oblRF <- predict.obliqueRF(oblRF, newdata = newdata, type = "prob")[,2]
  b <- Sys.time()
  pred.oblRF <- list(mod = pred.oblRF,
                     class = as.numeric(pred.oblRF > 0.5),
                     time = b - a)
  # 

  # Linear combination random forest
  a <- Sys.time()
  source("../functions/LCRandomForest_v1.R", local = TRUE)
  LCRF <- LCRandomForest(X, y,
                         # method = "glm",
                         subsetSize = nrow(X) / 5,
                         ngen = 1000,
                         ntree = ntree,
                         maxnodes = max_depth,
                         nodesize = nodesize)
  pred.lcRF <- predict.LCRandomForest(LCRF, newdata = newdata, type = "prob")[,2]
  b <- Sys.time()
  pred.lcRF_glm <- list(mod = pred.lcRF,
                        class = as.numeric(pred.lcRF > 0.5),
                        time = b - a)


  # output
  re <- list(
    las = pred.las, 
    xgb = pred.xgb,
    svm = pred.svm,
    RF = pred.RF,
    extraRF = pred.extraRF,
    ferns = pred.ferns,
    oblRF = pred.oblRF,
    lcRF_glm = pred.lcRF_glm
    # lcRF_ridge = pred.lcRF_ridge
  )
  
  return(re)
}




