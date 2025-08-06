
source("../functions/helpFun.R")
source("./fitFun.R")

library(glmnet)
library(dplyr)
library(doParallel)
library(doSNOW)
library(tidyr)
library(doRNG)

#-------------------------------------------------------------------------------
#- run it
#- Load data
#-------------------------------------------------------------------------------

nsim   <- 10
nfolds <- 5
nCores <- 25
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
# summarize prediction performance 
#-------------------------------------------------------------------------------

library(MLmetrics)

aucs <- list()
accuracies <- list()

for(nam in datNames) {
  
  auc <- c()
  brier <- c()
  accuracy <- c()
  F1 <- c()
  
  for(samID in 1:nsim) {
    
    for(j in 1:nfolds) {
      
      load(paste0("../../output/realData/real_", nam, "_", samID, "_", j, ".RData"))
      preds <- pred$pred
      # preds$rf <- preds$rf$probabilities
      tmp <- unlist(lapply(preds, function(x) getAUC(x$mod, pred$y.v)))
      auc <- rbind(auc, tmp)
      
      tmp <- unlist(lapply(preds, function(x) getAcc(x$class, pred$y.v, weights = NULL)))
      accuracy <- rbind(accuracy, tmp)
      
    }
  }
  aucs[[nam]] <- auc
  accuracies[[nam]] <- accuracy
}

# import clinical information
baseInfo <- read.csv("../../results/baseInfo.csv")

# output AUC
out_auc <- t(as.data.frame(lapply(aucs, colMeans))) %>% 
  as.data.frame() %>%
  mutate(across(everything(), ~ sprintf("%.3f",.x))) %>%
  tibble::rownames_to_column(var = "Data") %>%
  left_join(., 
            select(baseInfo, data, nSams), 
            by = join_by(Data == data)) %>%
  arrange(desc(-nSams))

openxlsx::write.xlsx(out_auc, asTable = FALSE, file = paste0("../../results/real_auc.xlsx"))


# output accuracy 
out_acc <- t(as.data.frame(lapply(accuracies, colMeans))) %>% 
  as.data.frame() %>%
  mutate(across(everything(), ~ sprintf("%.3f",.x))) %>%
  tibble::rownames_to_column(var = "Data") %>%
  left_join(., 
            select(baseInfo, data, nSams), 
            by = join_by(Data == data)) %>%
  arrange(desc(-nSams))

openxlsx::write.xlsx(out_acc, asTable = FALSE, file = paste0("../../results/real_acc.xlsx"))

#-------------------------------------------------------------------------------
# computational time
#-------------------------------------------------------------------------------

times <- list()

for(nam in datNames) {
  
  time <- c()
  
  for(samID in 1:nsim) {
    
    for(j in 1:nfolds) {
      
      load(paste0("../../output/realData/real_", nam, "_", samID, "_", j, ".RData"))
      preds <- pred$pred
      # preds$rf <- preds$rf$probabilities
      tmp <- unlist(lapply(preds, function(x) as.numeric(x$time, units = "secs")))
      time <- rbind(time, tmp)
    }
  }
  
  times[[nam]] <- time
}

# data BLCA
time1 <- colMeans(times$BLCA)

dat <- data.frame(method = names(time1), time = time1) %>% 
  arrange(time) %>%
  filter(method != "lcRF_ridge") %>%
  mutate(method = recode(method,
                         "las" = "Lasso",
                         "xgb" = "XGboosting",
                         "svm" = "SVM",
                         "RF" = "RF",
                         "extraRF" = "Extreme Trees",
                         "ferns" = "Random Ferns",
                         "oblRF" = "ORF",
                         "lcRF_glm" = "SLC-RF"),
         method = factor(method, method))
  

library(ggplot2)
p1 <- ggplot(dat, aes(x = method, y = time)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) +
  labs(x = NULL, title = "A, BLCA", y = "Time/Secs") + 
  coord_flip()

# data BRCA
time2 <- colMeans(times$BRCA)

dat <- data.frame(method = names(time2), time = time2) %>% 
  arrange(time) %>%
  filter(method != "lcRF_ridge") %>%
  mutate(method = recode(method,
                         "las" = "Lasso",
                         "xgb" = "XGboosting",
                         "svm" = "SVM",
                         "RF" = "RF",
                         "extraRF" = "Extreme Trees",
                         "ferns" = "Random Ferns",
                         "oblRF" = "ORF",
                         "lcRF_glm" = "SLC-RF"),
         method = factor(method, method))

p2 <- ggplot(dat, aes(x = method, y = time)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) +
  labs(x = NULL, title = "B, BRCA", y = "Time/Secs") + 
  coord_flip()


p <- gridExtra::grid.arrange(p1, p2, nrow = 1)

ggsave(filename = "../../results/time.tiff", plot = p, units = "in",
       width = 7, height = 2.5, dpi = 600, compression = "lzw")


# end
#------------------------------------------------------------------------------





