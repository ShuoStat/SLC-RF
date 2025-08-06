
# Summarize the simulation results

library(doParallel)
library(doRNG)
library(doSNOW)
library(tidyverse)

source("../functions/helpFun.R")

# simulation parameters
rep = 100
n_genes = 1000
n_path = 50
rho_b = 0.3
rho_w <- c(0.3, 0.5, 0.7)
n_obs <- c(100, 400)
nonLinear <- c(TRUE, FALSE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

aucs <- list()
accs <- list()
bries <- list()

scenario <- 1

for(rho_w in c(0.3, 0.5, 0.7)) {
    
    for(n_obs in c(100, 400)) {
        
        for(nonLinear in c(TRUE, FALSE)) {
            
            message(scenario, " in a total of 12 scenarios! \n",
                    "rho_w = ", rho_w, "\n", 
                    "n_obs = ", n_obs, "\n",
                    "nonLinear = ", nonLinear, "\n")
            
            scenario <- scenario + 1
            
            tmpAUC <- c()
            tmpACC <- c()
            tmpBrier <- c()
            
            for(times in 1:rep) {
                
                nam <- paste0(c(n_obs, rho_w, nonLinear, times), collapse = "_")
                load(paste0("../../output/sim/sim_", nam, ".RData"))
                
                haty <- lapply(mods$preds, `[[`, "mod")
                hatClass <- lapply(mods$preds, `[[`, "class")
                y <- mods$y
                
                tmpAUC <- rbind(tmpAUC, unlist(lapply(haty, getAUC, y = y)))
                tmpACC <- rbind(tmpACC, unlist(lapply(hatClass, getAcc, y = y)))
                tmpBrier <- rbind(tmpBrier, unlist(lapply(haty, getBrier, y = y)))
                
            }
            
            nam <- paste0(c("s", n_obs, rho_w, nonLinear), collapse = "_")
            aucs[[nam]] <- tmpAUC
            accs[[nam]] <- tmpACC
            bries[[nam]] <- tmpBrier
        }
    }
}


# re <- lapply(aucs, colMeans)
# round(t(as.data.frame(re)), 3)
# 
# ac <- lapply(accuracies, colMeans)
# round(t(as.data.frame(ac)), 3)
# 
# 
# br <- lapply(bries, colMeans)
# round(t(as.data.frame(br)), 3)
# 
# lapply(bries, colMeans)

#-------------------------------------------------------------------------------
# output table


out_auc <- t(as.data.frame(lapply(aucs, colMeans))) %>% 
    as.data.frame() %>%
    mutate(across(everything(), ~ sprintf("%.3f",.x))) %>%
    tibble::rownames_to_column(var = "Data") 

openxlsx::write.xlsx(out_auc, asTable = FALSE, file = paste0("../../results/sim_auc.xlsx"))

#-

out_acc <- t(as.data.frame(lapply(accuracies, colMeans))) %>% 
    as.data.frame() %>%
    mutate(across(everything(), ~ sprintf("%.3f",.x))) %>%
    tibble::rownames_to_column(var = "Data") 

openxlsx::write.xlsx(out_acc, asTable = FALSE, file = paste0("../../results/sim_acc.xlsx"))

#-------------------------------------------------------------------------------

# generate stacked barplot 
barPlot <- function(obj, ylim = NULL, ylab = "", title = "") {
    
    plotDat <- obj %>%
        as.data.frame() %>%
        pivot_longer(cols = everything(), names_to = "method") %>%
        mutate(method2 = case_when(method == "las" ~ "Lasso",
                                   method == "xgb" ~ "XGBoost", 
                                   method == "svm" ~ "SVM", 
                                   method == "extraRF" ~ "ExTrees",
                                   method == "ferns" ~ "rFerns",
                                   method == "oblRF" ~ "ORF",
                                   method == "lcRF_glm" ~ "SLC-RF",
                                   .default = as.character(method))) %>%
        mutate(method2 = factor(method2, levels = unique(method2))) %>%
        group_by(method2) %>%
        summarise(mean = mean(value),
                  se = sd(value) / sqrt(rep)) 
    
    
    if(is.null(ylim)) {
        ylim = range(plotDat$mean) * c(0.98, 1.02)
    }
    
    p <- ggplot(plotDat, aes(x = method2, y = mean, colour = method2))+ 
        geom_pointrange(aes(ymin = mean - se, ymax = mean + se), size = 0.3) +
        theme_bw() + 
        ylab(ylab) + xlab(NULL) + 
        theme(
            legend.title = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)) + 
        scale_y_continuous(limits = ylim) +
        scale_color_brewer(palette = "Set1") + 
        ggtitle(title)
    
    return(p)
}

#-------------------------------------------------------------------------------
# Compare sample size
# AUC

auc <- aucs[["s_100_0.5_FALSE"]]
p1 <- barPlot(auc, ylim = c(0.55, 0.65), ylab = "AUC", title = "A, N = 100")

auc <- aucs[["s_400_0.5_FALSE"]]
p2 <- barPlot(auc, ylim = c(0.60, 0.70), ylab = "AUC", title = "B, N = 400")

# g <- ggpubr::ggarrange(p1, p2, nrow = 1)
# ggsave(file = "../../../results/sim_auc_N.tiff", plot = g, width = 6, height = 3, units = "in", dpi = 600,
#        compression = "lzw")

# ACC

acc <- accs[["s_100_0.5_FALSE"]]
p3 <- barPlot(acc, ylim = c(0.55, 0.60), ylab = "Accuracy", title = "C, N = 100")

acc <- accs[["s_400_0.5_FALSE"]]
p4 <- barPlot(acc, ylim = c(0.55, 0.65), ylab = "Accuracy", title = "D, N = 400")

g <- ggpubr::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
ggsave(file = "../../results/sim_N.tiff", plot = g, width = 6, height = 5.5, units = "in", dpi = 600,
       compression = "lzw")


#-------------------------------------------------------------------------------
# linear vs non-linear
# AUC
auc <- aucs[["s_100_0.5_FALSE"]]
p1 <- barPlot(auc, ylim = c(0.575, 0.650), ylab = "AUC", title = "A, all linear effects")

auc <- aucs[["s_100_0.5_TRUE"]]
p2 <- barPlot(auc, ylim = c(0.660, 0.715), ylab = "AUC", title = "B, with non-linear effects")

# g <- ggpubr::ggarrange(p1, p2, nrow = 1)
# ggsave(file = "../../../results/sim_auc_Effects.tiff", plot = g, width = 6, height = 3, units = "in", dpi = 600,
#        compression = "lzw")

# ACC
acc <- accs[["s_100_0.5_FALSE"]]
p3 <- barPlot(acc, ylim = c(0.55, 0.60), ylab = "Accuracy", title = "C, all linear effects")

acc <- accs[["s_100_0.5_TRUE"]]
p4 <- barPlot(acc, ylim = c(.60, 0.675), ylab = "Accuracy", title = "D, with non-linear effects")

g <- ggpubr::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
ggsave(file = "../../results/sim_Effects.tiff", plot = g, width = 6, height = 5.5, units = "in", dpi = 600,
       compression = "lzw")


#-------------------------------------------------------------------------------
# Correlations

auc <- aucs[["s_100_0.3_FALSE"]]
p1 <- barPlot(auc, title = expression(paste("A, ", rho, "= 0.3")), ylab = "AUC", ylim = c(0.560, 0.630))

auc <- aucs[["s_100_0.5_FALSE"]]
p2 <- barPlot(auc, title = expression(paste("B, ", rho, "= 0.5")), ylab = "AUC", ylim = c(0.575, 0.650))

auc <- aucs[["s_100_0.7_FALSE"]]
p3 <- barPlot(auc, title = expression(paste("C, ", rho, "= 0.7")), ylab = "AUC", ylim = c(0.630, 0.700))

# g <- ggpubr::ggarrange(p1, p2, p3, nrow = 1)
# ggsave(file = "../../../results/sim_auc_Correlations.tiff", plot = g, width = 9, height = 3, units = "in", dpi = 600,
#        compression = "lzw")
# 
# 

acc <- accs[["s_100_0.3_FALSE"]]
p4 <- barPlot(acc, title = expression(paste("A, ", rho, "= 0.3")),  ylab = "Accuracy", ylim = c(0.500, 0.600))

acc <- accs[["s_100_0.5_FALSE"]]
p5 <- barPlot(acc, title = expression(paste("B, ", rho, "= 0.5")),  ylab = "Accuracy", ylim = c(0.550, 0.600))

acc <- accs[["s_100_0.7_FALSE"]]
p6 <- barPlot(acc, title = expression(paste("C, ", rho, "= 0.7")),  ylab = "Accuracy", ylim = c(0.550, 0.630))

g <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
ggsave(file = "../../results/sim_Correlations.tiff", plot = g, width = 9, height = 5.5, units = "in", dpi = 600,
       compression = "lzw")


#-------------------------------------------------------------------------------
# decomposition of bias and variance
#-------------------------------------------------------------------------------

deError <- function(haty, y) {
    
    ind <- abs(y - haty) < 0.5
    bias <- ifelse(ind, 0, 1)
    variance <- ifelse(ind, abs(y - haty), abs(y - haty) - 1)
    
    list(bias = bias,
         variance = abs(variance))
}

biasList <- list()
VarList <- list()

scenario = 1

for(rho_w in c(0.3, 0.5, 0.7)) {
    
    for(n_obs in c(100, 400)) {
        
        for(nonLinear in c(TRUE, FALSE)) {
            
            message(scenario, " in a total of 12 scenarios! \n",
                    "rho_w = ", rho_w, "\n",
                    "n_obs = ", n_obs, "\n",
                    "nonLinear = ", nonLinear, "\n")
            
            scenario <- scenario + 1
            
            for(times in 1:rep) {
                
                nam <- paste0(c(n_obs, rho_w, nonLinear, times), collapse = "_")
                load(paste0("../../output/sim/sim_", nam, ".RData"))
                
                # calculate AUC
                haty <- lapply(mods$pred, `[[`, "mod")
                y <- mods$y
                
                tmpError <- lapply(haty, deError, y = y)
                tmpBias  <- lapply(tmpError, `[[`, "bias")
                tmpVariance <- lapply(tmpError, `[[`, "variance")
                
                if (times == 1){
                    bias <- tmpBias
                    variance <- tmpVariance
                } else {
                    bias <- mapply(cbind, tmpBias, bias, SIMPLIFY = FALSE)
                    variance <- mapply(cbind, tmpVariance, variance, SIMPLIFY = FALSE)
                }
            }
            
            sumbias <- c()
            sumUnbiasedVar <- c()
            sumBiasedVar <- c()
            
            for(i in seq_along(bias)) {
                
                b <- bias[[i]]
                bs <- ifelse(rowMeans(b) > 0.5, 1, 0)
                sumbias <- c(sumbias, mean(bs))
                
                v <- variance[[i]]
                
                sumUnbiasedVar <- c(sumUnbiasedVar, mean(rowMeans(v * abs(1 - b))))
                sumBiasedVar <- c(sumBiasedVar, mean(rowMeans(v * b)))
                
            }
            
            names(sumbias) <- names(bias)
            names(sumUnbiasedVar) <- names(bias)
            
            names(sumBiasedVar) <- names(bias)
            netVar <- sumUnbiasedVar - sumBiasedVar
            
            nam <- paste0(c(n_obs, rho_w, nonLinear), collapse = "_")
            biasList[[nam]] <- sumbias
            VarList[[nam]] <- netVar
        }
    }
}


bias <- t(as.data.frame(biasList))

plotData <- bias %>%
    as.data.frame() %>% 
    rownames_to_column(var = "scenario")%>%
    pivot_longer(c("RF", "lcRF_glm"), names_to = "method") %>%
    mutate(scenario = gsub("X", "", scenario),
           scenario = gsub("_", ",", scenario),
           scenario = sub("(.*),1$", "\\1,Lin", scenario),
           scenario = sub("(.*),0$", "\\1,Non-Lin", scenario),
           scenario = paste0("(", scenario, ")"),
           method = case_when(method == "RF" ~ "Random Forest",
                              method == "lcRF_glm" ~ "SLC-RF")) %>%
    mutate(ifLinear = str_detect(scenario, "Non")) %>%
    arrange(-desc(ifLinear)) %>%
    mutate(scenario = factor(scenario, levels = unique(scenario)))

plotData_wide <- plotData %>% 
    pivot_wider(names_from = method, values_from = value)


# 3. Create the plot
biasPlot <- ggplot() +
    geom_point(data = plotData, aes(x = scenario, y = value, fill = method),
               shape = 21, size = 3, stroke = 1) +
    geom_segment(data = plotData_wide,
                 aes(x = scenario, xend = scenario, 
                     y = `SLC-RF`, yend = `Random Forest`,
                     color = ifelse(`Random Forest` > `SLC-RF`, "red", "blue")),
                 linewidth = 1) +
    # Use the actual color names ("red", "blue") directly for the segments.
    scale_color_identity() +
    # Manually set the fill colors for the points to match the Python plot.
    scale_fill_manual(values = c("Random Forest" = "darkorange", "SLC-RF" = "steelblue")) +
    scale_y_continuous(limits = c(0.28, 0.43)) + 
    
    # Add labels and a title.
    labs(
        title = "A, bias",
        x = element_blank(),
        y = "Bias",
        fill = "Method" # This will be the title of the legend for the points.
    ) + 
    # Apply a clean theme and customize text elements.
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom")


# vars plot


vars <- t(as.data.frame(VarList))

plotData <- vars %>%
    as.data.frame() %>% 
    rownames_to_column(var = "scenario")%>%
    pivot_longer(c("RF", "lcRF_glm"), names_to = "method") %>%
    mutate(scenario = gsub("X", "", scenario),
           scenario = gsub("_", ",", scenario),
           scenario = sub("(.*),1$", "\\1,Lin", scenario),
           scenario = sub("(.*),0$", "\\1,Non-Lin", scenario),
           scenario = paste0("(", scenario, ")"),
           method = case_when(method == "RF" ~ "Random Forest",
                              method == "lcRF_glm" ~ "SLC-RF")) %>%
    mutate(ifLinear = str_detect(scenario, "Non")) %>%
    arrange(-desc(ifLinear)) %>%
    mutate(scenario = factor(scenario, levels = unique(scenario)))

plotData_wide <- plotData %>% 
    pivot_wider(names_from = method, values_from = value)


# 3. Create the plot
varsPlot <- ggplot() +
    geom_point(data = plotData, aes(x = scenario, y = value, fill = method),
               shape = 21, size = 3, stroke = 1) +
    geom_segment(data = plotData_wide,
                 aes(x = scenario, xend = scenario, 
                     y = `SLC-RF`, yend = `Random Forest`,
                     color = ifelse(`Random Forest` > `SLC-RF`, "red", "blue")),
                 linewidth = 1) +
    # Use the actual color names ("red", "blue") directly for the segments.
    scale_color_identity() +
    # Manually set the fill colors for the points to match the Python plot.
    scale_fill_manual(values = c("Random Forest" = "darkorange", "SLC-RF" = "steelblue")) +
    scale_y_continuous(limits = c(0, 0.15)) + 
    # Add labels and a title.
    labs(
        title = "B, net variance",
        x = element_blank(),
        y = "Variance",
        fill = "Method" # This will be the title of the legend for the points.
    ) + 
    # Apply a clean theme and customize text elements.
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom")

g <- ggpubr::ggarrange(biasPlot, varsPlot, nrow = 2, common.legend  = TRUE, legend = "bottom")
ggsave(file = "../../../results/sim_Bias_vars.tiff", plot = g, width = 7, height = 7, units = "in", dpi = 600, 
       compression = "lzw")











































































