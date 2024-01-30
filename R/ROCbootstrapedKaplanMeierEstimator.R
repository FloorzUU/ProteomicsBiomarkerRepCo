# bootstraped kaplan meier estimator 

###

# bootstrapped (n = 1000) time regression curves at 3 years 
# the ROC curves were averaged and interpolated 
# A curve was given for an estimate with just the risk factors 
# another with risk factors and the proteins added
# p.values were calulated with a paired t.test and the roc plots are stored per 
# protein 
# This analysis was performed for all univariatly significant proteins (FDR < 0.05)
#
# Required input:
# - dataInterestMACE
# - pValues
#
# Resulting output:
# - averaged ROC curves for all included proteins 
# - a table in latex format showing the outcomes

###


### Required libraries 


library(prodlim)      # version 2023.03.31
#library(timereg)
library(pec)          # version 2023.04-12
library(timeROC)      # version 0.4
library(rms)          # version 6.7-1
library(survival)     # version 3.3-1
library(Matrix)       # version 1.6-1.1
library(polspline)    # version 1.1.24
#library(pROC)         
library(glmnet)       # version 4.1-8 
library(kableExtra)   # version 1.3.4

AUCresultsProtein <- data.frame()
proteins <- rownames(pValues[pValues$Univariate < 0.05, ])
for(protein in proteins){
  num_bootstrap <- 1000
  AUC <- c()
  AUCprotein <- c()
  TPprotein <- data.frame(matrix(data = NA, nrow = dim(dataInterestMACE)[1], ncol = num_bootstrap))
  FP <- data.frame(matrix(data = NA, nrow = dim(dataInterestMACE)[1], ncol = num_bootstrap))
  TP <- data.frame(matrix(data = NA, nrow = dim(dataInterestMACE)[1], ncol = num_bootstrap)) 
  
  #interpolated values
  x_values <- seq(0, 1, 1/300)
  TPapproxProtein <- data.frame(matrix(data = NA, nrow = length(x_values), ncol = num_bootstrap))
  TPapprox <- data.frame(matrix(data = NA, nrow = length(x_values), ncol = num_bootstrap))
  
  x <- summary(coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, sep = "")), data = dataInterestMACE))[7][[1]][1]
  
  for (i in 1:num_bootstrap){
    bootstrap_sample <- sample(1:dim(dataInterestMACE)[1], replace = TRUE)
    
    if (x > 0){
      ROCtimeCoxProtein <- timeROC(T=dataInterestMACE[bootstrap_sample,"timeMACE"], 
                                   delta=dataInterestMACE[bootstrap_sample, "eventMACE"],
                                   marker=dataInterestMACE[bootstrap_sample, protein],
                                   cause=1,
                                   weighting="marginal",
                                   other_markers= as.matrix(dataInterestMACE[bootstrap_sample,
                                                                         c("GFR_MDRD", "creat", "HDL", "LDL", "DM.composite", "hb", "Age", "Sex")]),
                                   times=c(2.99) ,ROC=TRUE, iid = TRUE) 
    } else {
      ROCtimeCoxProtein <- timeROC(T=dataInterestMACE[bootstrap_sample,"timeMACE"], 
                                   delta=dataInterestMACE[bootstrap_sample, "eventMACE"],
                                   marker=-dataInterestMACE[bootstrap_sample, protein],
                                   cause=1,
                                   weighting="marginal",
                                   other_markers= as.matrix(dataInterestMACE[bootstrap_sample,
                                                                         c("GFR_MDRD", "creat", "HDL", "LDL", "DM.composite", "hb", "Age", "Sex")]),
                                   times=c(2.99),ROC=TRUE, iid = TRUE) 
    }
    FPprotein[1:length(ROCtimeCoxProtein[["FP"]][ ,2]),i] <- ROCtimeCoxProtein[["FP"]][ ,2]
    TPprotein[1:length(ROCtimeCoxProtein[["TP"]][ ,2]),i] <- ROCtimeCoxProtein[["TP"]][ ,2]
    AUCprotein <- c(AUCprotein, ROCtimeCoxProtein[["AUC"]][2])
    
    ROCtimeCox <- timeROC(T=dataInterestMACE[bootstrap_sample,"timeMACE"], 
                          delta=dataInterestMACE[bootstrap_sample, "eventMACE"],
                          marker=-dataInterestMACE[bootstrap_sample, "GFR_MDRD"],
                          cause=1,
                          weighting="marginal",
                          other_markers=as.matrix(dataInterestMACE[bootstrap_sample,
                                                               c("creat", "HDL", "LDL", "DM.composite", "hb", "Age", "Sex")]),
                          times=c(2.9),ROC=TRUE, iid = TRUE) 
    AUC <- c(AUC, ROCtimeCox[["AUC"]][2])
    FP[1:length(ROCtimeCox[["FP"]][ ,2]),i] <- ROCtimeCox[["FP"]][ ,2]
    TP[1:length(ROCtimeCox[["TP"]][ ,2]),i] <- ROCtimeCox[["TP"]][ ,2]
    
    # interpolation 
    TPapproxProtein[ ,i] <- approx(ROCtimeCoxProtein[["FP"]][ ,2], 
                                   ROCtimeCoxProtein[["TP"]][ ,2], xout = x_values)$y
    TPapprox[ ,i] <- approx(ROCtimeCox[["FP"]][ ,2], 
                            ROCtimeCox[["TP"]][ ,2], xout = x_values)$y
  }
  # Calculate mean values at each point
  meanValuesProtein <- rowMeans(data.frame(TPapproxProtein), na.rm = FALSE)
  meanValues <- rowMeans(data.frame(TPapprox), na.rm = FALSE)
  
  # Calculate confidence intervals
  lowerCIProtein<- apply(TPapproxProtein, 1, function(x) quantile(x, 0.025))
  upperCIProtein <- apply(TPapproxProtein, 1, function(x) quantile(x, 0.975))
  lowerCI<- apply(TPapprox, 1, function(x) quantile(x, 0.025))
  upperCI <- apply(TPapprox, 1, function(x) quantile(x, 0.975))
  
  # Plot individual lines with inconsistent points
  pdf(file = paste(resultsFolder, "ROCproteinAndRisk/",protein,".pdf", sep = ""), width = 4.5, height = 4.5)
  par(mar = c(4,4,1,1))
  plot(x_values, meanValues, type = "l", col = "black", lwd = 2, xlab = "1-Specificity", ylab = "Sensitivity")
  lines(x_values, meanValuesProtein, type = "l", col = "red", lwd = 2)
  legend(x = "bottomright", 
         legend = c("Riskfactors + protein", 
                    paste("Mean AUC:", signif(mean(AUCprotein), 3)),
                    "Riskfactors", 
                    paste("Mean AUC:", signif(mean(AUC), 3))),
         col = c("red", "white", "black", "white"), lty=1, cex=0.8, bty = "n")
  title(data[data$OlinkID == protein, "Assay"][[1]][1])
  
  # Shade the area between confidence intervals
  polygon(c(x_values, rev(x_values)), c(lowerCI, rev(upperCI)), col = rgb(0,0,0, alpha = 0.2), border = NA)
  polygon(c(x_values, rev(x_values)), c(lowerCIProtein, rev(upperCIProtein)), col = rgb(0.8, 0.2, 0.2, alpha = 0.2), border = NA)
  
  AUCresultsProtein[protein, "mean"] <- mean(AUCprotein)
  AUCresultsProtein[protein, "CI-lower"] <- quantile(AUCprotein, 0.025)
  AUCresultsProtein[protein, "CI-upper"] <- quantile(AUCprotein, 0.975)
  AUCresultsProtein[protein, "Assay"] <- data[data$OlinkID == protein, "Assay"][[1]][1]
  
  AUCresultsProtein[protein, "p.value"] <- t.test(AUC, AUCprotein, alternative = "less", paired = TRUE)$p.value
  AUCresultsProtein[protein, "estimate"] <- t.test(AUC, AUCprotein, alternative = "less", paired = TRUE)$estimate
  AUCresultsProtein[protein, "CI upper"] <- t.test(AUC, AUCprotein, alternative = "less", paired = TRUE)$conf.int[1:2][2]
  
  AUCresultsProtein[protein, "mean risk"] <- mean(AUC)
  AUCresultsProtein[protein, "CI-lower risk"] <- quantile(AUC, 0.025)
  AUCresultsProtein[protein, "CI-upper risk"] <- quantile(AUC, 0.975)
  
  save(AUCresultsProtein, file = paste(resultsFolder, "ROCproteinAndRisk/AUCresultsProtein.RData", sep = ""))
  dev.off()
}

# preparing table making 
AUCtable <- AUCresultsProtein
AUCtable$mean <- paste(signif(AUCresultsProtein$mean, 4), "(", signif(AUCresultsProtein$`CI-lower`, 4), "-", signif(AUCresultsProtein$`CI-upper`, 4), ")")
AUCtable$meanW <- paste(signif(AUCresultsProtein$`mean risk`, 4), "(", signif(AUCresultsProtein$`CI-lower risk`, 4), "-", signif(AUCresultsProtein$`CI-upper risk`, 4), ")")
AUCtable$estimate <- AUCtable$estimate*-1
AUCtable$p.value.adjusted <- p.adjust(AUCtable$p.value)
AUCtable <- AUCtable[AUCtable$p.value.adjusted < 0.05, ]
kable(AUCtable[ ,c("Assay", "mean", "meanW", "estimate", "p.value", "p.value.adjusted")], format = "latex")