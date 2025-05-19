# bootstraped kaplan meier estimator 

###

# bootstrapped (n = 1000) time regression curves at 3 years 
# the ROC curves were averaged and interpolated 
# A curve was given for an estimate with just the risk factors 
# another with risk factors and the proteins added
# p.values were caluclated with timeROC compare function 
# protein 
# This analysis was performed for all univariatly significant proteins (FDR < 0.05)
#
# Required input:
# - dataInterestMACE
# - univariateResults
#
# Resulting output:
# - plot with top n protein timeROC AUC mean and 95% CI interval at time point indicated
# - results dataframe with timeROC calculates AUC for each univariately significant proteins

###


### Required libraries 
library(timeROC)
library(survival)
library(pROC)
library(splitTools)
library(ggplot2) 

# approach insprired by https://datascienceplus.com/time-dependent-roc-for-survival-prediction-models-in-r/?

time_points <- c(2.99)
n <- 34

### calculate timeROC - for risk, risk + protein and only single proteins
# only for univariately significant ones

data <- dataInterestMACE
proteins <- univariateProteins
results <- data.frame(matrix(data = NA, nrow = length(proteins)+1, ncol = 8))
colnames(results) <- c("protein: AUC",
                       "protein: AUC lower CI",
                       "protein: AUC upper CI",
                       "protein + risk: AUC",
                       "protein + risk: AUC lower CI",
                       "protein + risk: AUC upper CI",
                       "p-value",
                       "adjusted p-value")
rownames(results) <- c("risk.only", proteins)

model.risk <- coxph(Surv(timeMACE, eventMACE) ~ GFR_MDRD + creat + HDL + LDL + DM.composite + hb + Age + Sex, data = data)
pred.risk <- predict(model.risk, type = "lp")
roc.risk <- timeROC(T = data$timeMACE,
                    delta = data$eventMACE,
                    marker = pred.risk,
                    cause = 1,
                    weighting = "marginal",
                    times = time_points,
                    ROC = TRUE,
                    iid = TRUE)
results["risk.only", c("protein + risk: AUC", "protein + risk: AUC lower CI", "protein + risk: AUC upper CI")] <- c(round(roc.risk$AUC[2], digits = 3),
                                                                                                                  confint(roc.risk, level = 0.95, n.sim = 1000)$CI_AUC[1:2]/100)

for (protein in proteins){
  model.risk.protein <- coxph(as.formula(paste("Surv(timeMACE, eventMACE) ~", protein, "+ GFR_MDRD + creat + HDL + LDL + DM.composite + hb + Age + Sex")), data = data)
  pred.risk.protein <- predict(model.risk.protein, type = "lp")
  roc.risk.protein <- timeROC(T = data$timeMACE,
                 delta = data$eventMACE,
                 marker = pred.risk.protein,
                 cause = 1,
                 weighting = "marginal",
                 times = time_points,
                 ROC = TRUE,
                 iid = TRUE) 
  confint(roc.risk.protein, level = 0.95, n.sim = 1000)$CI_AUC[1:2]
  results[protein, c("protein + risk: AUC", "protein + risk: AUC lower CI", "protein + risk: AUC upper CI")] <- c(round(roc.risk.protein$AUC[2], digits = 3),
                                                                                             confint(roc.risk.protein, level = 0.95, n.sim = 1000)$CI_AUC[1:2]/100)

  results[protein, "p-value"] <- compare(roc.risk.protein, roc.risk)$p_values_AUC[2]  # calculate significant differences in risk.only and riskfactors+single protein models 

  model.protein <- coxph(as.formula(paste("Surv(timeMACE, eventMACE) ~", protein)), data = data)
  pred.protein <- predict(model.protein, type = "lp")
  roc.protein <- timeROC(T = data$timeMACE,
                 delta = data$eventMACE,
                 marker = pred.protein,
                 cause = 1,
                 weighting = "marginal",
                 times = time_points,
                 ROC = TRUE,
                 iid = TRUE)

  results[protein, c("protein: AUC", "protein: AUC lower CI", "protein: AUC upper CI")] <- c(round(roc.protein$AUC[2], digits = 3),
                                                                                             confint(roc.protein, level = 0.95, n.sim = 1000)$CI_AUC[1:2]/100)
}
 
### plot results (mean AUC and 95% CI for top n protein measurements)
results$protein <- c("risk.only", univariateResults[proteins, "Assay"])
results <- results[order(results$`protein + risk: AUC`, decreasing = TRUE), ]
color <- c(rep("#9B8EDA", n), "#6699FF")
 
# create dataframe used as subset to visualize in plot
df <- results[1:n, ]
df[n+1, ] <- results["risk.only", ]

df$protein<- factor(df$protein, levels = make.unique(as.vector(df$protein)))
df$n <- 1:nrow(df) # preserve original order 

plotAUC <- ggplot(df) +
  # First set (shift slightly left)
  geom_point(aes(x = n - 0.1, y = `protein + risk: AUC`),
             size = 2, alpha = 0.75, color = color) +
  geom_errorbar(aes(x = n - 0.1,
                    ymin = `protein + risk: AUC lower CI`,
                    ymax = `protein + risk: AUC upper CI`),
                width = 0.1, color = color) +

  # Second set (centered or right)
  geom_point(aes(x = n + 0.1 , y = `protein: AUC`),
             size = 2, alpha = 0.75, color = "#52B6A8") +
  geom_errorbar(aes(x = n + 0.1,
                    ymin = `protein: AUC lower CI`,
                    ymax = `protein: AUC upper CI`),
                width = 0.1, color = "#52B6A8") + 
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.05)) +
  scale_x_continuous(breaks = as.numeric(df$n), labels = df$protein) +
  labs(title = "AUC with 95% Confidence Intervals by Protein",
       x = "Protein", y = "AUC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

plotAUC

 

 

 
