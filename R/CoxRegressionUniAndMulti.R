# Cox regressions univariate and with different confounders 

### 

# In this script univariate cox regressions were performed for all proteins
# as well as cox regressions with different combinations of confounders 
# the outcomes are depicted univariatly in a vulcano plot 
# the outcomes (p.values, beta and hazard ratios) are depicted multivariatly 
# in heatmaps one with top proteins selected and one for all proteins
#  
# Required Input: 
# - dataInterestMACE
#
# Resulting output:
# - Vulcano plot with outcomes in the Cox-regression 
# - heatmap beta values in different combinations of confounders (selection proteins and all)
# - heatmap adjusted p.values in different combinations of confounders (selection proteins and all)
# - heatmap hazard ratio's in different combinations of confounders (selection proteins and all)

###

### setting computer path 
if (Sys.info()[4][[1]] == "TURBO-PC-WESTER"){
  pathStart <- "C:/Users/Bas/Documents/FloorsFiles/MASTER/MajorInternship/"
} else if (Sys.info()[4][[1]] == "DESKTOP-2DADF21"){
  pathStart <- "C:/Users/Floor van der Zalm/Documents/MASTER/MajorInternship/"
} else {
  warning ("Please adjust file paths before executing")
}
# file locations
resultsFolder <- paste0(pathStart, "Pipeline1/Results/pipelineLargeDataset/predictionPreperation/MACE/")


### import relevant libraries 
library(survival)         # version 3.3-1
library(pheatmap)         # version 1.0.12


### Perform univariate cox regression for all proteins 
proteins <- colnames(dataInterestMACE)[grepl("^OID", colnames(dataInterestMACE))]
univariateFormulas <- sapply(proteins, function(x) as.formula(paste('Surv(timeMACE, eventMACE)~', x)))
univariateModels <- lapply(univariateFormulas, function(x){coxph(x, data = dataInterestMACE)})
univariateResults <- lapply(univariateModels, function(x){ 
  x <- summary(x)
  p.value <- signif(x$coef[5], digits = 4)
  wald.test <- signif(x$wald["test"], digits = 4)
  beta <- signif(x$coef[1], digits = 3);                 # coefficient beta
  HR <- signif(x$coef[2], digits = 3);                   # exp(beta)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  Assay <- data[data$OlinkID %in% rownames(x$coef), "Assay"][[1,1]]
  results <- c(Assay, beta, HR, wald.test, p.value, "NA")
  names(results) <- c("Assay", "beta", "HR (95% CI for HR)", "wald.test", "p.value", "adjusted.p.value")
  return(results)
})
univariateResults <- t(as.data.frame(univariateResults))
univariateResults[ ,"adjusted.p.value"] <- signif(p.adjust(as.numeric(univariateResults[ ,"p.value"]), method = "BH"), 4)
univariateResults <- univariateResults[order(univariateResults[ ,"p.value"]), ] 
univariateProteins <- names(which(univariateResults[ ,"adjusted.p.value"] < 0.05))

### plotting vulcano plot univariate Cox regression 
colors <- c()
univariateResults <- as.data.frame(univariateResults)
colors[univariateResults[ ,"adjusted.p.value"] < 0.05] <- "red" 
colors[univariateResults[ ,"adjusted.p.value"] > 0.05] <- "#00C7E1"
pdf(file = paste(resultsFolder, "volcanoPlot.pdf", sep = ""), 
    width = 5, height = 6.5)
plot(univariateResults$beta, -log(as.numeric(univariateResults$p.value), base = 10), 
     xlab = "Beta", ylab = "-log10(p.values)", pch = 20, col = colors, 
     cex = 0.7) 
abline(h = -log10(0.05), col = "grey", lty = 2)
legend("topleft", legend = c("Significant", "Not Significant"), 
       col = c("red", "#00C7E1"), pch = 20, cex = 0.7, bg = "transparent", bty = "o")
dev.off()


### Perform cox regression with different confounders for all proteins 
# the different of combinations to run for confounding variables
confounders <- c("", "+ LDL", "+ HDL", "+ Age", "+ Sex", "+ creat", "+ hb", 
                 "+ DM.composite", "+ GFR_MDRD", 
                 "+ LDL + HDL + Age + Sex + GFR_MDRD + creat + hb + DM.composite")

# set up dataframes pValues and Beta values 
pValues <- as.data.frame(matrix(ncol = length(confounders), nrow = length(proteins)))
rownames(pValues) <- proteins
colnames(pValues) <- c("Univariate", "LDL", "HDL", "Age", "Sex", "Creatinin", "Hemoglobin", 
                       "Diabetes mellitus", "eGFR", 
                       "Multivariate")
betaValues <- pValues
HRValues <- pValues

i <- 1
for (confounder in confounders){
  multivariableFormulas <- sapply(colnames(dataInterestMACE)[grepl("^OID", colnames(dataInterestMACE))],
                                  function(x) as.formula(paste('Surv(timeMACE, eventMACE)~', x, confounder)))
  multivariableModels <- lapply(multivariableFormulas, function(x){coxph(x, data = dataInterestMACE)})
  multivariableResults <- lapply(multivariableModels, function(x){ 
    x <- summary(x)
    p.value <- signif(x$coef[1,5], digits = 4)
    wald.test <- signif(x$wald["test"], digits = 4)
    beta <- signif(x$coef[1, 1], digits = 3)                  # coefficient beta
    HR <- signif(x$coef[1, 2], digits = 3)                    # exp(beta)
    HR.confint.lower <- signif(x$conf.int[1, "lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[1, "upper .95"], 2)
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    Assay <- data[data$OlinkID %in% rownames(x$coef), "Assay"][[1,1]]
    HRshort <- signif(x$coef[1, 2], digits = 3)
    results <- c(Assay, beta, HR, HRshort, wald.test, p.value, "NA")
    names(results) <- c("Assay", "beta", "HR (95% CI for HR)", "HR", "wald.test", "p.value", "adjusted.p.value")
    return(results)
  })
  multivariableResults <- t(as.data.frame(multivariableResults))
  multivariableResults[ ,"adjusted.p.value"] <- signif(p.adjust(as.numeric(multivariableResults[ ,"p.value"]), method = "BH"), 4)
  
  pValues[ ,colnames(pValues)[i]] <- signif(p.adjust(as.numeric(multivariableResults[ ,"p.value"]), method = "BH"), 4)
  betaValues[ ,colnames(betaValues)[i]] <- multivariableResults[ ,"beta"]
  HRValues[ ,colnames(HRValues[i])] <- multivariableResults[ ,"HR"]
  
  multivariableResults <- multivariableResults[order(multivariableResults[ ,"p.value"]), ] 
  i <- i + 1
} 

assayNames <- character(length(proteins))
for (i in 1:length(proteins)) {
  assayNames[i] <- data[data$OlinkID == proteins[i], "Assay"][[1]][1]
}
append_numbers <- function(x) {
  ave(x, x, FUN = function(y) if(length(y) > 2) paste0(y, "_", seq_along(y)) else y)
}

# heatmap of protein values
x <- list()
rownames(pValues) <- append_numbers(assayNames)
pValuesSel <- pValues[apply(pValues, 1, FUN = min) < 0.05, ]
breaks <- seq(floor(min(-log(pValues, base = 10))), ceiling(max(-log(pValues, base = 10))), length.out = 101)
pdf(file = paste(resultsFolder, "pvaluesConfounders.pdf", sep = ""),width = 5, height = 9)
x[[1]] <- pheatmap(-log(pValues, base = 10), annotation_names_row = FALSE, show_rownames = FALSE, 
                   cutree_rows = 5, legend = TRUE, annotation_legend = TRUE, 
                   breaks = breaks, angle_col = "45")
print(x[[1]])
dev.off()

# zoomed in map selestion of proteins
selection <- order(pValuesSel$Multivariate)[1:30]
pValuesSelZoomed <- pValuesSel[selection,!colnames(pValues) %in% c("Multivariate")]
namesZoomed <- sub("^[^ ]+ ", "", rownames(pValuesSel[selection, ]))
namesZoomed <- rownames(pValuesSel[selection, ])
rownames(pValuesSelZoomed) <- namesZoomed
pValuesSelZoomed <- pValuesSelZoomed[ ,colnames(pValuesSelZoomed) %in% 
                                        c("Univariate","LDL","HDL","Age","Sex","Creatinin","Diabetes mellitus","eGFR", "Hemoglobin","Multivariate")]
asterisk_matrix <- pValuesSelZoomed
asterisk_matrix[pValuesSelZoomed < 0.05 & pValuesSelZoomed >= 0.01] <- "*"
asterisk_matrix[pValuesSelZoomed < 0.01 & pValuesSelZoomed >= 0.001] <- "**"
asterisk_matrix[pValuesSelZoomed < 0.001] <- "***"
asterisk_matrix[pValuesSelZoomed > 0.05 & pValuesSelZoomed < 0.1] <- "~"
asterisk_matrix[pValuesSelZoomed > 0.1 ] <- ""

pdf(file = paste(resultsFolder, "pvaluesConfoundersZoomed.pdf", sep = ""), 
    width = 6.5, height = 6.5)
x[[2]] <- pheatmap(-log(pValuesSelZoomed, base = 10), breaks = breaks, 
                   display_numbers = asterisk_matrix, legend = TRUE, cluster_rows = T,
                   fontsize_row = 9, fontsize_col = 9, angle_col = "45")
print(x[[2]])
dev.off()

# heatmap of beta values 
mat_colors <- colorRampPalette(c("navy", "white", "red"))(100)
betaValues <- as.data.frame(sapply(betaValues, as.numeric))
rownames(betaValues) <- append_numbers(assayNames)
betaValuesSel <- betaValues[rownames(pValuesSel), ]
breaks <- seq(-max(abs(betaValuesSel)), max(abs(betaValuesSel)), length.out = 101)
pdf(file = paste(resultsFolder, "betaConfounders.pdf", sep = ""), 
    width = 6.5, height = 6.5)
x[[3]] <- pheatmap(betaValuesSel, show_rownames = FALSE, color = mat_colors, 
                   breaks = breaks, cutree_rows = 3, scale = "none")
print(x[[3]])
dev.off()

# heatmap beta values
betaValuesSelZoomed <- betaValuesSel[selection, ]
pdf(file = paste(resultsFolder, "betaConfoundersZoomed.pdf", sep = ""), 
    width = 6.5, height = 6.5)
x[[4]] <- pheatmap(betaValuesSelZoomed, show_rownames = TRUE, 
                   cutree_rows =  2, breaks = breaks, color = mat_colors, 
                   cluster_rows = F, cluster_cols = F, 
                   fontsize_row = 10, fontsize_col = 10, angle_col = "45") 
print(x[[4]])
dev.off()

# heatmap of HR values 
mat_colors <- colorRampPalette(c("navy", "white", "red"))(100)
HRValues <- as.data.frame(sapply(HRValues, as.numeric))
rownames(HRValues) <- append_numbers(assayNames)
HRValuesSel <- HRValues[rownames(pValuesSel), ]
breaks <- c(seq(0, 1, length.out = 51), seq(1.01, max(HRValues), length.out = 50))
pdf(file = paste(resultsFolder, "HRConfounders.pdf", sep = ""), 
    width = 5, height = 9)
x[[5]] <- pheatmap(HRValues, show_rownames = FALSE, color = mat_colors, 
                   cutree_rows = 1, breaks = breaks, angle = "45")
print(x[[5]])
dev.off()

# heatmap HR values
HRValuesSelZoomed <- HRValuesSel[selection, ]
breaks <- c(seq(0, 1, length.out = 51), seq(1.01, max(HRValuesSelZoomed), length.out = 50))
pdf(file = paste(resultsFolder, "HRConfoundersZoomed.pdf", sep = ""), 
    width = 6.5, height = 6.5)
x[[6]] <- pheatmap(HRValuesSelZoomed, show_rownames = TRUE, 
                   cutree_rows =  2, color = mat_colors, 
                   cluster_rows = F, cluster_cols = F, scale = "none", 
                   fontsize_row = 10, fontsize_col = 10,
                   breaks=breaks, colsep = 1, angle = "45")
print(x[[6]])
dev.off()

rm(breaks, mat_colors, asterisk_matrix, x, HRValuesSelZoomed, HRValuesSel, 
   betaValuesSel, pValuesSel, betaValuesSelZoomed, pValuesSelZoomed)

