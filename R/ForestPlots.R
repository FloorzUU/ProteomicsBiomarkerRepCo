### creating forest plot of top proteins of the cox regression 

### 

# In this script two different forest plots are created. 
# The first plot is one that shows the univariate and multivariate cox regression
# outcomes. For this plot a number (30) of proteins with the smallest p.values
# in the univariate cox regression. 
# Outcomes of the cox regression were depicted
# 
# The second forest plot shows a selected group of proteins -> selectionAssay
# For each of these proteins outcomes of cox regression with different confounders 
# are created. This forest highlight te effects of different confounders on the 
# hazard ratio's of the selected proteins 
#
# Required Input: 
# - dataInterestMACE
# - pValues (as part of the outcomes from the cox regression script)
#
# Resulting output: 
# - Forest plot of multi- and univariate cox regression of the most significant proteins 
# - Forest plot of a selection of proteins with outcomes of cox regression under 
# the addition of different confounders
 
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

### import required libraries 
library(forestplot)       # version 3.1.3
library(RColorBrewer)     # version 1.1-3

### forest plot 
numberProteins <- 30
rownames(pValues) <- colnames(dataInterestMACE)[grepl("^OID", colnames(dataInterestMACE))]
topProteins <- rownames(pValues[order(pValues$Multivariate)[1:numberProteins], ])
resultsForest <- as.data.frame(matrix(nrow = 2*numberProteins, ncol = 7))
colnames(resultsForest) <- c("labeltext", "Assay", "mean", "lower", "upper", "group", "p.value")

i = 1
for (protein in topProteins){
  for (j in c(1, 2)){
    if (j == 1){
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ LDL + HDL + Age + Sex + GFR_MDRD + creat + hb + DM.composite", sep = "")), data = dataInterestMACE)
      status <- "multivariate"
    } else{ 
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, sep = "")), data = dataInterestMACE)
      status <- "univariate"
    }
    x <- summary(x)
    p.value <- as.numeric(signif(x$coef[1, 5], digits = 4))
    HR <- as.numeric(signif(x$coef[1, 2], digits = 3))                   
    HR.confint.lower <- as.numeric(signif(x$conf.int[1,"lower .95"], 2))
    HR.confint.upper <- as.numeric(signif(x$conf.int[1,"upper .95"], 2))
    Assay <- data[data$OlinkID %in% rownames(x$coef), "Assay"][[1,1]]
    
    resultsForest[i,] <- c(protein, Assay, HR, HR.confint.lower, HR.confint.upper, status, p.value)
    i = i + 1
  }
}
# prepare plot for use of forest plot package
resultsForest$mean <- as.numeric(resultsForest$mean)
resultsForest$lower <- as.numeric(resultsForest$lower)
resultsForest$upper <- as.numeric(resultsForest$upper)
resultsForest$pvalue.multi <- rep(resultsForest$p.value[c(TRUE, FALSE)], each = 2)
resultsForest$pvalue.uni <- rep(resultsForest$p.value[c(FALSE, TRUE)], each = 2)
resultsForest$pvalue.uni.adjusted <- rep(p.adjust(pValues$Univariate)[topProteins], each = 2)

# add LOD 
for (protein in resultsForest$labeltext){
  resultsForest[resultsForest$labeltext == protein, "LOD"] <- data[data$OlinkID == protein, "MissingFreq"][[1]][1]
} # I end up not printing LOD anymore 

# order on p.value multi
resultsForest <- resultsForest[order(as.numeric(resultsForest$pvalue.multi)), ]

p <- resultsForest |>
  group_by(group) |>
  forestplot(labeltext = c(Assay, pvalue.uni, pvalue.uni.adjusted, pvalue.multi),
             legend_args = fpLegend(pos = list(x = .2, y = 0.98)),
             clip = c(0.1, 4.5),
             ci.vertices = TRUE,
             ci.vertices.height = 0.1,
             boxsize = .25,
             xlab = "Hazard ratio", 
             cex = 2, 
             lineheight = "auto", 
             logodds = TRUE,
             axes = gpar(cex = 0.6), 
             txt_gp = fpTxtGp(cex=0.75)) |> 
  fp_add_lines("steelblue") |> 
  fp_add_header(Assay = c("", "Protein"),
                pvalue.uni = c("p.value", "univariate"), 
                pvalue.uni.adjusted = c("adjusted p.value", "univariate"),
                pvalue.multi = c("p.value", "multivariate")) |> 
  fp_set_style(box = c("#98C7D4", "#002060") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) |> 
  fp_decorate_graph(grid = structure(c(1), 
                                     gp = gpar(lty = 2, col = "black")), graph.pos = 5) |> 
  fp_set_zebra_style("#f9f9f9")
p[["align"]] <- c('l', 'l', 'l', 'l', 'l')
p[["txt_gp"]][["xlab"]][["cex"]] <- 0.7
p[["txt_gp"]][["ticks"]][["cex"]] <- 0.6

pdf(file = paste(resultsFolder, "ForestPlot.pdf", sep = ""), width = 9, height = 0.2*numberProteins)
plot(p)
dev.off()

rm(i, j, Assay, HR, HR.confint.lower, HR.confint.upper, resultsForest, x, status, 
   topProteins, p, numberProteins)


### forest plot infividual confounders
n <- 10
selectionAssay <- c("CHCHD10", "IFI30", "AREG", "FABP5", "SPRR3", 
                    "OSM", "PTGES2", "FGFBP2", "WFDC2", "GDF15", "CXCL8", "NT-proBNP")
selectionOlinkID <- c()
for(sel in selectionAssay){
  selectionOlinkID <- c(selectionOlinkID, data[grepl(paste("^",sel,"$", sep = ""), data$Assay), "OlinkID"][[1]][1])
}
rownames(pValues) <- colnames(dataInterestMACE)[grepl("^OID", colnames(dataInterestMACE))]
resultsForest <- as.data.frame(matrix(nrow = n*length(selectionOlinkID), ncol = 7))
colors <- brewer.pal(n, "Spectral")

colnames(resultsForest) <- c("labeltext", "Assay", "mean", "lower", "upper", "group", "p.value")

i = 1
for (protein in selectionOlinkID){
  for (j in 1:n){
    if (j == 1){
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ LDL + HDL + Age + Sex + GFR_MDRD + creat + hb + DM.composite", sep = "")), data = dataInterestMACE)
      status <- "Multivariate"
    } else if (j == 2){
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ creat", sep = "")), data = dataInterestMACE)
      status <- "Creatinine"
    } else if (j == 3){
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ GFR_MDRD", sep = "")), data = dataInterestMACE)
      status <- "eGFR"
    } else if (j == 10) { 
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, sep = "")), data = dataInterestMACE)
      status <- "Univariate"
    } else if (j == 7){
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ LDL", sep = "")), data = dataInterestMACE)
      status <- "LDL"
    } else if (j == 8){
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ HDL", sep = "")), data = dataInterestMACE)
      status <- "HDL"
    } else if (j == 4) { 
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ hb", sep = "")), data = dataInterestMACE)
      status <- "Hemoglobin"
    } else if (j == 5) { 
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ Age", sep = "")), data = dataInterestMACE)
      status <- "Age"
    } else if (j == 9) { 
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ Sex", sep = "")), data = dataInterestMACE)
      status <- "Sex"
    } else if (j == 6) { 
      x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "+ DM.composite", sep = "")), data = dataInterestMACE)
      status <- "Diabetus mellitus"
    }
    
    x <- summary(x)
    p.value <- as.numeric(signif(x$coef[1, 5], digits = 4))
    HR <- as.numeric(signif(x$coef[1, 2], digits = 3))                   
    HR.confint.lower <- as.numeric(signif(x$conf.int[1,"lower .95"], 2))
    HR.confint.upper <- as.numeric(signif(x$conf.int[1,"upper .95"], 2))
    Assay <- data[data$OlinkID %in% rownames(x$coef), "Assay"][[1,1]]
    
    resultsForest[i,] <- c(protein, Assay, HR, HR.confint.lower, HR.confint.upper, status, p.value)
    i = i + 1
  }
}

# prepare plot for use of forest plot package
resultsForest$mean <- as.numeric(resultsForest$mean)
resultsForest$lower <- as.numeric(resultsForest$lower)
resultsForest$upper <- as.numeric(resultsForest$upper)
resultsForest$pvalue.multi <- rep(resultsForest$p.value[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)], each = n)
resultsForest$pvalue.uni <- rep(resultsForest$p.value[c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)], each = n)
resultsForest$creat <- rep(resultsForest$p.value[c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)], each = n)
resultsForest$LDL <- rep(resultsForest$p.value[c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)], each = n)
resultsForest$HDL <- rep(resultsForest$p.value[c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)], each = n)
resultsForest$DM.composite <- rep(resultsForest$p.value[c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)], each = n)
resultsForest$hb <- rep(resultsForest$p.value[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)], each = n)
resultsForest$GFR_MDRD <- rep(resultsForest$p.value[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)], each = n)
resultsForest$Age <- rep(resultsForest$p.value[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)], each = n)
resultsForest$Sex <- rep(resultsForest$p.value[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)], each = n)
resultsForest$pvalue.uni.adjusted <- rep(pValues[selectionOlinkID, "Univariate"], each = n)

               
p <- resultsForest |>
  group_by(group) |>
  forestplot(labeltext = c(Assay, pvalue.uni, pvalue.uni.adjusted, pvalue.multi),
             legend_args = fpLegend(pos = list("topright"),),
             clip = c(0.1, 4.5),
             ci.vertices = TRUE,
             ci.vertices.height = 0.038, 
             boxsize = .06,
             xlab = "Hazard ratio", 
             cex = 2, 
             lineheight = "auto", 
             axes = gpar(cex = 0.6), 
             txt_gp = fpTxtGp(cex=0.75), 
             logodds = TRUE) |> 
  fp_add_lines("steelblue") |> 
  fp_add_header(Assay = c("", "Protein"),
                pvalue.uni = c("p.value", "univariate"), 
                pvalue.uni.adjusted = c("adjusted p.value", "univariate"),
                pvalue.multi = c("p.value", "multivariate")) |> 
  fp_set_style(box = colors |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) |> 
  fp_decorate_graph(grid = structure(c(1), 
                                     gp = gpar(lty = 6, col = "black")), graph.pos = 5) |> 
  fp_set_zebra_style("#f9f9f9")
p[["align"]] <- c('l', 'l', 'l', 'l', 'l')
p[["txt_gp"]][["xlab"]][["cex"]] <- 0.7
p[["txt_gp"]][["ticks"]][["cex"]] <- 0.6

pdf(file = paste(resultsFolder, "ForestPlot_extra.pdf", sep = ""), width = 9, height = 16)
plot(p)
dev.off()
