# table 1 

### 

# In this script a table 1 for MACE is created 
# additionally correlations between the confounders are calculated and depicted 
# this together allowed us to make desisions on which confounders to include 
# in later multivariate analysis
#
# based on the p value outcomes in table 1 and the correlation 
# between these confounders the 
#
# Required input: 
# - dataInterestMACE 
# - includedPatientData 
# - MACEdefinition 
#
# Resulting output: 
# - in viewer: html based table 1 
# - printed: latex table
# - csv files used for ML in python, one with and one without the selected proteins
# - correlation plot depicting the correlation between confounders
# - selectedRiskFactors     #confounders to be included in multivariate analysis  

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
library(dplyr)            # version 1.1.3
library(table1)           # version 1.4.3
library(tidyr)            # version 1.3.0
library(kableExtra)       # version 1.3.4
library(corrplot)         # version 0.92


### selecting the patient data to included 
orderedMACECol <- c("eventMACE", "Age", "Sex", "BMI", "Hypertension1",
                    "SmokerCurrent", "DM.composite", "creat", "GFR_MDRD", "LDL", "HDL", 
                    "totalchol", "hb", "Stroke_history", "CAD_history", "PAOD", "risk614", "Symptoms.4g", "Med.statin.derived")
patientDataMACE <- dataInterestMACE[ ,orderedMACECol]

# display data in moderate and severe grouped
patientDataMACE[patientDataMACE$Symptoms.4g %in% c("2", "3"), "Symptoms.2g"] <- "severe"
patientDataMACE[patientDataMACE$Symptoms.4g %in% c("0", "1"), "Symptoms.2g"] <- "moderate"

# exclude 
patientDataMACE <- patientDataMACE[ ,!colnames(patientDataMACE) %in% 
                                   c("Symptoms.4g", "ht", "homocys", "systolic", "diastoli")]

# Setting data types
patientDataMACE$MACE <- factor(patientDataMACE$MACE, levels = c(0, 1), labels = c("No MACE","MACE"))
patientDataMACE$Age <- as.numeric(patientDataMACE$Age)
patientDataMACE$Sex <- factor(patientDataMACE$Sex, levels = c(0, 1), labels = c("Female","Male"))
patientDataMACE$BMI <- as.numeric(patientDataMACE$BMI)
patientDataMACE$Hypertension1 <- factor(patientDataMACE$Hypertension1, levels = c("1", "0"), labels = c("yes", "no"))
patientDataMACE$creat <- as.numeric(patientDataMACE$creat)
patientDataMACE$GFR_MDRD <- as.numeric(patientDataMACE$GFR_MDRD)
patientDataMACE$SmokerCurrent <- factor(patientDataMACE$SmokerCurrent, levels = c(1, 0), labels = c("yes", "no"))
patientDataMACE$DM.composite <- factor(patientDataMACE$DM.composite, levels = c(1, 0), labels = c("yes", "no"))
patientDataMACE$PAOD <- factor(patientDataMACE$PAOD, levels = c(1, 0), labels = c("yes", "no"))
patientDataMACE$LDL <- as.numeric(patientDataMACE$LDL)
patientDataMACE$HDL <- as.numeric(patientDataMACE$HDL)
patientDataMACE$totalchol <- as.numeric(patientDataMACE$totalchol)
patientDataMACE$hb <- as.numeric(patientDataMACE$hb)
patientDataMACE$CAD_history <- factor(patientDataMACE$CAD_history, levels = c(1, 0), labels = c("yes", "no"))
patientDataMACE$Stroke_history <- factor(patientDataMACE$Stroke_history, levels = c(1, 0), labels = c("yes", "no"))
patientDataMACE$Med.statin.derived <- factor(patientDataMACE$Med.statin.derived, levels = c(1, 0), labels = c("yes", "no"))
patientDataMACE$risk614 <- factor(patientDataMACE$risk614, levels = c(1, 0), labels = c("yes", "no"))

# Setting labels  
label(patientDataMACE$Age) <- "Age (years)"
label(patientDataMACE$BMI) <- "BMI (kg/m\u00B2)"
label(patientDataMACE$Hypertension1) <- "Hypertension"
label(patientDataMACE$creat) <- "Creatinin (μmol/L)"
label(patientDataMACE$GFR_MDRD) <- "eGFR (mL/min/1.73m2)"
label(patientDataMACE$SmokerCurrent) <- "Current smoking"
label(patientDataMACE$DM.composite) <- "Diabetus mellitus"
label(patientDataMACE$Symptoms.2g) <- "Preoperative symptoms"
label(patientDataMACE$hb) <- "Hemoglobin (mmol/L)"
label(patientDataMACE$Stroke_history) <- "History of stroke"
label(patientDataMACE$CAD_history) <- "History of CAD"
label(patientDataMACE$PAOD) <- "Peripheral arterial occlusive disease"
label(patientDataMACE$totalchol) <- "Total cholesterol level (mmol/L)"
label(patientDataMACE$LDL) <- "LDL (mmol/L)"
label(patientDataMACE$HDL) <- "HDL (mmol/L)"
label(patientDataMACE$Med.statin.derived) <- "Statin use"
label(patientDataMACE$risk614) <- "Hypercholesterolemia"

caption = "Table 1: Baseline charateristics of 429 patients included grouped by MACE"
footnote = "Data are presented as n (%), mean ± standard deviation, or number patients (percentage patients). 
            P values for numeric variables are calulated with a standard two-sample t-test, 
            for categorical values a chi-square of independence is performed."

my.render.cont <- function(x, digits = 3) {
  stats_result <- stats.default(x)
  mean_str <- sprintf("%s", signif(stats_result$MEAN, digits))
  sd_str <- sprintf("%s", signif(stats_result$SD, digits))
  result <- sprintf("%s (&plusmn; %s)", mean_str, sd_str)
  return(result)
}
my.render.cat <- function(x, digits = 3) {
  if (length(unique(x)) == 2) {
    if ("yes" %in% levels(x)) {
      q <- list()
      q$yes <- stats.default(x)$yes
      result <- c(sapply(q, function(y) {
        freq_str <- sprintf("%s", signif(y$FREQ, digits))
        pct_str <- sprintf("%s", signif(y$PCT, digits))
        return(sprintf("%s (%s %%)", freq_str, pct_str))
      }))
    } else {
      result <- c("", sapply(stats.default(x), function(y) {
        freq_str <- sprintf("%s", signif(y$FREQ, digits))
        pct_str <- sprintf("%s", signif(y$PCT, digits))
        return(sprintf("%s (%s %%)", freq_str, pct_str))
      }))
    }
  } else {
    result <- c("", sapply(stats.default(x), function(y) {
      freq_str <- sprintf("%s", signif(y$FREQ, digits))
      pct_str <- sprintf("%s", signif(y$PCT, digits))
      return(sprintf("%s (%s %%)", freq_str, pct_str))
    }))
  }
  return(result)
}
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c(sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

table1(~.|eventMACE , data = patientDataMACE, caption=caption, footnote=footnote, 
       render.continuous=my.render.cont, render.categorical=my.render.cat, 
       extra.col=list(`P-value`=pvalue), topclass="Rtable1-zebra", overall=F)

tableLatex <- table1(~.|eventMACE , data = patientDataMACE, caption=caption, footnote=footnote, 
                     render.continuous=my.render.cont, render.categorical=my.render.cat, 
                     extra.col=list(`P-value`=pvalue), topclass="Rtable1-zebra", overall=F)
tableLatex <- kable(tableLatex, format = "latex", booktabs = TRUE)  #adjusments are still required 
tableLatex <- gsub("&plusmn", "pm", tableLatex)
print(tableLatex)


### checking for multicollinearity 
# correlation continous confounders with p.values < 0.01
cor_matrix <- signif(cor(dataInterestMACE[ ,c("LDL", "Age", "creat", "HDL", "totalchol", "hb", "GFR_MDRD")]), 
                     digits = 3)
colnames(cor_matrix) <- c("LDL", "Age", "Creatinin", "HDL", "Cholesterol", "Hemoblobin", "eGFR")
rownames(cor_matrix) <- colnames(cor_matrix)

# heatmap confounders correlation
pdf(file = paste(resultsFolder, "confoudersCor.pdf", sep = ""), 
    width = 6.5, height = 6.5)
corrplot(cor_matrix, method = 'shade', pch.col = 'black', tl.col = "black", 
         order = "hclust", cl.pos = 'n', addCoef.col = 'white')
dev.off()

selectedRiskFactors <- c("LDL", "HDL", "GFR_MDRD", "Age", "Sex", "DM.composite", "creat", "hb")


### saving data for ML in python 
# risk factors to be included are selected based on table1 
proteins <- colnames(dataInterestMACE)[grepl("^OID", colnames(dataInterestMACE))]
write.csv(dataInterestMACE[ ,c("eventMACE", proteins)], 
          paste0(resultsFolder, "MACEproteins.csv"), row.names=TRUE)
write.csv(dataInterestMACE[ ,c("eventMACE", selectedRiskFactors, proteins)], 
          paste0(resultsFolder, "MACEproteinsRisk.csv"), row.names=TRUE)

rm(footnote, caption, orderedMACECol, cor_matrix)


