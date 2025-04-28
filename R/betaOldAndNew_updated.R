# Cox regression MACE old data  

### requires Panel 96 data, MACEdefinition, dataInterestMACE and univariateResults


MACEdefinition <- "major"

### importing libraries
library(OlinkAnalyze)
library(haven)
library(ggplot2)
library(reshape2)
library(mice)
library(dplyr)
library(tidyr)
library(survival)        

### setting computer path 
if (Sys.info()[4][[1]] == "TURBO-PC-WESTER"){
  pathStart <- "C:/Users/Bas/Documents/FloorsFiles/MASTER/MajorInternship/"
} else if (Sys.info()[4][[1]] == "DESKTOP-2DADF21"){
  pathStart <- "C:/Users/Floor van der Zalm/Documents/MASTER/MajorInternship/"
} else {
  warning ("Please adjust file paths before executing")
}

resultsFolder <- paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/results/MACEprotein/")


# file locations
olinkLongFormatPath <- "RawData/NPXfiles/olink2022-012-080_Explore3072_EXTENDED_NPX_2023-07-10.csv"

# import data Proteins old and new
CVD2 <- read.table("data/raw/DataRocheProcessed/AEOS12_CVD2_merged.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
CVD2$SampleID <- tolower(CVD2$SampleID)

CVD3 <- read.table("data/raw/DataRocheProcessed/AEOS12_CVD3_merged.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
CVD3$SampleID <- tolower(CVD3$SampleID)

CM <- read.table("data/raw/DataRocheProcessed/AEOS1_CM.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
CM$SampleID <- tolower(CM$SampleID)

# OLINK large dataset 
dataUnproccesedExplore <- read_NPX(paste(pathStart, olinkLongFormatPath, sep = ""))

# clinical data
dataRiskFactors <- as.data.frame(read_sav("data/raw/ClinicalData/2022-21-03 AtheroExpress Database.sav"))


### setting 
overlapping_patients <- FALSE

# clean olink explore dataset 
controlFind <- function(x) !grepl("control", x)
dataExplore <- dataUnproccesedExplore[controlFind(dataUnproccesedExplore$Assay), ]  # removes control assays
numbersOnly <- function(x) !grepl("\\D", x)
dataExplore <- dataExplore[numbersOnly(dataExplore$SampleID), ]          # removes negative and positive controls 
dataExplore <- dataExplore[dataExplore$QC_Warning != "EXCLUDED", ]       # removes EXCLUDED assays 
dataExplore <- dataExplore[dataExplore$MissingFreq < 0.25, ]
patientsExcludeExplore <- c("32", "1746", "1921")          
dataExplore <- dataExplore[!dataExplore$SampleID %in% patientsExcludeExplore, ]

# patients with incomplete MACE data 
patientsMACEcomplete <- as.character(dataRiskFactors[complete.cases(dataRiskFactors[, c("epmajor.3years", "ep_major_t_3years")]), "STUDY_NUMBER"])
dataExplore <- dataExplore[dataExplore$SampleID %in% patientsMACEcomplete, ]

# totals EXPLORE
patientsExplore <- unique(dataExplore$SampleID)
proteinsExplore <- unique(dataExplore$Assay)

# clean olink old data
patientsPanel <-  intersect(intersect(unique(CM$SampleID), unique(CVD2$SampleID)), unique(CVD3$SampleID))  # same patients per plate
# remove patients due to considerations about informed consent
patientsExcludePanel <- c("1721", "2423", "2427", "2439", "2617", "3844", "2083", "3220")  
patientsPanel <- setdiff(patientsPanel, paste0("ae",patientsExcludePanel)) 
patientsPanel <- intersect(patientsPanel, paste0("ae", patientsMACEcomplete))  ### examine overlap in patients 


# filter on QC warning 
CM <- CM[CM$QC_Warning == "Pass", ]
CVD2 <- CVD2[CVD2$QC_Warning == "Pass", ]
CVD3 <- CVD3[CVD3$QC_Warning == "Pass", ]

CM <- CM[CM$SampleID %in% patientsPanel, ]
CVD2 <- CVD2[CVD2$SampleID %in% patientsPanel, ]
CVD3 <- CVD3[CVD3$SampleID %in% patientsPanel, ]

# proteins in panels
proteinsCVD2 <- unique(CVD2$Assay)
proteinsCVD3 <- unique(CVD3$Assay)
proteinsCM <- unique(CM$Assay)
proteinsPanel <- intersect(unique(dataExplore$Assay), unique(c(proteinsCVD2, proteinsCVD3, proteinsCM))) 

### patients and proteins intersecting 
proteins <- intersect(proteinsPanel, proteinsExplore)
patientsOverlap <- intersect(paste0("ae", patientsExplore), patientsPanel)
patientsDifferent <- setdiff(patientsPanel, paste0("ae", patientsExplore))

## dataframe with overlapp in patients 
dataPanelpatientsOverlap <- as.data.frame(matrix(nrow = length(patientsOverlap), ncol = length(proteins)))
rownames(dataPanelpatientsOverlap) <- patientsOverlap
colnames(dataPanelpatientsOverlap) <- proteins

for (protein in proteins){
  if (protein %in% proteinsCVD2){
    subset <- CVD2[CVD2$Assay == protein, ]
    for (patient in patientsOverlap){
      dataPanelpatientsOverlap[patient,protein] <- subset[subset$SampleID == patient, 'NPX'][1]
    }
  } 
  else if (protein %in% proteinsCVD3) {
    subset <- CVD3[CVD3$Assay == protein, ]
    for (patient in patientsOverlap){
      dataPanelpatientsOverlap[patient,protein] <- subset[subset$SampleID == patient, 'NPX'][1]
    } 
  }
  else if (protein %in% proteinsCM){
    subset <- CM[CM$Assay == protein, ]
    for (patient in patientsOverlap){
      dataPanelpatientsOverlap[patient, protein] <- subset[subset$SampleID == patient, 'NPX'][1]
    }
  }
}


## dataframe without Different in patients 
dataPanelpatientsDifferent <- as.data.frame(matrix(nrow = length(patientsDifferent), ncol = length(proteins)))
rownames(dataPanelpatientsDifferent) <- patientsDifferent
colnames(dataPanelpatientsDifferent) <- proteins

for (protein in proteins){
  if (protein %in% proteinsCVD2){
    subset <- CVD2[CVD2$Assay == protein, ]
    for (patient in patientsDifferent){
      dataPanelpatientsDifferent[patient,protein] <- subset[subset$SampleID == patient, 'NPX'][1]
    }
  } 
  else if (protein %in% proteinsCVD3) {
    subset <- CVD3[CVD3$Assay == protein, ]
    for (patient in patientsDifferent){
      dataPanelpatientsDifferent[patient,protein] <- subset[subset$SampleID == patient, 'NPX'][1]
    } 
  }
  else if (protein %in% proteinsCM){
    subset <- CM[CM$Assay == protein, ]
    for (patient in patientsDifferent){
      dataPanelpatientsDifferent[patient, protein] <- subset[subset$SampleID == patient, 'NPX'][1]
    }
  }
}



screen_missingness <- function(data, resultsfolder, name, row_thresh = 0.20, col_thresh = 0.05, plot = TRUE) {
  # order rows & columns by missingness 
  na_rows <- rowSums(is.na(data))
  na_cols <- colSums(is.na(data))
  
  ord_rows <- order(na_rows, decreasing = TRUE)
  ord_cols <- names(sort(na_cols, decreasing = FALSE))
  data_ordered  <- data[ord_rows, ord_cols]
  
  ## screen based on threshold (filter patients out first then select for proteins)
  good_rows <- names(na_rows[na_rows / ncol(data) <= row_thresh])
  na_cols_after <- colSums(is.na(data[good_rows, ]))
  good_cols <- names(na_cols_after[na_cols_after / length(good_rows)  <= col_thresh])
  
  v_break <-  length(good_cols) + 0.5
  h_break <- nrow(data_ordered) -  length(good_rows) + 0.5
  
  
 if (plot == TRUE) {
    df_melt <- melt(is.na(data_ordered))
    
    pdf(file = paste0(resultsFolder, name, ".pdf", sep = ""), width = 6, height = 1+nrow(data)*0.01)
    p <- ggplot(df_melt, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile(color = NA) +
      scale_fill_manual(values = c("TRUE"  = "gray60",
                                   "FALSE" = "grey90"),
                        name = "Missing?") +
      geom_vline(xintercept = v_break, colour =  "#4477AA", linewidth = 0.5) +
      geom_hline(yintercept = h_break, colour = "red", linewidth = 0.5) +
      labs(x = paste0("Proteins (n = ", ncol(data_ordered),")"), 
           y = paste0("Patients (n = ", nrow(data_ordered),")"), 
           title = ) +
     # labs(x = "Proteins", y = "Patients") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks  = element_blank(),
            panel.grid  = element_blank(),
            legend.position = "bottom") 

    print(p)
    dev.off()
    
  }
  
  list(
    clean_data   = data[good_rows, good_cols, drop = FALSE],
    kept_rows    = good_rows,
    kept_cols    = good_cols,
    row_thresh   = row_thresh,
    col_thresh   = col_thresh
  )
}

add_MACE <- function(data, MACEdefinition) {
  patient_ids <- sub("ae", "", rownames(data))
  col_event <- ifelse(MACEdefinition == "composite", "epcom.3years", "epmajor.3years")
  col_time <- ifelse(MACEdefinition == "composite", "ep_com_t_3years", "ep_major_t_3years")
  
  match_idx <- match(patient_ids, dataRiskFactors$STUDY_NUMBER)
  data$eventMACE <- dataRiskFactors[match_idx, col_event]
  data$timeMACE  <- dataRiskFactors[match_idx, col_time]
  
  return(data)
}

cox_regressions <- function(data){
  proteins <- setdiff(colnames(data),c("timeMACE", "eventMACE"))
  univariateFormulas <- sapply(proteins, function(x) as.formula(paste('Surv(timeMACE, eventMACE)~', x)))
  univariateModels <- lapply(univariateFormulas, function(x){coxph(x, data = data)})
  univariateResults <- lapply(univariateModels, function(x){ 
    x <- summary(x)
    p.value <- signif(x$coef[5], digits = 4)
    wald.test <- signif(x$wald["test"], digits = 4)
    beta <- signif(x$coef[1], digits = 3)                 # coefficient beta
    HR <- signif(x$coef[2], digits = 3)                   # exp(beta)
    HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
    HRextended <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    Assay <- data[data$Assay %in% rownames(x$coef), "Assay"][[1,1]]
    results <- c(Assay, beta, HR, HRextended, wald.test, p.value, "NA")
    names(results) <- c("beta", "HR", "HR (95% CI for HR)", "wald.test", "p.value", "adjusted.p.value")
    return(results)
  })
  univariateResults <- t(as.data.frame(univariateResults))
  univariateResults[ ,"adjusted.p.value"] <- signif(p.adjust(as.numeric(univariateResults[ ,"p.value"]), method = "BH"), 4)
  univariateResults <- univariateResults[order(univariateResults[ ,"p.value"]), ] 
  
  return(univariateResults)
}


add_explore_results <- function(df_results_panel, df_results_explore){
  for (assay in rownames(df_results_panel)){
    matched_row <- data.frame(df_results_explore)[df_results_explore[ ,"Assay"] == assay, ][1, ]
    df_results_panel[assay, "Explore_HR"] <- matched_row$HR
    df_results_panel[assay, "Explore_adjusted.p.value"] <- matched_row$adjusted.p.value
    df_results_panel[assay, "Explore_p.value"] <- matched_row$adjusted.p.value
  }

  df_results_panel$color <- "#BBBBBB" # Assuming default color is gray
  df_results_panel$color[df_results_panel$adjusted.p.value < 0.05 & df_results_panel$Explore_adjusted.p.value < 0.05] <- "red" 
  df_results_panel$color[df_results_panel$Explore_adjusted.p.value < 0.05 & df_results_panel$adjusted.p.value >= 0.05] <- "#AA3377"
  df_results_panel$color[df_results_panel$Explore_adjusted.p.value >= 0.05 & df_results_panel$adjusted.p.value < 0.05] <- "#4477AA"
  
  return(df_results_panel)
}

plot_HR <- function(data, resultsfolder, name, title){
  cor <- round(cor(as.numeric(data$Explore_HR), as.numeric(data$HR), method="spearman"), digits = 3)
  
  pdf(file = paste0(resultsFolder, name, ".pdf"), 
      width = 7, 
      height = 7.5)
  
  plot(as.numeric(data$Explore_HR), as.numeric(data$HR), 
     xlab = "Hazard Ratio OLINK Explore 3072", 
     ylab = "Hazard Ratio OLINK Target 96",  
     col = data$color,  
     pch = 20,  
     log = "xy", 
     xlim = c(0.3, 15), 
     ylim = c(0.3, 15), 
  )
  title(main = paste(title, cor), cex = 0.8, font.main = 3,)
  abline(h = 1, col = "darkgrey", lty = 2)
  abline(v = 1, col = "darkgray", lty = 2)
  
  # Create legend
  legend("topright", legend = c("Explore 3072: adj.p < 0.05,   Panel 96: adj.p < 0.05", 
                                   "Explore 3072: adj.p < 0.05,   Panel 96 adj.p >= 0.05", 
                                   "Explore 3072: adj.p >= 0.05, Panel 96 adj.p < 0.05", 
                                   "Explore 3072: adj.p >= 0.05, Panel 96 adj.p >= 0.05 "), 
         col = c("red", "#AA3377", "#4477AA", "#BBBBBB"), 
         pch = 20, pt.cex = 1.4, cex = 0.8, bty = "n")
  dev.off()
}


resultOverlap <- screen_missingness(dataPanelpatientsOverlap,
                                    resultsfolder, 
                                    "Missingness overlap cohort",
                                    row_thresh = 0.20,   # drop patients >20% missing
                                    col_thresh = 0.1)    # drop proteins >10% 

dataPanelOverlapClean <- mice(resultOverlap$clean_data, m=10, maxit=5, seed = 15, method="pmm")$data
dataPanelOverlapClean <- add_MACE(dataPanelOverlapClean, MACEdefinition)
resultsUnivariateCoxOverlap <- data.frame(cox_regressions(dataPanelOverlapClean))
resultsUnivariateCoxOverlap <- add_explore_results(resultsUnivariateCoxOverlap, univariateResults)
plot_HR(resultsUnivariateCoxOverlap, resultsFolder, "HRvalidationOverlap", "Spearmans correlation of hazard ratios overlapping datasets: ")

resultDifferent <- screen_missingness(dataPanelpatientsDifferent,
                                      resultsfolder, 
                                      "Missingness disjoint cohort",
                                      row_thresh = 0.20,   # drop patients >20% missing
                                      col_thresh = 0.1)    # drop proteins >10% missing
dataPanelDifferentClean <- mice(resultDifferent$clean_data, m=10, maxit=5, seed = 15, method="pmm")$data
dataPanelDifferentClean <- add_MACE(dataPanelDifferentClean, MACEdefinition)
resultsUnivariateCoxDifferent <- data.frame(cox_regressions(dataPanelDifferentClean))
resultsUnivariateCoxDifferent <- add_explore_results(resultsUnivariateCoxDifferent, univariateResults)
plot_HR(resultsUnivariateCoxDifferent, resultsFolder, "HRvalidationDisjoint", "Spearmans correlation of hazard ratios disjoint datasets: ")

