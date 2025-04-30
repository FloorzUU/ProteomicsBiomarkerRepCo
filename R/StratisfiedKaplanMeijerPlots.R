### Construct Kaplan-Meijer curves (stratisfied by terciles)

# In this script we construct by NPX stratified kaplan-Meijer curves 
# using survminer. We stratified the patients in three groups based on NPX value 
# the low, mid and high group. 
#
# required input
# - dataInterestMACE 
#
# Edit selected proteins by changing variables:
# selectionAssay and the corresponding: selectionOlinkID
# 
# resulting output
# - a figure with kaplan-meijer plots in pdf format

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

### import libraries 
library(survminer)    # version 0.4.9
library(ggplot2)      # version 3.5.1

### Kaplan-Meier plot
selectionAssay <- c("CHCHD10", "IFI30", "AREG", "FABP5", "WFDC2", "GDF15", "NT-proBNP", "HAVCR1", "CHCHD6" )
selectionOlinkID <- c("OID31116","OID31505","OID21357","OID21043","OID21505","OID20251","OID20125","OID21422", "OID31426")

plots_list <- list()
dataInterestMACE$timeMACE <- as.numeric(dataInterestMACE$timeMACE)
# Loop through each protein in selectionOlinkID
for (i in seq_along(selectionOlinkID)){
  protein <- selectionOlinkID[i]
  
  quantiles <- c(0.33, 0.66)
  quantile_values <- quantile(dataInterestMACE[ ,protein], probs = quantiles)
  dataQuant <- as.numeric(findInterval(dataInterestMACE[ ,protein], quantile_values))
  dataQuant <- dataQuant + 1
  dataQuant <- as.factor(dataQuant)
  dataKaplanMeier <- dataInterestMACE[ ,c("timeMACE", "eventMACE")]
  dataKaplanMeier[ ,paste("Q_", selectionOlinkID[i], sep = "")] <- dataQuant
  dataKaplanMeier[ ,selectionOlinkID[i]] <- dataInterestMACE[ ,protein]
  
  
  # Kaplan-Meier curve
  fit <- survfit(as.formula(paste('Surv(timeMACE, eventMACE) ~', paste("Q_", selectionOlinkID[i], sep = ""))), data = dataKaplanMeier)
  p <- ggsurvplot(fit,
                  conf.int = FALSE,
                  ggtheme = theme_bw(),
                  palette = c("#924BBD", "#2E9FDF", "#60A890"),
                  legend.title = selectionAssay[i],
                  legend.labs = c("Tercile 1", "Tercile 2", "Tercile 3"),
                  xlab = "Time (Years)",
                  pval = TRUE, pval.coord = c(0.05, 0.05),
                  main = "Kaplan-Meier survival curve",
                  pval.size = 4, legend = c(0.21, 0.32))
  
  # Add the Kaplan-Meier plot to the list
  plots_list[[i]] <- p$plot
}

# Arrange plots in a 4 by 3 grid using cowplot
final_plot <- plot_grid(plotlist = plots_list, ncol = 4)
plot(final_plot)

# Display or save the final plot
ggsave(file = paste(resultsFolder, "KaplanMeierPlot_3.pdf", sep = ""), 
       final_plot, width = 12, height = 9.6)



