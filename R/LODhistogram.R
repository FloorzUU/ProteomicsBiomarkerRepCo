# LOD histograms per panel 

###

# Plotting histograms depicting the fraction of LOD per proteins in a panel 
#
# Required Input: 
# - data
# - dataInterestMACE
#  
# Resulting output: 
# - An 8 panel plot with histograms depicting the number of proteins with a specific fraction of LOD

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
# None are required just base R 

### add MACE in each row of data 
dataSelected <- dataUprocessed
dataSelected$MACE <- NA
for (patient in rownames(dataInterestMACE)){
  dataSelected[dataSelected$SampleID == patient, "MACE"] <- dataInterestMACE[patient, "MACE"]
}
dataSelected$MACE <- as.factor(dataSelected$MACE)
dataSelected <- dataSelected[!is.na(dataSelected$MACE), ]

### LOD histograms per panel
LODpanel <- list()
proteins <- colnames(dataInterestMACE)[grepl("^OID", colnames(dataInterestMACE))]
# preperation list with missing value freq per panel
for (protein in proteins){
  # previously created a dataset called dataSelected with only data with MACE values 
  panel <- dataSelected[dataSelected$OlinkID == protein, "Panel"][[1]][1]
  value <- dataSelected[dataSelected$OlinkID == protein, "MissingFreq"][[1]][1]
  LODpanel[[panel]] <- c(LODpanel[[panel]], value)
} 

# plot
pdf(file = paste(resultsFolder, "LOD.pdf", sep = ""), width = 9, height = 6)
par(mfrow = c(2, 4))
for (panel in names(LODpanel)) {
  hist_data <- LODpanel[[panel]]
  breaks <- 100
  
  hist_result <- hist(hist_data, breaks = breaks, plot = FALSE)                 # Compute the histogram
  max_freq_index <- which.max(hist_result$counts)                               # Find the bin with the maximum frequency
  
  # Plot the histogram
  hist(hist_data, breaks = breaks, xlab = "Fraction NPX values below LOD",
       ylab = "Number of proteins", main = panel, ylim = c(0, 25),
       col = ifelse(seq_along(hist_result$counts) == max_freq_index, "darkgreen", "grey"),
       border = ifelse(seq_along(hist_result$counts) == max_freq_index, "darkgreen", "grey"))
  # Annotate the bar with the maximum frequency
  text(hist_result$mids[max_freq_index] + 0.095, 25 - 2,
       labels = paste("-", hist_result$counts[max_freq_index], sep = ""), pos = 3, col = "darkgreen", cex = 0.9)
}
dev.off()

rm(LODpanel, panel, value, hist_data, breaks, max_freq_index)
