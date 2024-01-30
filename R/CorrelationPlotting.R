# Plotting correlations between proteins 

###

# Plotting two correlation plots 
# The first plot shows the correlation of the top 30 significant proteins
# significance was determined with the pvalues of the multivariate cox regression 
#
# The second plot depicts the correlations on all proteins with an 
# adjusted p.value < 0.05 
#
# Required Input: 
# - dataInterestMACE
# - pValues 
#  
# Resulting output: 
# - plot correlations of the top significant multivariatly proteins (n = 30)
# - plot correlations of all univariatly significant proteins (n = 319)

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
library(corrplot)         # version 0.92
library(ggplotify)        # version 0.1.2
library(caret)            # version 6.0-94
library(cowplot)          # version 1.1.1

### plotting correlation top proteins
numberProteins <- 30
cor_threshold <- 0.5
selectedProteins <- order(pValues$Multivariate)[1:numberProteins]
cor_matrix <- signif(cor(dataInterestMACE[ ,selectedProteins]), digits = 3)
length(findCorrelation(cor_matrix, cutoff = cor_threshold))

# set rownames to readable Assay names 
assayNames <- c()
for (i in 1:length(selectedProteins)) {
  assayNames <- c(assayNames,data[data$OlinkID == selectedProteins[i], "Assay"][[1]][1])
}
append_numbers <- function(x) {
  ave(x, x, FUN = function(y) if(length(y) > 1) paste0(y, "_", seq_along(y)) else y)
}
assayNames <- append_numbers(assayNames)
colnames(cor_matrix) <- assayNames
rownames(cor_matrix) <- assayNames
pdf(file = paste(resultsFolder, "TopProteinsCor.pdf", sep = ""), width = 4.5, height = 4.5)
par(mar = c(1,15,150,10))
p <- corrplot(cor_matrix, method = 'shade', pch.col = 'black', tl.col = "black", 
              order = "hclust", tl.cex = 0.6)
dev.off()

### plotting correlation all proteins significant univariate
selectedProteinsIdx <- as.vector(pValues$Univariate < 0.05)
dataSelected <- dataInterestMACE[1:dim(pValues)[1]]
cor_matrix <- cor(signif(dataSelected[ ,selectedProteinsIdx], digits = 3))
selectedProteins <- rownames(pValues[pValues$Univariate < 0.05, ])
colnames(cor_matrix) <- selectedProteins
rownames(cor_matrix) <- selectedProteins
selectedFeaturesCor <- length(findCorrelation(cor_matrix, cutoff = cor_threshold))
pdf(file = paste(resultsFolder, "significantProteinsCor.pdf", sep = ""), width = 9, height = 8)
p <- corrplot(cor_matrix, type = "full", order = "hclust", method = "color", tl.pos = 'n')
dev.off()


### corrplot double measured assays 
doubles <- names(which(table(data$Assay) > length(unique(data$SampleID))))
breaks <- seq(0.6, 1, length.out = 101)
plot_list <- list()
i <- 1
for (double in doubles){
  cor_matrix <- cor(dataInterestMACE[ ,unique(data[data$Assay == double, "OlinkID"])[[1]]])
  p <- as.ggplot(pheatmap(cor_matrix, breaks = breaks, type = "full", order = "hclust", method = "color",
                          tl.pos = 'n',cluster_rows = F,cluster_cols = F, main = double, legend = FALSE))
  plot_list[[i]] <- p
  i <- i + 1
} 
# Arrange plots in a 4 by 3 grid using cowplot
final_plot <- plot_grid(plotlist = plot_list)
# Display or save the final plot
ggsave(file = paste(resultsFolder, "corDoubles.pdf", sep = ""), 
       final_plot, width = 6, height = 4)

rm(i, doubles, breaks, p, final_plot, plot_list, cor_matrix)

