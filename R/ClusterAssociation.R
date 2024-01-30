# Analysis association proteins and clusters

###

# Kruskal Wallis testing was performed with the olink package for the 
# transcriptome defined clusters. 
# Also a pie chart was made to visualize the number of patients per 
# cluster in the cohort 
#
# Required Input: 
# - dataPatientXProteinClusterMICE
# - data
#
# Resulting output:
# - pie chart per cluster
# - kruskal_resultsCluster (outcomes kruskal-wallis for all proteins) 
# - boxplots with a few top proteins 

###

### Required libraries 
library(OlinkAnalyze)     # version 3.5.1

### pie chart number of patients per cluster 
custom_colors <- c("#00C7E1", "#FE1F04", "#00559E", "#FFC700", "#006273")
pie(table(dataPatientXProteinRiskClusterMICE$cluster), radius = 1, 
    main = "Patients seperated by cluster", cex = 2, col = custom_colors)
rm(custom_colors)

### clusters
dataSelected <- data
dataSelected$Clusters <- NA
for (patient in patients){
  x <- dataPatientXProteinRiskClusterMICE[patient, "cluster"]
  dataSelected[dataSelected$SampleID  == sub("^ae", "", patient), "Clusters"] <- as.character(x)
}
dataSelected <- dataSelected[!is.na(dataSelected$Clusters), ]

kruskal_resultsCluster <- olink_one_non_parametric(df = dataSelected, variable = "Clusters")
x <- kruskal_resultsCluster
x <- x[order(x$Adjusted_pval), "OlinkID"]
x <- unique(x[ ,1])[1:4, ]
p1 <- olink_boxplot(df = dataSelected, 
                    variable = "Clusters", 
                    olinkid_list = x[[1]],
                    number_of_proteins_per_plot = 4)
p1[[1]]