### symptoms

###

# Associations between MACE and symptoms 
# Tested with Kruskal-Wallis testing with the olink package 
# NEFL springs out so we look if it has similarities with troponin 
#
# Required input: 
# - data
# - dataRiskFactors 
# - dataInterestMACE
# 
# Resulting output: 
# - Variable: kruskal_results for all proteins 
# - Variable: kruskal_results-post_hoc for all proteins 
# - boxplots of the most significant proteins in the R plots viewer 
# - NEFL plot - time OR in the R viewer

###

### Required libraries 
library(OlinkAnalyze)     # version 3.5.1

### Statistical testing symptoms and proteins 
patients <- unique(data$SampleID)
dataSelected <- data
dataSelected$Symptoms.4g <- NA
for (patient in patients){
  x <- dataRiskFactors[dataRiskFactors$STUDY_NUMBER == paste("ae", patient, sep = ""), "Symptoms.4g"]
  dataSelected[dataSelected$SampleID  == sub("^ae", "", patient), "Symptoms.4g"] <- as.character(x)
}
dataSelected <- dataSelected[!is.na(dataSelected$Symptoms.4g), ]
dataSelected[dataSelected$Symptoms.4g == "0", "Symptoms.4g"] <- "asymptomatic"
dataSelected[dataSelected$Symptoms.4g == "1", "Symptoms.4g"] <- "ocular"
dataSelected[dataSelected$Symptoms.4g == "2", "Symptoms.4g"] <- "TIA"
dataSelected[dataSelected$Symptoms.4g == "3", "Symptoms.4g"] <- "stroke"

kruskal_results <- olink_one_non_parametric(df = dataSelected, variable = "Symptoms.4g")
significantProteins  <- kruskal_results[kruskal_results$Threshold == "Significant", "OlinkID"]
kruskal_results_post_hoc <- olink_one_non_parametric_posthoc(df = dataSelected, 
                                                             olinkid_list = significantProteins$OlinkID, 
                                                             test = "kruskal", 
                                                             variable = "Symptoms.4g")
# Create a custom order for the "Symptoms.4g" variable
custom_order <- c("asymptomatic", "ocular", "TIA", "stroke")  # Replace with your desired order
# Reorder the levels of the "Symptoms.4g" variable
dataSelected$Symptoms.4g <- factor(dataSelected$Symptoms.4g, levels = custom_order)

x <- kruskal_results#_post_hoc[kruskal_results_post_hoc$contrast == "stroke - TIA", ]
x <- x[order(x$Adjusted_pval), "OlinkID"]
x <- unique(x[ ,1])[1:12, ]
p1 <- olink_boxplot(df = dataSelected, 
                    variable = "Symptoms.4g", 
                    olinkid_list = x[[1]],
                    number_of_proteins_per_plot = 4, 
                    posthoc_results = kruskal_results_post_hoc)
p1[[1]]
p1[[2]]


### NEFL closeup - plot OR time vs NPX value 
proteinInterest <- "OID20871"   # which protein to plot 
dataNEFL <- dataInterestMACE[!is.na(dataInterestMACE$Symptoms.4g),
                                            colnames(dataInterestMACE) %in% c(proteinInterest, "Symptoms.4g")]
dataOR <- dataRiskFactors[dataRiskFactors$STUDY_NUMBER %in% paste("ae", rownames(dataNEFL), sep = ""),
                          colnames(dataRiskFactors) %in% c("STUDY_NUMBER", "Time_event_OR")]

rownames(dataOR) <- sub("^ae", "", dataOR$STUDY_NUMBER)
dataOR <- dataOR[-1]
dataOR[-which(is.na(dataOR$Time_event_OR)), ]
dataNEFL[ ,"Time_event_OR"] <- dataOR[rownames(dataNEFL), "Time_event_OR"]
dataNEFL <- dataNEFL[dataNEFL$Time_event_OR >= 0, ]
dataNEFL <- dataNEFL[-which(is.na(dataNEFL$Time_event_OR)), ]
dataNEFL <- dataNEFL[!dataNEFL$Symptoms.4g == "0", ]
dataNEFL <- dataNEFL[rowSums(is.na(dataNEFL)) == 0, ]


ggplot(dataNEFL, aes(x = Time_event_OR, y = OID20871, color = factor(Symptoms.4g))) +
  geom_point(shape = 19, size = 1) +
  stat_smooth(, se = T) + 
  labs(x = "Time event OR (years)", y = "NPX value", title = "Time OR vs NPX value of NEFL") +
  theme_minimal() +
  scale_color_discrete(name = "Symptoms.4g") + 
  coord_cartesian(xlim = c(0, 150)) + 
  scale_color_discrete(name = "Symptoms.4g", labels = c("Ocular", "TIA", "Stroke"))  # Specify legend labels


