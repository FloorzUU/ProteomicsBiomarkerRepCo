# Interaction Cox regression MACE ~ Protein * Sex

###

# Performing Cox regression with interaction 
# formula: MACE ~ Protein * Sex 
# The outcomes of the top protiens by smallest p.value are depicted (n = nProteinsTable = 30)
# The output is given in the form of a latex table 
#
# Required Input: 
# - dataInterestMACE
#  
# Resulting output: 
# - An printed code for a latex table presenting the top 30 outcomes of the interaction studies 

###

### import required libraries 
library(survival)         # version 3.3-1
library(kableExtra)       # version 1.3.4

### Number of proteins printed in latex table 
nProteinsTable <- 30 

### Cox regression for all proteins 
proteins <- colnames(dataInterestMACE)[grepl("^OID", colnames(dataInterestMACE))]
# setting up dataframe for results 
resultsInteraction <- as.data.frame(matrix(nrow = length(proteins), ncol = 6))
colnames(resultsInteraction) <- c("OlinkID", "Assay", "Hazard ratio", "p.value interaction", "adjusted.p.value", "LOD")

i = 1
for (protein in proteins){
  x <- coxph(as.formula(paste('Surv(timeMACE, eventMACE)~', protein, "*Sex", sep = "")), data = dataInterestMACE)
  x <- summary(x)
  p.value <- as.numeric(signif(x$coef[3, 5], digits = 4))
  HR <- as.numeric(signif(x$coef[3, 2], digits = 3))                   
  HR.confint.lower <- signif(x$conf.int[3, "lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[3, "upper .95"], 2)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  Assay <- data[data$OlinkID %in% protein, "Assay"][[1,1]]
  LOD <- data[data$OlinkID %in% protein, "MissingFreq"][[1]][1]
  resultsInteraction[i,] <- c(protein, Assay, HR, p.value, NA, LOD)
  i = i + 1
}
resultsInteraction$adjusted.p.value <- p.adjust(resultsInteraction$p.value, method = "BH")
  
### Create table with top results interaction 
# order on p.value multi
resultsInteraction <- resultsInteraction[order(as.numeric(resultsInteraction$p.value)), ]
resultsInteractionTop <- resultsInteraction[1:nProteinsTable,!colnames(resultsInteraction) %in% "OlinkID"]
rownames(resultsInteractionTop) <- resultsInteractionTop$Assay
resultsInteractionTop <- resultsInteractionTop[ ,!colnames(resultsInteractionTop) %in% "Assay"]
kable(resultsInteractionTop, "latex")

