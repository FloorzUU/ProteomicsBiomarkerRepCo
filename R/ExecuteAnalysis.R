

if (Sys.info()[4][[1]] == "TURBO-PC-WESTER"){
  pathStart <- "C:/Users/Bas/Documents/FloorsFiles/MASTER/MajorInternship/"
} else if (Sys.info()[4][[1]] == "DESKTOP-2DADF21"){
  pathStart <- "C:/Users/Floor van der Zalm/Documents/MASTER/MajorInternship/"
} else {
  warning ("Please adjust file paths before executing")
}

NotRemove <- c("dataInterestMACE", "includedPatientData", "data", "MACEdefinition", 
               "pValues", "univariateResults", "pathStart", "gsea_results", 
               "kruskal_results", "kruskal_results_post_hoc", "kruskal_resultsCluster", 
               "dataPatientXProteinRiskClusterMICE", "NotRemove", "dataUnproccesed", "dataRiskFactors")

### prep 
source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/ImportAndPrep.R"))
rm(list=ls()[! ls() %in% NotRemove])

### LOD
source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/LODhistogram.R"))
rm(list=ls()[! ls() %in% c(NotRemove)])

### Clusters 
source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/ClusterAssociation.R"))
rm(list=ls()[! ls() %in% NotRemove])

### Symptoms
source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/SymptomsAssociationAndNEFL.R"))
rm(list=ls()[! ls() %in% NotRemove])

### MACE analysis 
source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/Table1AndConfounders.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/CoxRegressionUniAndMulti.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/CorrelationPlotting.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/ForestPlots.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/InteractionMaceProteinAndSex.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/GSEApathwayAnalysis.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/ROCbootstrapedKaplanMeierEstimator.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "StudentAssistentWork/ResearchCompendium/ProteomicsBiomarkerRepCo/R/StratisfiedKaplanMeijerPlots.R"))
rm(list=ls()[! ls() %in% NotRemove])


