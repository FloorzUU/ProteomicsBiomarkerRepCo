

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
               "dataPatientXProteinRiskClusterMICE", "NotRemove", "dataUnprocessed", "dataRiskFactors")

### prep 
source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/ImportAndPrep.R"))
rm(list=ls()[! ls() %in% NotRemove])

### LOD
source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/LODhistogram.R"))
rm(list=ls()[! ls() %in% c(NotRemove, "dataUnproccesed")])

### Clusters 
source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/ClusterAssociation.R"))
rm(list=ls()[! ls() %in% NotRemove])

### Symptoms
source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/SymptomsAssociationAndNEFL.R"))
rm(list=ls()[! ls() %in% NotRemove])


### MACE analysis 
source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/Table1AndConfounders.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/CoxRegressionUniAndMulti.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/CorrelationPlotting.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/ForestPlots.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/InteractionMaceProteinAndSex.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/GSEApathwayAnalysis.R"))
rm(list=ls()[! ls() %in% NotRemove])

source(paste0(pathStart, "Pipeline1/GitHub/ScriptFinal/ROCbootstrapedKaplanMeierEstimator.R"))
rm(list=ls()[! ls() %in% NotRemove])




