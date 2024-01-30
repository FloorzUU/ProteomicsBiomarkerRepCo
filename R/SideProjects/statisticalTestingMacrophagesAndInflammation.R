### statistical testing Magcrophages and inflammation

### 

# Statistical testing with Wilcox of dichotomized data
# 
# Required input: 
# - long format protein data (in olinkLongFormatPath)
# - Clinical risk factors to be tested (in riskFactorPathSav)
# - FactorsOfInterest
#
# Resulting output: 
# - variable: outcomes (contains outcomes Wilxoc rank sum test all FactorsOfInterest)

###

# setting computer path 
if (Sys.info()[4][[1]] == "TURBO-PC-WESTER"){
  pathStart <- "C:/Users/Bas/Documents/FloorsFiles/MASTER/MajorInternship/"
} else if (Sys.info()[4][[1]] == "DESKTOP-2DADF21"){
  pathStart <- "C:/Users/Floor van der Zalm/Documents/MASTER/MajorInternship/"
} else {
  warning ("Please adjust file paths before executing")
}

# file locations
olinkLongFormatPath <- "RawData/NPXfiles/olink2022-012-080_Explore3072_EXTENDED_NPX_2023-07-10.csv"
riskFactorPathSav <- "RawData/ClinicalData/2022-21-03 AtheroExpress Database.sav"
factorsInterest <- c("hsCRP_plasma", "Macrophages.bin", "IPH.bin", "vessel_density_averaged")

### importing relevant libraries 
library(dplyr)
library(OlinkAnalyze)
library(haven)
library(sjmisc)

### importing data
# OLINK large dataset 
dataUnproccesed <- read_NPX(paste(pathStart, olinkLongFormatPath, sep = ""))

# data Risk factors 
dataRiskFactorsSav <- as.data.frame(read_sav(paste(pathStart, riskFactorPathSav, sep = "")))
rownames(dataRiskFactorsSav) <- dataRiskFactorsSav$STUDY_NUMBER
dataRiskFactorsSav <- dataRiskFactorsSav %>% select(factorsInterest)

### remove negative and positive controls, control assays and EXCLUDED assays
controlFind <- function(x) !grepl("control", x)
data <- dataUnproccesed[controlFind(dataUnproccesed$Assay), ]  # removes control assays
numbersOnly <- function(x) !grepl("\\D", x)
data <- data[numbersOnly(data$SampleID), ]          # removes negative and positive controls 
data <- data[data$QC_Warning != "EXCLUDED", ]       # removes EXCLUDED assays 
data <- data[data$MissingFreq < 0.25, ]             # removes assays with over 25% measurement below LOD

### patients 
patientsExclude <- c("32", "1746", "1921")
patients <- intersect(data$SampleID, rownames(dataRiskFactorsSav))
patients <- patients[!patients %in% patientsExclude]
for (patient in patients){
  data[data[ ,"SampleID"] == patient, colnames(dataRiskFactorsSav)] <- dataRiskFactorsSav[patient, colnames(dataRiskFactorsSav)]
}
print("Number of NA vlaues per column:")
colSums(is.na(data[ ,colnames(dataRiskFactorsSav)]))/length(unique(data$OlinkID))

outcomes <- list()
for (colTest in factorsInterest){
  dataTest <- data[!is.na(data[ ,colTest][1]), ]
  if (length(unique(dataTest[ ,colTest][[1]])) > 2){
    dataTest[ ,colTest][[1]] <- dicho(dataTest[ ,colTest][[1]])
  }
  dataTest[ ,colTest][[1]] <- as.factor(dataTest[ ,colTest][[1]])
  outcomes[[colTest]] <- olink_wilcox(dataTest, variable = colTest)
}

