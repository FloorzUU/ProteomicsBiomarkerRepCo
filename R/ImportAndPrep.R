# pipeline MACE and cluster import and prep data 

### 

# Reformat OLINK dataset and merge with MACE & cluster data 
# format will become a dataframe with patient X protein 
# Uniprot codes will be used to identify proteins 
#
# An Additional dataframe created with riskfactors for all patients included 
# The number of missing datapoints per riskfactor is checked 
# The threshold for inclusion was set at 0.1
#
# Missing data for included risk factors are imputed using MICE 
# A dataframe without the imputed data is saved aswell
#
# In this script a table 1 for MACE is also created 
# And a setup of the COX regression is made though the exploration of relevant 
# column reporting on MACE and the time at which MACE occurs 
# For now it is mainly focused on the within 3 year columns, but this still requires 
# some more attention. As timing has some inconsistent columns
#
# Required input: 
# - long format OLINK data
# - risk factor dataset to determine MACE and the other risk factors
# - dataset cluster per patient
# - input list of risk factors
# 
# Resulting output: 
# - reformatted protein data with MACE column 
# - reformatted protein data with MACE column and risk factors  
# - reformatted protein data with MACE column and risk factors without na values 
# - reformatted protein and risk factors MACE with only patients with time and event data present: dataInterestMACE
# - reformatted protein data with cluster column 
# - reformatted protein data with cluster column and risk factors  
# - reformatted protein data with cluster column and risk factors without na values 
# - report on number missing values
# - report distribution risk factor and risk factor properties 
# All these values will be stored in the resultsfolder: predictionPreperation 

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
olinkLongFormatPath <- "RawData/NPXfiles/olink2022-012-080_Explore3072_EXTENDED_NPX_2023-07-10.csv"
riskFactorPath <- "RawData/ClinicalData/Atheroexpress_Database_New_2022_MB_2.xlsx"
riskFactorPathSav <- "RawData/ClinicalData/2022-21-03 AtheroExpress Database.sav"
clusterPath <- "RawData/Clusters/batch2_batch1_combined_v4_5clusters.txt"
resultsFolder <- paste(pathStart, "Pipeline1/Results/pipelineLargeDataset/predictionPreperation/", sep = "")

### risk factors to include, some may still be excluded due to missing values
includedPatientData <- c("Gender", "SmokerCurrent", "DM.composite", "CAD_history", "Stroke_history", "PAOD", "Symptoms.4g",
                         "GFR_MDRD", "Hypertension1", "systolic", "diastoli", "Age", "BMI", "hb", "ht", "creat", 
                         "homocys", "totalchol", "TG_clinic", "LDL", "HDL", "glucose", "hsCRP_plasma", "risk614", "Med.statin.derived")
percentagaNaExclude <- as.numeric(0.1)

### patients excluded 
patientsExclude <- c("32", "1746", "1921")

### import required libraries 
library(OlinkAnalyze)     # version 3.4.1
library(svDialogs)        # version 1.1.0
library(excel.link)       # version 0.9.11
library(haven)            # version 2.5.3
library(mice)             # version 3.15.0
library(dplyr)            # version 1.1.3

### importing data
# OLINK large dataset 
dataUnproccesed <- read_NPX(paste(pathStart, olinkLongFormatPath, sep = ""))

# RiskFactors 
password <- dlgInput("Enter password protected excel files: ", Sys.info()["user"])$res
dataRiskFactors <- xl.read.file(paste(pathStart, riskFactorPath, sep = ""), password = password, write.res.password = password)

# try different database
dataRiskFactorsSav <- as.data.frame(read_sav(paste(pathStart, riskFactorPathSav, sep = "")))

# clusters assigned to patients 
dataClusters <- read.table(paste(pathStart, clusterPath, sep = ""), header = TRUE)
rownames(dataClusters) <- sub("^ae", "", dataClusters$STUDY_NUMBER)
dataClusters <- dataClusters[2]
colnames(dataClusters) <- "cluster"

rm(password, olinkLongFormatPath, riskFactorPath)

### merge dataRiskFactors and dataRiskFactorsSav
# dataRiskFactorSav is better updated and is used as basis 
# it does miss a few columns which are present dataRiskFactors
columnsMissing <- dataRiskFactors[ ,!colnames(dataRiskFactors) %in% intersect(colnames(dataRiskFactors), colnames(dataRiskFactorsSav))]
columnsMissingInterest <- intersect(colnames(columnsMissing), includedPatientData)
# add missing columns 
dataRiskFactorsSav$STUDY_NUMBER <- paste("ae", dataRiskFactorsSav$STUDY_NUMBER, sep = "")
dataRiskFactors <- left_join(dataRiskFactorsSav, 
                             dataRiskFactors[ ,colnames(dataRiskFactors) %in% c("STUDY_NUMBER", columnsMissingInterest)], 
                             by = "STUDY_NUMBER")

rm(columnsMissing, columnsMissingInterest, dataRiskFactorsSav)


### remove negative and positive controls, control assays and EXCLUDED assays
controlFind <- function(x) !grepl("control", x)
data <- dataUnproccesed[controlFind(dataUnproccesed$Assay), ]  # removes control assays
numbersOnly <- function(x) !grepl("\\D", x)
data <- data[numbersOnly(data$SampleID), ]          # removes negative and positive controls 
data <- data[data$QC_Warning != "EXCLUDED", ]       # removes EXCLUDED assays 
data <- data[data$MissingFreq < 0.25, ]

### long format to patient x protein matrix 
patients <- sort(unique(data$SampleID))
proteins <- names(which(table(data$OlinkID) == length(patients)))

dataPatientXProtein <- matrix(nrow = length(patients), ncol = length(proteins))
rownames(dataPatientXProtein) <- patients 
colnames(dataPatientXProtein) <- proteins
# reformat data 
for (col in colnames(dataPatientXProtein)){
  df <- data[data$OlinkID == col, "NPX"][[1]]
  names(df) <- data[data$OlinkID == col, "SampleID"][[1]]
  df <- df[sort(names(df))]                                     # sorting all values in order of patients in dataframe
  dataPatientXProtein[ ,col] <- df                              # fill in matrix
}

dataPatientXProtein <- as.data.frame(dataPatientXProtein)
print(paste("Number Na's in dataPatientXProtein",  sum(is.na(dataPatientXProtein))))
# remove patients with questions about consent
dataPatientXProtein <- dataPatientXProtein[!rownames(dataPatientXProtein)
                                           %in% patientsExclude, ]   
patients <- patients[!patients %in% patientsExclude]            # remove excluded patients to prevent issues later

rm(df, controlFind, numbersOnly)


### add MACE data to dataframe and select for those who have it 
MACEdefinition <- "major"
if (MACEdefinition == "composite"){
  MACE <- dataRiskFactors$EP_composite
} else if (MACEdefinition == "major"){
  MACE <- dataRiskFactors$EP_major
}
names(MACE) <- sub("^ae", "", dataRiskFactors$STUDY_NUMBER)
dataPatientXProteinMACE <- dataPatientXProtein
dataPatientXProteinMACE$MACE <- MACE[patients]  
dataPatientXProteinMACE <- na.omit(dataPatientXProteinMACE)
save(dataPatientXProteinMACE, file = paste(resultsFolder, "MACE/dataProteinMACE.Rdata", sep = "")) # save data


### add cluster data to dataframe and select for those who have it 
patientsCluster <- intersect(rownames(dataPatientXProtein), rownames(dataClusters))
dataPatientXProteinCluster <- cbind(dataPatientXProtein[patientsCluster, ], dataClusters[patientsCluster, "cluster"])
colnames(dataPatientXProteinCluster)[dim(dataPatientXProteinCluster)[2]] <- "cluster"
dataPatientXProteinCluster <- na.omit(dataPatientXProteinCluster)
save(dataPatientXProteinCluster, file = paste(resultsFolder, "Cluster/dataProteinCluster.Rdata", sep = "")) # save data


### add risk factors 
mergeDataAndRiskfactors <- function(data, dataRiskFactors, includedPatientData, percentagaNaExclude, resultsPath, y){
  # risk factors 
  rownames(dataRiskFactors) <- sub("^ae", "", dataRiskFactors$STUDY_NUMBER)
  dataRiskFactorsSelected <- dataRiskFactors[rownames(data), colnames(dataRiskFactors) %in% includedPatientData]
  
  # Change term Gender to Sex as the term sex better represents the meaning of the data 
  colnames(dataRiskFactorsSelected)[colnames(dataRiskFactorsSelected) == "Gender"] <- "Sex"
  
  # include those risk factors filled for more than 10% of values in dataset of patients with cluster assigned 
  # from the list of included risk factors
  include <- names(which(colSums(is.na(dataRiskFactorsSelected)) <= percentagaNaExclude*dim(dataRiskFactorsSelected)[1]))
  dataRiskFactorsIncluded <- dataRiskFactorsSelected[ ,include]
  
  # report on risk factor that couldn't be included 
  exclude <- includedPatientData[!includedPatientData %in% include]
  write("The following risk factors were excluded:", file = paste(resultsPath, y, "/excludedRiskFactor.txt", sep = ""), append = FALSE)
  write(exclude, file = paste(resultsPath, y, "/excludedRiskFactor.txt", sep = ""), append = TRUE)
  # number of na values in included risk factors
  naCountAll <- dataRiskFactorsIncluded %>% is.na() %>% colSums() 
  write("Na counts for all patients in dataset used as input \nNumber na values per risk factor", file = paste(resultsPath, y, "/excludedRiskFactor.txt", sep = ""), append = TRUE)
  write.table(naCountAll, file = paste(resultsPath, y, "/excludedRiskFactor.txt", sep = ""), append = TRUE, col.names = F)
  
  dataAllNoMICE <- cbind(data, dataRiskFactorsIncluded)
  save(dataAllNoMICE, file = paste(resultsPath, y, "/dataProteinRiskNoMice", y, ".Rdata", sep = ""))
  
  return(dataAllNoMICE)
}


### execute add risk factors 
# MACE 
dataPatientXProteinRiskMACE <- mergeDataAndRiskfactors(dataPatientXProteinMACE, dataRiskFactors, includedPatientData, percentagaNaExclude, resultsFolder, "MACE")

# cluster
dataPatientXProteinRiskCluster <- mergeDataAndRiskfactors(dataPatientXProteinCluster, dataRiskFactors, includedPatientData, percentagaNaExclude, resultsFolder, "Cluster")


### impute missing data 
MultipleImputation <- function(data, resultsPath, y){
  set.seed(111)
  print("It is assumed the data to be Missing At Random (MAR), as mice shouldn't be used on MNAR data. Keep this in mind")
  print("For more info if I want to improve my knowledge: https://stefvanbuuren.name/fimd/sec-idconcepts.html")
  
  # making factors with more than 6 unique values numeric 
  colToNum <- names(data)[sapply(data, function(x) length(unique(x))) > 5]
  data[colToNum] <- lapply(data[colToNum], as.numeric)
  colToFact <- names(data)[sapply(data, function(x) length(unique(x))) <= 5]
  data[colToFact] <- lapply(data[colToFact], as.factor)
  if(all(sapply(data, function(x) is.numeric(x) || is.character(x)))) {
    print("check the classes of your input data") 
  } else {
    print("input classes all good")
  }
  
  columns_with_na <- names(which(colSums(is.na(data)) > 0))
  # execute imputation of MICE with default methods settings
  imputateModel <- mice(data = data[ ,columns_with_na], m = 20, maxit = 100, seed = 13, method = "pmm")
  # extract a imputed and complete data set 
  dataComplete <- complete(imputateModel)
  dataAll <- data
  dataAll[ ,columns_with_na] <- dataComplete
  
  # save data 
  save(dataAll, file = paste(resultsPath, y, "/dataProteinRiskMICE", y,".Rdata", sep = ""))
  
  return(dataAll)
}

### execute imputation missing data with MICE
# MACE
dataPatientXProteinRiskMACEMICE <- MultipleImputation(dataPatientXProteinRiskMACE, resultsFolder, "MACE")

# Cluster 
dataPatientXProteinRiskClusterMICE <- MultipleImputation(dataPatientXProteinRiskCluster, resultsFolder, "Cluster")


### MACE data with time 
MACE <- data.frame()
if (MACEdefinition == "composite"){
  MACE <- dataRiskFactors[sub("^ae", "", dataRiskFactors$STUDY_NUMBER) %in% patients, c("ep_com_t_3years", "epcom.3years", "STUDY_NUMBER")]
} else if (MACEdefinition == "major"){
  MACE <- dataRiskFactors[sub("^ae", "", dataRiskFactors$STUDY_NUMBER) %in% patients, c("ep_major_t_3years","epmajor.3years", "STUDY_NUMBER")]
}
rownames(MACE) <- MACE$STUDY_NUMBER
MACE <- MACE[ ,!colnames(MACE) %in% c("STUDY_NUMBER")]
MACE <- MACE %>% drop_na()
MACE <- mutate_all(MACE, function(x) as.numeric(as.character(x)))

# select data with MACE columns and prepare for COX regression 
dataInterestMACE <- dataPatientXProteinRiskMACEMICE[sub("^ae", "", rownames(MACE)), ]
dataInterestMACE[sub("^ae", "", rownames(MACE)), "timeMACE"] <- MACE[ ,1]
dataInterestMACE[sub("^ae", "", rownames(MACE)), "eventMACE"] <- MACE[ ,2]

