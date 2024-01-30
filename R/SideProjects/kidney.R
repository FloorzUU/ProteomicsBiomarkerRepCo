library(rentrez)
library(httr)
library(jsonlite)
library(survival)

weightProteins <- c()

get_molecular_weight <- function(uniprot_code) {
  base_url <- "https://www.ebi.ac.uk/proteins/api/proteins"
  full_url <- paste0(base_url, "/", uniprot_code)
  
  response <- GET(full_url, accept("application/json"))
  
  protein_data <- content(response, "text", encoding = "UTF-8")
  protein_data <- fromJSON(protein_data)
  molecular_weight <- protein_data[["sequence"]][["mass"]]
  return(molecular_weight)
}

resultsKindey <- data.frame()
for (protein in colnames(dataPatientXProtein)){
  UniProt <- data[data$OlinkID == protein, "UniProt"][[1]][1]
  if (!UniProt %in% c("NTproBNP", "Q29980_Q29983", "P21217_Q11128", "P29459_P29460", 
                      "P0DUB6_P0DTE7_P0DTE8", "Q8IXS6")){
    resultsKindey[protein, "weight"] <- get_molecular_weight(UniProt)
  }
}

confounders <- c("",  "+ GFR_MDRD")
confoundersName <- c("univariate", "eGFR")
i <- 1
for (confounder in confounders){
  multivariableFormulas <- sapply(rownames(resultsKindey),
                                  function(x){as.formula(paste('Surv(timeMACE, eventMACE)~', x, confounder))})
  multivariableModels <- lapply(multivariableFormulas, function(x){coxph(x, data = dataInterest)})
  # Extract data 
  multivariableResults <- lapply(multivariableModels, function(x){ 
    x <- summary(x)
    p.value <- signif(x$coef[1,5], digits = 4)
    wald.test <- signif(x$wald["test"], digits = 4)
    beta <- signif(x$coef[1, 1], digits = 3)                  # coefficient beta
    HR <- signif(x$coef[1, 2], digits = 3)                    # exp(beta)
    HR.confint.lower <- signif(x$conf.int[1, "lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[1, "upper .95"], 2)
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    Assay <- data[data$OlinkID %in% rownames(x$coef), "Assay"][[1,1]]
    HRshort <- signif(x$coef[1, 2], digits = 3)
    results <- c(Assay, beta, HR, HRshort, wald.test, p.value, "NA")
    names(results) <- c("Assay", "beta", "HR (95% CI for HR)", "HR", "wald.test", "p.value", "adjusted.p.value")
    return(results)
  })
  multivariableResults <- t(as.data.frame(multivariableResults))
  multivariableResults[ ,"adjusted.p.value"] <- signif(p.adjust(as.numeric(multivariableResults[ ,"p.value"]), method = "BH"), 4)
  
  resultsKindey[ ,paste0(confoundersName[i], ".p.value")] <- as.numeric(multivariableResults[ ,"p.value"])
  resultsKindey[ ,paste0(confoundersName[i], ".HR")] <- as.numeric(multivariableResults[ ,"HR"])
  resultsKindey[ ,paste0(confoundersName[i], ".beta")] <- as.numeric(multivariableResults[ ,"beta"])
  
  i <- i + 1
} 

resultsKindey <- resultsKindey[order(resultsKindey$weight), ]
plot(resultsKindey$weight, resultsKindey[ ,"univariate.HR"] - resultsKindey[ ,"eGFR.HR"])
plot(resultsKindey$weight, resultsKindey[ ,"univariate.beta"] - resultsKindey[ ,"eGFR.beta"])

plot(resultsKindey[resultsKindey$univariate.p.value < 0.05, "weight"], resultsKindey[resultsKindey$univariate.p.value < 0.05 ,"eGFR.beta"])
points(resultsKindey[resultsKindey$univariate.p.value < 0.05, "weight"], resultsKindey[resultsKindey$univariate.p.value < 0.05 ,"univariate.beta"], col = "blue")
# Define a color palette from red to blue
color_palette <- colorRampPalette(c("red", "white", "blue"))
num_colors <- 100
colors <- color_palette(num_colors)
values <- resultsKindey[resultsKindey$univariate.p.value < 0.05 ,"univariate.beta"] - resultsKindey[resultsKindey$univariate.p.value < 0.05 ,"eGFR.beta"]
scaled_values <- scale(values, center = min(values), scale = max(values) - min(values))
color_indices <- cut(scaled_values, breaks = num_colors, labels = FALSE)

i <- 1
for (protein in rownames(resultsKindey[resultsKindey$univariate.p.value < 0.05, ])){
 segments(resultsKindey[protein, "weight"], resultsKindey[protein, "univariate.beta"], 
          resultsKindey[protein, "weight"], resultsKindey[protein, "eGFR.beta"], 
          col = colors[color_indices[i]], lwd = 2)
  i <- i + 1
}


effect <- rownames(resultsKindey[resultsKindey$univariate.p.value < 0.05 & resultsKindey$eGFR.p.value > 0.05, ])
NOeffect <- rownames(resultsKindey[resultsKindey$univariate.p.value < 0.05 & resultsKindey$eGFR.p.value < 0.05, ])
summary(resultsKindey[effect,"weight"])
summary(resultsKindey[NOeffect,"weight"])

library(vioplot)
x1 <- resultsKindey$weight[resultsKindey$univariate.p.value < 0.05 & resultsKindey$eGFR.p.value > 0.05]
x2 <- resultsKindey$weight[resultsKindey$univariate.p.value < 0.05 & resultsKindey$eGFR.p.value < 0.05]
vioplot(x1, x2, names=c("confounder effect", "confounder no effect"),
        col="gold")
boxplot(x1, x2)



