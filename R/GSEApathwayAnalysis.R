# Pathway analysis 

###

# Performs Gene set enrichtment analysis with GO biological processes database
# First data is preped (code nearly the same as in Olink Vingette)
# then the analysis is performed, the ranked gene list is made with the beta  
# values from univariate cox regressions
#
# Required input:
# - data
# - UnivariateResults 
#
# Resulting output: 
# - variable: GSEA results
# - dotplot all significant pathways 
# - cnet plot 
# - treeplot all significant pathways, grouped in 5 groups 

###


### setting computer path 
if (Sys.info()[4][[1]] == "TURBO-PC-WESTER"){
  pathStart <- "C:/Users/Bas/Documents/FloorsFiles/MASTER/MajorInternship/"
} else if (Sys.info()[4][[1]] == "DESKTOP-2DADF21"){
  pathStart <- "C:/Users/Floor van der Zalm/Documents/MASTER/MajorInternship/"
} else {
  warning ("Please adjust file paths before executing")
}
resultsFolder <- paste(pathStart, "Pipeline1/Results/pipelineLargeDataset/predictionPreperation/MACE/", sep = "")


### import required libraries 
library(clusterProfiler)    # version 4.4.4
library(msigdbr)            # version 7.5.1
library(dplyr)              # version 1.1.3
library(DOSE)               # version 3.22.1
library(enrichplot)         # version 1.16.2
library(scales)             # version 1.2.1
library(ggplot2)            # version 3.4.2

# plot cnet? 
cnet <- F


GSEA_analysis <- function(data, test_results, subcategory = NULL){
  # checks for na's and picks the best duplicates
  data <- prep_input(data)
  # removes removed duplicates in test results
  test_results <- test_results %>% filter(OlinkID %in% unique(data$OlinkID))
  
  # create Genelist with Assay names
  estimate <- as.numeric(test_results$beta) # beta is the estimate given by the cox regression
  names(estimate) <- test_results$Assay
  geneList <- sort(estimate, decreasing = TRUE)
  
  # retrieve gene sets from MsigDBr - GO specifically for biological processes
  # possible subcategories may be: GO:BP, GO:CC, GO:MF, HPO
  df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = subcategory)
  # assays not represented in the selected msidgbdr database
  message(paste0(length(setdiff(names(geneList), df$gene_symbol)),
                 " assays are not found in the database. Please check the Assay names for the following assays:\n ",
                 toString(setdiff(names(geneList), df$gene_symbol))))
  df <- as.data.frame(df[ ,colnames(df) %in% c("gs_name", "gene_symbol")])
  
  GSEA <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = df, pvalueCutoff = 0.05)
  return(GSEA)
}

# Data and test results filtered for highest detectibility in duplicate assay names
# A function written in OLINK vingette and slightly adjusted
prep_input <- function(data) {
  # Remove NA data and filter for valid Olink IDs
  npx_check <- data
  data <- data %>% filter(stringr::str_detect(OlinkID, "OID[0-9]{5}"))
  # Filter highest detectibility for repeated IDs
  olink_ids <- data %>%
    dplyr::filter(!stringr::str_detect(string = SampleID, pattern = "CONTROL*.")) %>%
    dplyr::filter(toupper(QC_Warning) == "PASS") %>%
    dplyr::mutate(Detected = as.numeric(NPX > LOD)) %>%
    dplyr::select(OlinkID, Assay, Detected) %>%
    dplyr::group_by(OlinkID, Assay) %>%
    dplyr::summarise(Sum = (sum(Detected) / n())) %>%
    dplyr::arrange(dplyr::desc(Sum)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Assay, .keep_all = TRUE) %>%
    dplyr::pull(OlinkID)
  
  data <- data %>% dplyr::filter(OlinkID %in% olink_ids)
  
  return(data)
}

### prep univariate results
# So it can be used to make the gene list 
univariateResults <- as.data.frame(univariateResults)
univariateResults[ ,"OlinkID"] <- rownames(univariateResults)

### data prep and pathway analysis 
gsea_results <- GSEA_analysis(data = data, test_results = univariateResults,
                              subcategory = "GO:BP")

### dotplot - activated and suppresed seperated
pdf(file = paste(resultsFolder, "gseaPlot.pdf", sep = ""),
    width = 12, height = 9)
dotplot(gsea_results, showCategory = dim(gsea_results@result)[1], title = "Enriched Pathways" , split=".sign") +
  facet_grid(.~.sign) + 
  scale_y_discrete(labels = label_wrap(100)) + 
  theme(axis.text.y=element_text(size=9)) + 
  theme(axis.text.x=element_text(size=9))
dev.off()  

### cnet plot
if (cnet == T){
  gsea_resultsPairwise <- pairwise_termsim(gsea_results)
  emapplot(gsea_resultsPairwise, showCategory = 15, cex_label_category = 0.5)
  ridgeplot(gsea_results, showCategory = 15) + 
    labs(x = "enrichment distribution") +
    theme(axis.text.y=element_text(size=9)) 
  cnetplot(gsea_resultsPairwise, node_label="gene", 
           cex_label_gene = 0.5, showCategory = 15,
           colorEdge = FALSE)#, foldChange = geneList)
  
  p <- cnetplot(gsea_resultsPairwise, node_label="gene", cex_label_gene = 0.9, showCategory = 15) 
  pdf(paste(resultsFolder, "cnetPlot.pdf", sep = ""), width=12, height=9)
  print(p)
  graphics.off()
}


### Treeplot 
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

pdf(file = paste(resultsFolder, "gseaTreePlot.pdf", sep = ""),
    width = 12, height = 5)
p <- treeplot(gsea_resultsPairwise, cex_category = 0.7, 
              group_color = c("#52B6A8", "#8ADBE6", "#6699FF", "#4472C4", "#9B8EDA"),
              fontsize = 2.5, label_format_tiplab = 200, cex_label_gene = 0.5, nWords = 1, 
              label_format = 0.6, legend = FALSE, showCategory = dim(gsea_results@result)[1]) 
p$layers[[7]]$aes_params$size <- 2.7
plot(p)
dev.off()