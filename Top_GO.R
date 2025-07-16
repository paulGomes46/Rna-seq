# Load necessary libraries
library(topGO)
library(readxl)
library(writexl)
library(Rgraphviz)

# Load background gene-to-GO mapping
geneID2GO <- readMappings(file = "C:/Users/gomes/OneDrive/Documents/PhD/GO_enrichment/GO_background.txt")
str(head(geneID2GO))

# Load all genes with binary presence/absence (1/0) from the Excel file
all_genes <- read_excel("C:/Users/gomes/OneDrive/Documents/PhD/GO_enrichment/gene_list_binary.xlsx")

# Ensure the column names match the expected format
str(head(all_genes))  # Expected columns: "gene_ID" and "significant"

# Prepare a named binary vector (1 for significant, 0 for not significant)
all_genes_named <- setNames(all_genes$significant, all_genes$gene_ID)
str(all_genes_named)

# Create a selection function to choose significant genes (1 = significant)
topDiffGenes <- function(allScore) {
  return(allScore == 1)  # Select significant genes based on binary presence/absence
}

# Create the topGOdata object
GOdata <- new("topGOdata",
              description = "GO enrichment",
              ontology = "BP",  # Change to "MF" or "CC" if needed
              allGenes = all_genes_named,  # Named vector with 1/0
              geneSel = topDiffGenes,  # Selection function for binary data
              nodeSize = 5,  # Minimum GO term size
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)  # Background GO mapping

# Run GO enrichment analysis using Fisher's exact test (binary data compatible)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher_01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

  # Generate a summary table with the top 20 enriched GO terms
  allRes <- GenTable(GOdata,
                     classicFisher = resultFisher,
                     elimFisher = resultElim,
                     fisher01 = resultFisher_01,
                     orderBy = "fisher01",
                     ranksOf = "elimFisher",
                     topNodes = 200)

# Display the results
print(allRes)

# Save the results to a CSV file
write_xlsx(allRes, "C:/Users/gomes/OneDrive/Documents/PhD/GO_enrichment/topGO_results_binary.xlsx")

# Visualize significant GO terms
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = "def")

# Save the graph to a PDF file
printGraph(GOdata, resultFisher, firstSigNodes = 15, 
           fn.prefix = "C:/Users/gomes/OneDrive/Documents/PhD/GO_enrichment/GO_graph_binary",
           useInfo = "def", pdfSW = TRUE)


GO_terms <- c("GO:0009698", "GO:0009699", "GO:0009809")
ann.genes <- genesInTerm(GOdata,GO_terms)
ann.genes_total <- genesInTerm(GOdata) 


ann.genes.df <- do.call(rbind, lapply(names(ann.genes_total), function(go_id) {
  data.frame(GO_ID = go_id, Gene_ID = ann.genes_total[[go_id]], stringsAsFactors = FALSE)
}))

# Save the data frame to an Excel file
write_xlsx(ann.genes.df, "C:/Users/gomes/OneDrive/Documents/PhD/GO_enrichment/GO_genes.xlsx")
