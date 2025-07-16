# Load necessary libraries
#BiocManager::install("edgeR")

library(edgeR)
library(dplyr)
library(ggplot2)

# Define sample metadata
sample_names <- c("OE1", "OE2", "OE3", "WT1", "WT2", "WT3")
condition <- factor(c("OE", "OE", "OE", "WT", "WT", "WT"))  # Experimental conditions
files <- c(
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467110_quantif/abundance.tsv",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467111_quantif/abundance.tsv",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467112_quantif/abundance.tsv",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467113_quantif/abundance.tsv",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467114_quantif/abundance.tsv",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467115_quantif/abundance.tsv"
)

# Read in abundance files (Extract transcript IDs and estimated counts)
read_abundance <- function(file) {
  data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  return(data[, c("target_id", "est_counts")])  # Select only relevant columns
}

# Read all abundance.tsv files
count_data_list <- lapply(files, read_abundance)

# Merge into a single matrix
count_data <- Reduce(function(x, y) merge(x, y, by = "target_id"), count_data_list)
colnames(count_data) <- c("GeneID", sample_names)  # Rename columns

# Convert to matrix format for EdgeR
rownames(count_data) <- count_data$GeneID
count_data <- count_data[, -1]  # Remove GeneID column

# Convert to DGEList object
dge <- DGEList(counts = as.matrix(count_data), group = condition)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
dge <- calcNormFactors(dge)

# Create experimental design matrix
design <- model.matrix(~ condition)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmQLFit(dge, design)

# Perform differential expression analysis
qlf <- glmQLFTest(fit, coef = 2)  # coef=2 compares OE vs WT

# Get results table
results <- topTags(qlf, n = Inf)$table

# Adjusted p-values (False Discovery Rate)
results <- results %>% arrange(FDR)

# Save DEGs to CSV file
write.csv(results, "DEG_results_EdgeR.csv", row.names = TRUE)

# Filter significant DEGs (FDR < 0.05)
significant_DEGs <- results %>% filter(FDR < 0.05)

# Save significant DEGs
write.csv(significant_DEGs, "Significant_DEGs_EdgeR.csv", row.names = TRUE)

# Print top DEGs
head(significant_DEGs, 20)

# Create a Volcano Plot
ggplot(results, aes(x = logFC, y = -log10(FDR))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  theme_minimal() +
  labs(title = "Volcano Plot - EdgeR", x = "Log2 Fold Change", y = "-log10(FDR)")
