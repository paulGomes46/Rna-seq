# Load necessary libraries
library(DESeq2)
library(tximport)
library(readxl)
library(writexl)

# Read tx2gene mapping (already done)
tx2gene <- read_excel("C:/Users/gomes/OneDrive/Documents/PhD/quantifications/convertion_transcript_gene.xlsx")

# Define the base directory for each sample
base_dir <- "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/my_data"

sample <- c("L_N_1A", "L_N_2A", "L_N_3A", 
            "L_P_1A", "L_P_2A", "L_P_3A", 
            "S_N_1A", "S_N_2A", "S_N_3A", 
            "S_P_1A", "S_P_2A", "S_P_3A")
dirs <- file.path(base_dir, sample)

condition <- c("L_N", "L_N", "L_N",
               "L_P", "L_P", "L_P",
               "S_N", "S_N", "S_N",
               "S_P", "S_P", "S_P")

# Step 3: Construct file paths for each sample's abundance.h5 file
files <- file.path(dirs, "abundance.h5")
names(files) <- sample

# Step 4: Import Kallisto data using tximport
txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene = tx2gene)

# Check the first few rows of the counts data
head(txi.kallisto$counts)

# Step 5: Create sampleTable with experimental design (condition of each sample)
sampleTable <- data.frame(condition = condition)
rownames(sampleTable) <- colnames(txi.kallisto$counts)  # Ensure the rownames match the sample names


############ VERIFICATION
# Verify that sampleTable matches txi.kallisto$counts
print("Checking sampleTable alignment...")
print(rownames(sampleTable))
print(colnames(txi.kallisto$counts))

# If they match perfectly, you're good to go!
if (all(rownames(sampleTable) == colnames(txi.kallisto$counts))) {
  print("✅ SampleTable is correctly aligned!")
} else {
  print("⚠️ WARNING: SampleTable does NOT match txi.kallisto$counts! Check your assignments.")
}


# Step 6: Create DESeqDataSet from tximport object and sampleTable
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

# Step 7: Perform differential expression analysis
dds <- DESeq(dds)


# Perform pairwise comparisons (you can add more as needed)
res_1 <- results(dds, contrast = c("condition", "L_P", "L_N"))  # A vs B
res_2 <- results(dds, contrast = c("condition", "S_P", "S_N"))  # A vs C
res_3 <- results(dds, contrast = c("condition", "L_P", "S_P"))  # A vs D
res_4 <- results(dds, contrast = c("condition", "L_N", "S_N"))  # B vs C


# Export results for each comparison, including gene IDs
res_1_df <- as.data.frame(res_1)
res_1_df$gene_id <- rownames(res_1_df)  # Add the gene IDs to a new column

res_2_df <- as.data.frame(res_2)
res_2_df$gene_id <- rownames(res_2_df)

res_3_df <- as.data.frame(res_3)
res_3_df$gene_id <- rownames(res_3_df)

res_4_df <- as.data.frame(res_4)
res_4_df$gene_id <- rownames(res_4_df)


# Export the results to Excel, including the gene IDs
write_xlsx(list(
  "L_P vs L_N" = res_1_df,
  "S_P vs S_N" = res_2_df,
  "L_P vs S_P" = res_3_df,
  "L_N vs S_N" = res_4_df
), path = "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/my_data/DESeq2_results_pairwise_comparisons.xlsx")


##### extract raw data
# Extract the counts from txi.kallisto
expression_data <- as.data.frame(txi.kallisto$counts)

# Optional: Add gene IDs as a column (if you want to make it explicit)
expression_data$gene_id <- rownames(expression_data)

# export to excel
library(writexl)
write_xlsx(expression_data, path = "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/my_data/raw_data_my_data.xlsx")


##### extract normalized data
# Extract normalized counts
foo <- counts(dds, normalized = TRUE)

# Convert to data frame and add gene IDs as a column (optional)
foo_df <- as.data.frame(foo)
foo_df$gene_id <- rownames(foo_df)

# Export to CSV, including gene IDs
write_xlsx(foo_df, path="C:/Users/gomes/OneDrive/Documents/PhD/quantifications/my_data/normalized_data_my_data.xlsx")



