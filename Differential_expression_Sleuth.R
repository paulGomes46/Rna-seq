library(sleuth)

# Define your sample information
sample_id <- c("OE1", "OE2", "OE3", "wt1", "wt2", "wt3")
condition <- c("OE", "OE", "OE", "wt", "wt", "wt")
#path <- file.path(S_N_1A, S_N_2A, S_N_3A, S_P_1A, S_P_2A, S_P_3A)

path <- c(
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467110_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467111_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467112_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467113_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467114_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467115_quantif"
)

# Create the sample metadata table
sample_table <- data.frame(sample = sample_id, condition = condition, path = path)

# Prepare the data (normalization and log-transformation)
so <- sleuth_prep(sample_table, extra_bootstrap_summary = TRUE)

# Fit models
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')

# Perform differential expression analysis using Likelihood Ratio Test (LRT)
so <- sleuth_lrt(so, 'reduced', 'full')

# Check models
models(so)

# Extract results
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# Filter significant DEGs (qval ≤ 0.05)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

# Export DEGs to CSV file
write.csv(sleuth_significant, "DEG_Lv.csv", row.names = FALSE)

# Plot the abundance of a gene in the 6 samples
plot_bootstrap(so, "Cs5g_pb029960.1", units = "est_counts", color_by = "condition")