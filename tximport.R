#BiocManager::install("tximport")
install.packages("readxl")

library(tximport)
library(readxl)


tx2gene <- read_excel("C:/Users/gomes/OneDrive/Documents/PhD/quantifications/convertion_transcript_gene.xlsx")
head(tx2gene)



# Step 1: Define the base directory for each sample
dirs <- c(
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467110_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467111_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467112_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467113_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467114_quantif",
  "C:/Users/gomes/OneDrive/Documents/PhD/quantifications/Lv/SRR22467115_quantif"
)

# Step 2: Construct file paths for each sample's abundance.h5 file
files <- file.path(dirs, "abundance.h5")
names(files) <- paste0("sample", 1:6)

txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene=tx2gene)

head(txi.kallisto$counts)
