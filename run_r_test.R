
source("R/iq.R")

# Read data
d <- read.delim("test_data.tsv", stringsAsFactors = FALSE)

# Preprocess
# Note: default na_string is "0". But our file has "NA".
# read.delim reads "NA" as NA value in R.
# iq.R preprocess checks:
# if (is.null(intensity_col_sep)) { if (!is.numeric(...)) ... }
# If F.PeakArea is read as numeric (with NAs), it enters the first block.
# In the first block (lines 34-51), it does:
# quant = log2(quant_table[, intensity_col])
# This handles NAs correctly (log2(NA) is NA).

norm_data <- preprocess(d, pdf_out = NULL)

# Create protein list
p_list <- create_protein_list(norm_data)

# Create protein table
res <- create_protein_table(p_list, method = "maxLFQ")

# Write result
write.table(res$estimate, "r_result.tsv", sep = "\t", quote = FALSE)
