
source("R/iq.R")

start_time <- Sys.time()

# Read data
message("Reading data...")
d <- read.delim("large_test_data.tsv", stringsAsFactors = FALSE)

# Preprocess
message("Preprocessing...")
norm_data <- preprocess(d, pdf_out = NULL, show_boxplot = FALSE)

# Create protein list
message("Creating protein list...")
p_list <- create_protein_list(norm_data)

# Create protein table
message("Quantifying...")
res <- create_protein_table(p_list, method = "maxLFQ")

end_time <- Sys.time()
message("Total time: ", end_time - start_time)

# Write result
write.table(res$estimate, "r_result_large.tsv", sep = "\t", quote = FALSE)
