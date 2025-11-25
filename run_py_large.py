
import pandas as pd
import numpy as np
import time
from iq import preprocess, create_protein_list, create_protein_table

start_time = time.time()

# Read data
print("Reading data...")
d = pd.read_csv("large_test_data.tsv", sep="\t")

# Preprocess
print("Preprocessing...")
norm_data = preprocess(d, pdf_out=None, show_boxplot=False)

# Create protein list
print("Creating protein list...")
p_list = create_protein_list(norm_data)

# Create protein table
print("Quantifying...")
res = create_protein_table(p_list, method="maxLFQ")

end_time = time.time()
print(f"Total time: {end_time - start_time:.2f} seconds")

# Write result
res['estimate'].to_csv("py_result_large.tsv", sep="\t")
