
import pandas as pd
import numpy as np
from iq import preprocess, create_protein_list, create_protein_table

# Read data
d = pd.read_csv("test_data.tsv", sep="\t")

# Preprocess
# Note: pandas reads "NA" as NaN.
# iq.py preprocess expects numeric column.
norm_data = preprocess(d, pdf_out=None, show_boxplot=False)

# Create protein list
p_list = create_protein_list(norm_data)

# Create protein table
res = create_protein_table(p_list, method="maxLFQ")

# Write result
res['estimate'].to_csv("py_result.tsv", sep="\t")
