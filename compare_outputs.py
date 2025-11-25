
import pandas as pd
import numpy as np
import sys

def compare_files(r_file, py_file, tolerance=1e-10):
    print(f"Comparing {r_file} and {py_file}...")
    
    # Read R file
    # R write.table with row.names=T, col.names=T produces a header line with n cols and data with n+1 cols.
    # Pandas handles this if we don't specify header=0? Or we need to be careful.
    # Let's try reading with index_col=0.
    try:
        # If the header has one fewer field than data, pandas read_csv with index_col=0 usually aligns correctly.
        df_r = pd.read_csv(r_file, sep="\t", index_col=0)
    except Exception as e:
        print(f"Error reading R file: {e}")
        sys.exit(1)
        
    # Read Python file
    # Python to_csv with index=True produces a header with an empty first field if index name is None.
    # This matches the data column count.
    try:
        df_py = pd.read_csv(py_file, sep="\t", index_col=0)
    except Exception as e:
        print(f"Error reading Python file: {e}")
        sys.exit(1)
        
    print("R shape:", df_r.shape)
    print("Py shape:", df_py.shape)
    
    # Align columns and index
    # Ensure same columns and rows
    common_cols = df_r.columns.intersection(df_py.columns)
    common_index = df_r.index.intersection(df_py.index)
    
    if len(common_cols) != len(df_r.columns) or len(common_cols) != len(df_py.columns):
        print("Warning: Columns mismatch")
        print("R cols:", df_r.columns)
        print("Py cols:", df_py.columns)
        
    if len(common_index) != len(df_r.index) or len(common_index) != len(df_py.index):
        print("Warning: Index mismatch")
        print("R index:", df_r.index)
        print("Py index:", df_py.index)
        
    df_r = df_r.loc[common_index, common_cols]
    df_py = df_py.loc[common_index, common_cols]
    
    # Compare values
    diff = np.abs(df_r.values - df_py.values)
    max_diff = np.nanmax(diff)
    
    print(f"Max difference: {max_diff}")
    
    if max_diff < tolerance:
        print("SUCCESS: Results match within tolerance.")
    else:
        print("FAILURE: Results differ.")
        print("Differences:")
        print(diff)
        sys.exit(1)

if __name__ == "__main__":
    compare_files("r_result.tsv", "py_result.tsv")
