
import pandas as pd
import numpy as np
import sys
from scipy.stats import pearsonr

def load_data(filename):
    try:
        # Try reading with index_col=0
        df = pd.read_csv(filename, sep="\t", index_col=0)
        return df
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        sys.exit(1)

def compare_r_py(r_file, py_file, tolerance=1e-8):
    print(f"Comparing {r_file} and {py_file}...")
    df_r = load_data(r_file)
    df_py = load_data(py_file)
    
    # Align
    common_index = df_r.index.intersection(df_py.index)
    common_cols = df_r.columns.intersection(df_py.columns)
    
    df_r = df_r.loc[common_index, common_cols]
    df_py = df_py.loc[common_index, common_cols]
    
    diff = np.abs(df_r.values - df_py.values)
    max_diff = np.nanmax(diff)
    
    print(f"Max difference (R vs Py): {max_diff}")
    
    if max_diff < tolerance:
        print("PASS: R and Python results match.")
    else:
        print("FAIL: R and Python results differ significantly.")
        # Check if it's just a few
        n_diff = np.sum(diff > tolerance)
        print(f"Number of differing values: {n_diff} / {diff.size}")

def compare_truth(est_file, truth_file):
    print(f"Comparing {est_file} and {truth_file}...")
    df_est = load_data(est_file)
    df_truth = load_data(truth_file)
    
    # Align
    common_index = df_est.index.intersection(df_truth.index)
    common_cols = df_est.columns.intersection(df_truth.columns)
    
    df_est = df_est.loc[common_index, common_cols]
    df_truth = df_truth.loc[common_index, common_cols]
    
    # Calculate correlation per protein
    correlations = []
    rmses = []
    
    for protein in common_index:
        est = df_est.loc[protein].values
        truth = df_truth.loc[protein].values
        
        # Filter NaNs
        mask = ~np.isnan(est) & ~np.isnan(truth)
        if np.sum(mask) < 3:
            continue
            
        est = est[mask]
        truth = truth[mask]
        
        # Correlation
        if np.std(est) > 0 and np.std(truth) > 0:
            corr, _ = pearsonr(est, truth)
            correlations.append(corr)
        
        # RMSD after centering (since absolute values might differ)
        est_centered = est - np.mean(est)
        truth_centered = truth - np.mean(truth)
        rmsd = np.sqrt(np.mean((est_centered - truth_centered)**2))
        rmses.append(rmsd)
        
    avg_corr = np.mean(correlations)
    avg_rmsd = np.mean(rmses)
    
    print(f"Average Pearson Correlation: {avg_corr:.4f}")
    print(f"Average RMSD (centered): {avg_rmsd:.4f}")
    
    if avg_corr > 0.9:
        print("PASS: High correlation with ground truth.")
    else:
        print("FAIL: Low correlation with ground truth.")

if __name__ == "__main__":
    compare_r_py("r_result_large.tsv", "py_result_large.tsv")
    print("-" * 20)
    compare_truth("py_result_large.tsv", "ground_truth.tsv")
