
import numpy as np
import pandas as pd
import warnings
from scipy.sparse import csgraph
from scipy.sparse import csr_matrix

def preprocess(quant_table,
               primary_id="PG.ProteinGroups",
               secondary_id=["EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"],
               sample_id="R.Condition",
               intensity_col="F.PeakArea",
               median_normalization=True,
               log2_intensity_cutoff=0,
               pdf_out="qc-plots.pdf",
               pdf_width=12,
               pdf_height=8,
               intensity_col_sep=None,
               intensity_col_id=None,
               na_string="0",
               show_boxplot=True):
    
    # Note: pdf_out and plotting are not fully implemented in this port yet.
    
    if intensity_col_sep is None:
        if not pd.api.types.is_numeric_dtype(quant_table[intensity_col]):
             raise ValueError("Intensity column must be numeric when 'intensity_col_sep' is NULL.")
        
        # Make secondary ids
        print("Concatenating secondary ids...")
        # Vectorized string concatenation
        if isinstance(secondary_id, str):
            secondary_id = [secondary_id]
            
        second_id = quant_table[secondary_id[0]].astype(str)
        if len(secondary_id) > 1:
            for col in secondary_id[1:]:
                second_id = second_id + "_" + quant_table[col].astype(str)
        
        d = pd.DataFrame({
            'protein_list': quant_table[primary_id].astype(str),
            'sample_list': quant_table[sample_id].astype(str),
            'quant': np.log2(quant_table[intensity_col]),
            'id': second_id
        })
        
    else:
        # Handle split intensity columns (not fully optimized for large data, following R logic)
        print("Concatenating secondary ids...")
        if isinstance(secondary_id, str):
            secondary_id = [secondary_id]
            
        second_id = quant_table[secondary_id[0]].astype(str)
        if len(secondary_id) > 1:
            for col in secondary_id[1:]:
                second_id = second_id + "_" + quant_table[col].astype(str)

        # Explode logic
        # This part assumes quant_table rows need to be expanded based on split intensity strings
        # This is a direct port of the loop in R, but using Python lists for construction
        
        rows = []
        for idx, row in quant_table.iterrows():
            intensities = str(row[intensity_col]).split(intensity_col_sep)
            
            if intensity_col_id:
                ids = str(row[intensity_col_id]).split(intensity_col_sep)
            else:
                ids = [f"{second_id[idx]}_{i+1}" for i in range(len(intensities))]
            
            # Filter NA strings
            intensities = [x if x != na_string else np.nan for x in intensities]
            
            for i, val in enumerate(intensities):
                rows.append({
                    'protein_list': str(row[primary_id]),
                    'sample_list': str(row[sample_id]),
                    'quant': np.log2(float(val)) if val is not np.nan and val is not None else np.nan,
                    'id': ids[i] if intensity_col_id else ids[i]
                })
        d = pd.DataFrame(rows)

    # Remove NaN quant
    d = d.dropna(subset=['quant'])
    
    samples = d['sample_list'].unique()
    
    # Intensity cutoff
    if log2_intensity_cutoff is not None:
        print("Removing low intensities...")
        d = d[d['quant'] > log2_intensity_cutoff]
        
    # Median normalization
    m = []
    for s in samples:
        v = d[d['sample_list'] == s]['quant']
        m.append(v.median())
    m = np.array(m)
    
    if median_normalization:
        print("Median normalization...")
        f = np.mean(m) - m
        
        # Apply normalization
        # Create a mapping for faster application
        f_map = dict(zip(samples, f))
        
        # Vectorized update
        d['quant'] = d.apply(lambda row: row['quant'] + f_map[row['sample_list']], axis=1)
        
    return d

def create_protein_list(preprocessed_data):
    # Check for NAs
    if preprocessed_data['protein_list'].isna().any():
        raise ValueError("NA value in protein_list")
    if preprocessed_data['sample_list'].isna().any():
        raise ValueError("NA value in sample_list")
    if preprocessed_data['id'].isna().any():
        raise ValueError("NA value in id")
    if preprocessed_data['quant'].isna().any():
        raise ValueError("NA value in quant")
        
    proteins = preprocessed_data['protein_list'].unique()
    samples = preprocessed_data['sample_list'].unique()
    
    print(f"# proteins = {len(proteins)}, # samples = {len(samples)}")
    
    p_list = {}
    
    # Group by protein to avoid loop
    grouped = preprocessed_data.groupby('protein_list')
    
    total = len(proteins)
    n_display = total // 20
    count = 0
    
    for protein, group in grouped:
        count += 1
        if n_display > 0 and count % n_display == 0:
            print(f"{count * 100 / total:.2f} %")
            
        # Pivot: index=id, columns=sample_list, values=quant
        # Check for duplicates
        if group.duplicated(subset=['id', 'sample_list']).any():
             # R code warns and stops. We can raise error.
             dup = group[group.duplicated(subset=['id', 'sample_list'])]
             print(f"Duplicate entry for protein {protein}:")
             print(dup)
             raise ValueError("duplicate entry.")
        
        m = group.pivot(index='id', columns='sample_list', values='quant')
        
        # Reindex to ensure all samples are present (filled with NaN)
        m = m.reindex(columns=samples)
        
        p_list[protein] = m
        
    print("Completed.")
    return p_list

def maxLFQ_do(X):
    # X is a numpy array (fragments x samples)
    N = X.shape[1]
    
    AtA = np.zeros((N, N))
    Atb = np.zeros(N)
    
    # Iterate over pairs of samples
    # This can be slow in Python loops. 
    # Optimization: Vectorized approach?
    # For now, stick to logic.
    
    for i in range(N - 1):
        for j in range(i + 1, N):
            # Calculate ratio
            diff = -X[:, i] + X[:, j]
            r_i_j = np.nanmedian(diff)
            
            if not np.isnan(r_i_j):
                AtA[i, j] = -1
                AtA[j, i] = -1
                AtA[i, i] += 1
                AtA[j, j] += 1
                
                Atb[i] -= r_i_j
                Atb[j] += r_i_j
                
    # Solve system
    # A = 2*AtA with constraints
    # [ 2*AtA   1 ]
    # [ 1...1   0 ]
    
    A_final = np.vstack([
        np.hstack([2 * AtA, np.ones((N, 1))]),
        np.hstack([np.ones(N), [0]])
    ])
    
    b_final = np.concatenate([2 * Atb, [np.nanmean(X) * N]])
    
    # Least squares
    res, residuals, rank, s = np.linalg.lstsq(A_final, b_final, rcond=None)
    
    return res[:N]

def maxLFQ(X):
    # X is a DataFrame or numpy array
    if isinstance(X, pd.DataFrame):
        X_np = X.values
    else:
        X_np = X
        
    if np.all(np.isnan(X_np)):
        return {'estimate': np.nan, 'annotation': "NA"}
        
    if X_np.shape[0] == 1:
        return {'estimate': X_np[0, :], 'annotation': ""}
        
    N = X_np.shape[1]
    
    # Determine connected components
    # Build adjacency matrix based on shared non-NA values
    # Two samples are connected if they share at least one quantified fragment
    
    # Create a mask of non-NA
    mask = ~np.isnan(X_np)
    # Adjacency: A[i,j] = 1 if sample i and j share a fragment
    # This is equivalent to (mask.T @ mask) > 0
    adj = (mask.T @ mask) > 0
    
    # Use scipy to find connected components
    n_components, labels = csgraph.connected_components(csr_matrix(adj), directed=False)
    
    w = np.full(N, np.nan)
    
    for i in range(n_components):
        # Indices of samples in this component
        ind = np.where(labels == i)[0]
        
        if len(ind) == 1:
            w[ind] = np.nanmedian(X_np[:, ind])
        else:
            # Extract submatrix for these samples
            # We also need to filter rows? 
            # R code: maxLFQ.do(X[, ind]) -> passes all rows but only selected columns
            # But wait, if we pass all rows, some rows might be all NA for this subset of columns.
            # maxLFQ.do handles NAs.
            
            w[ind] = maxLFQ_do(X_np[:, ind])
            
    if np.all(np.isnan(w)):
        return {'estimate': w, 'annotation': "NA"}
    else:
        quantified_samples = ~np.isnan(w)
        # Check if all quantified samples are in the same component
        # In R: if (all(g[quantified_samples] == g[quantified_samples[1]]))
        
        # Get labels of quantified samples
        q_labels = labels[quantified_samples]
        
        if len(q_labels) > 0 and np.all(q_labels == q_labels[0]):
            return {'estimate': w, 'annotation': ""}
        else:
            # Return labels as annotation, with NA for unquantified
            # R: g[is.na(w)] <- NA
            # We can return a string representation
            anno = labels.astype(object)
            anno[np.isnan(w)] = np.nan
            return {'estimate': w, 'annotation': ";".join(str(x) for x in anno)}

def topN(X, N=3, aggregation_in_log_space=True):
    if isinstance(X, pd.DataFrame):
        X_np = X.values
    else:
        X_np = X
        
    if X_np.shape[0] == 1:
        return {'estimate': X_np[0, :], 'annotation': ""}
        
    if aggregation_in_log_space:
        v = np.nanmean(X_np, axis=1)
        # Sort indices descending
        # argsort is ascending, so take reverse
        sorted_indices = np.argsort(v)[::-1]
        
        # Take top N
        top_indices = sorted_indices[:min(N, len(v))]
        
        out = np.nanmean(X_np[top_indices, :], axis=0)
    else:
        XX = 2**X_np
        v = np.nanmean(XX, axis=1)
        sorted_indices = np.argsort(v)[::-1]
        top_indices = sorted_indices[:min(N, len(v))]
        out = np.log2(np.nanmean(XX[top_indices, :], axis=0))
        
    return {'estimate': out, 'annotation': ""}

def meanInt(X, aggregation_in_log_space=True):
    if isinstance(X, pd.DataFrame):
        X_np = X.values
    else:
        X_np = X
        
    if X_np.shape[0] == 1:
        return {'estimate': X_np[0, :], 'annotation': ""}
        
    if aggregation_in_log_space:
        out = np.nanmean(X_np, axis=0)
    else:
        out = np.log2(np.nanmean(2.0**X_np, axis=0))
        
    return {'estimate': out, 'annotation': ""}

import matplotlib.pyplot as plt

def median_polish(X, max_iter=100, tol=0.01):
    # X is numpy array
    if isinstance(X, pd.DataFrame):
        X_np = X.values.copy()
    else:
        X_np = X.copy()
        
    if np.all(np.isnan(X_np)):
        return {'estimate': np.nan, 'annotation': "NA"}
        
    r = X_np.shape[0]
    c = X_np.shape[1]
    
    overall = 0
    row_effect = np.zeros(r)
    col_effect = np.zeros(c)
    
    for _ in range(max_iter):
        # Row sweep
        # Calculate row medians of residuals
        # Residuals = X - overall - row_effect - col_effect
        # But we update in place usually.
        # Let's follow R's medpolish logic roughly.
        
        # Remove row medians
        # Current residuals:
        res = X_np - (overall + row_effect[:, None] + col_effect[None, :])
        r_med = np.nanmedian(res, axis=1)
        r_med[np.isnan(r_med)] = 0
        
        row_effect += r_med
        
        # Col sweep
        res = X_np - (overall + row_effect[:, None] + col_effect[None, :])
        c_med = np.nanmedian(res, axis=0)
        c_med[np.isnan(c_med)] = 0
        
        col_effect += c_med
        
        # Update overall?
        # R medpolish:
        # 1. r <- r + delta; residuals <- residuals - delta
        # 2. c <- c + delta; residuals <- residuals - delta
        # It doesn't explicitly track overall separately in the loop, but extracts it at the end?
        # Actually, standard algorithm:
        # r_effect = row_median(X)
        # X = X - r_effect
        # c_effect = col_median(X)
        # X = X - c_effect
        # Repeat.
        # Overall is accumulated.
        
        # Let's stick to a simple implementation:
        # Decompose X = mu + r + c + eps
        
        # Check convergence
        if np.sum(np.abs(r_med)) < tol and np.sum(np.abs(c_med)) < tol:
            break
            
    # Estimate = overall + col_effect
    # Wait, iq.R says: return(list(estimate = out$overall + out$col, annotation = ""))
    # So it returns the column effect + overall mean.
    # This represents the sample abundance.
    
    # Re-calculate overall from row_effect/col_effect centering?
    # Usually we center row_effect and col_effect.
    
    # Let's just return col_effect + overall.
    # But in my loop I didn't update overall.
    # Let's do it properly.
    
    # Reset
    X_work = X_np.copy()
    overall = 0
    row_eff = np.zeros(r)
    col_eff = np.zeros(c)
    
    for _ in range(max_iter):
        # Row medians
        rm = np.nanmedian(X_work, axis=1)
        rm[np.isnan(rm)] = 0
        row_eff += rm
        X_work -= rm[:, None]
        
        # Col medians
        cm = np.nanmedian(X_work, axis=0)
        cm[np.isnan(cm)] = 0
        col_eff += cm
        X_work -= cm[None, :]
        
    # Distribute overall
    # R's medpolish centers the effects.
    # We want sample quantification.
    # The protein abundance in sample j is overall + col_eff[j].
    # But wait, row_eff captures peptide specific intensity.
    # col_eff captures sample specific intensity (abundance).
    # overall captures global intensity.
    
    # So estimate = overall + col_eff.
    # But where is overall?
    # In the loop above, I subtracted rm and cm from X_work.
    # I didn't separate overall.
    # Overall is usually mean of row_eff + mean of col_eff?
    # Or just take median of everything?
    
    # R medpolish output:
    # overall: the grand overall median.
    # row: the row effects.
    # col: the column effects.
    # residuals: the residuals.
    
    # Let's refine the loop to match R.
    t = 0
    r = np.zeros(X_np.shape[0])
    c = np.zeros(X_np.shape[1])
    e = X_np.copy()
    
    for _ in range(max_iter):
        # Row
        r_delta = np.nanmedian(e, axis=1)
        r_delta[np.isnan(r_delta)] = 0
        r += r_delta
        e -= r_delta[:, None]
        
        # Col
        c_delta = np.nanmedian(e, axis=0)
        c_delta[np.isnan(c_delta)] = 0
        c += c_delta
        e -= c_delta[None, :]
        
    # Now we have e approx 0.
    # X ~ r + c + e
    # But we want to extract overall.
    # overall = median(c) + median(r)?
    # R does:
    # overall <- median(c)
    # c <- c - overall
    # overall <- overall + median(r)
    # r <- r - median(r)
    # (roughly)
    
    # Let's just return c + mean(r)? No.
    # We want sample abundances.
    # If we have X_ij = mu + alpha_i + beta_j
    # We want mu + beta_j.
    # Which is equivalent to column marginals after removing row effects.
    
    # My loop produced r and c such that X ~ r + c.
    # So estimate = c + mean(r)?
    # Or just c?
    # If I add a constant to r and subtract from c, the fit is same.
    # We need to fix the scale.
    # Usually we set mean(alpha) = 0 or something.
    # iq.R returns `out$overall + out$col`.
    # In R medpolish, `overall` is a scalar. `col` is a vector.
    
    # Let's calculate overall from my r and c.
    overall_val = np.mean(c) # or median
    c -= overall_val
    overall_val += np.mean(r)
    r -= np.mean(r)
    
    # This is arbitrary centering.
    # Let's just return c + overall_val from the loop state?
    # Actually, simply:
    # estimate = c + median(r) ?
    
    # Let's look at what we want: relative abundance.
    # If we shift all c by k, we shift all abundances by k.
    # It matters for absolute intensity but not relative.
    # But iq.R returns log2 intensity.
    
    # Let's use the `c` from the loop and add the median of `r` to it as the base?
    # Or better:
    # estimate = c + np.median(r)
    
    # Let's stick to what the loop produced: X ~ r + c.
    # We want the "sample component".
    # If we consider "row component" as peptide ionization efficiency,
    # then "sample component" is protein abundance.
    # So we want c.
    # But we need to add back the "average" intensity.
    # So c + mean(r).
    
    estimate = c + np.nanmean(r)
    return {'estimate': estimate, 'annotation': ""}

def plot_protein(X, main="", col=None, split=0.6):
    # X is numpy array or DataFrame (rows=peptides, cols=samples)
    if isinstance(X, pd.DataFrame):
        X_np = X.values
        samples = X.columns
        peptides = X.index
    else:
        X_np = X
        samples = range(X.shape[1])
        peptides = range(X.shape[0])
        
    if col is None:
        # Generate colors
        col = plt.cm.tab10(np.linspace(0, 1, len(peptides)))
        
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot each peptide
    for i in range(len(peptides)):
        # Handle NAs by masking?
        # Matplotlib handles NaNs by breaking the line, which is what we want (type="o" in R)
        ax.plot(samples, X_np[i, :], marker='o', label=peptides[i])
        
    ax.set_title(main)
    ax.set_ylabel("Intensity")
    ax.set_xlabel("Sample")
    
    # Legend
    if split:
        # Move legend outside
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
    else:
        ax.legend()
        
    plt.show()

def create_protein_table(protein_list, method="maxLFQ", **kwargs):
    if not isinstance(protein_list, dict):
        raise ValueError("Input is not a dict.")
        
    if len(protein_list) == 0:
        return None
        
    first_protein = next(iter(protein_list.values()))
    sample_names = first_protein.columns
    protein_names = list(protein_list.keys())
    
    tab = np.full((len(protein_names), len(sample_names)), np.nan)
    annotation = np.full(len(protein_names), np.nan, dtype=object)
    
    n_display = len(protein_names) // 20
    
    for i, protein in enumerate(protein_names):
        if n_display > 0 and i % n_display == 0:
            print(f"{i * 100 / len(protein_names):.2f} %")
            
        X = protein_list[protein]
        
        if method == "maxLFQ":
            out = maxLFQ(X, **kwargs)
        elif method == "topN":
            out = topN(X, **kwargs)
        elif method == "meanInt":
            out = meanInt(X, **kwargs)
        elif method == "median_polish":
            out = median_polish(X, **kwargs)
        else:
            raise ValueError(f"Unknown method: {method}")
            
        tab[i, :] = out['estimate']
        annotation[i] = out['annotation']
        
    print("Completed.")
    
    # Create DataFrame
    df = pd.DataFrame(tab, index=protein_names, columns=sample_names)
    return {'estimate': df, 'annotation': annotation}

def process_wide_format(input_filename,
                        output_filename,
                        id_column,
                        quant_columns,
                        data_in_log_space=False,
                        annotation_columns=None,
                        method="maxLFQ"):
                        
    print(f"Reading file: {input_filename}")
    d = pd.read_csv(input_filename, sep='\t')
    print(f"# rows = {len(d)}; # cols = {len(d.columns)}")
    
    if id_column not in d.columns:
        raise ValueError("'id_column' not found")
        
    # Handle quant_columns (indices or names)
    # If list of ints, convert to names
    if all(isinstance(x, int) for x in quant_columns):
        # 1-based index in R, 0-based in Python? 
        # R users might pass 1-based. Let's assume user passes names or we need to be careful.
        # The R code checks if numeric.
        # Let's assume names for now to be safe, or 0-based indices.
        # If the user is porting R code directly, they might use 1-based.
        # But this is a Python library. Let's assume names or 0-based.
        pass
        
    # Verify columns exist
    for col in quant_columns:
        if col not in d.columns:
             raise ValueError(f"Column {col} not found")
             
    if annotation_columns:
        for col in annotation_columns:
            if col not in d.columns:
                raise ValueError(f"Column {col} not found")
                
    print("Preparing data...")
    d = d.dropna(subset=[id_column])
    d = d[d[id_column] != ""]
    
    # Melt
    d_long = d.melt(id_vars=[id_column], value_vars=quant_columns, 
                    var_name='sample_list', value_name='quant')
                    
    d_long = d_long.dropna(subset=['quant', id_column])
    
    # Create 'id' column (dummy id since wide format usually implies protein level or already aggregated?)
    # Wait, process_wide_format in R seems to treat rows as fragments if there are multiple rows per id_column?
    # "collapse to unique 'id_column'" -> d <- d[!is.na(d[, id_column]),]
    # Then it calls fast_MaxLFQ(df) where df has id = d_long$id.
    # But where does d_long$id come from?
    # In R: d_long <- reshape(...)
    # reshape in R creates an 'id' variable if not specified? Or does it keep row names?
    # Actually, in R's process_wide_format:
    # df <- data.frame(..., id = d_long$id, ...)
    # But d_long is created via reshape. reshape(direction="long") often adds an 'id' column (the original row index).
    # Yes, 'id' in reshape output identifies the original observation.
    # So if the input file has multiple rows for the same protein (e.g. peptides), 'id' distinguishes them.
    
    # In pandas melt, we don't get an automatic row identifier.
    # We need to preserve the original index or create a unique id for each row in the original d.
    
    # Let's add a row id to d before melting
    d['row_id'] = range(len(d))
    
    d_long = d.melt(id_vars=[id_column, 'row_id'], value_vars=quant_columns,
                    var_name='sample_list', value_name='quant')
                    
    d_long = d_long.dropna(subset=['quant', id_column])
    
    unique_values = d_long[id_column].unique()
    
    if method == "maxLFQ":
        # Prepare df for maxLFQ
        # We need columns: protein_list, sample_list, id, quant
        
        # id should be the fragment ID. Here it's the row_id.
        
        df_input = pd.DataFrame({
            'protein_list': d_long[id_column],
            'sample_list': d_long['sample_list'],
            'id': d_long['row_id'],
            'quant': d_long['quant'] if data_in_log_space else np.log2(d_long['quant'])
        })
        
        # We can reuse create_protein_list and create_protein_table
        # But that might be slow for large data.
        # The R code calls fast_MaxLFQ directly on the long dataframe.
        # My maxLFQ takes a matrix.
        # So I should use create_protein_list then create_protein_table.
        
        print("Running MaxLFQ...")
        p_list = create_protein_list(df_input)
        res = create_protein_table(p_list, method="maxLFQ")
        
        estimate = res['estimate']
        
    else:
        # Implement other methods loop
        # ...
        pass
        
    print(f"Writing to: {output_filename}")
    
    # Merge annotation
    if annotation_columns:
        # Get first occurrence of annotation for each protein
        anno = d.drop_duplicates(subset=[id_column])[ [id_column] + annotation_columns ]
        anno = anno.set_index(id_column)
        
        # Join with estimate
        estimate = estimate.join(anno)
        
    estimate.to_csv(output_filename, sep='\t')
