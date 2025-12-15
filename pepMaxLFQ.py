#!/usr/bin/env python3
"""pepMaxLFQ.py

Compute MaxLFQ-style protein intensities from a MaxQuant-style peptide table.

This is a self-contained, robust implementation intended for interactive use
and small-to-medium datasets. It implements the important algorithmic
aspects from Cox et al. (2014): pairwise medians, minimum shared-peptide
thresholds, large-ratio stabilization, weighting by shared-peptide counts,
and rescaling of solved profiles to preserve summed peptide intensity.

Usage:
  python iq/pepMaxLFQ.py peptides.txt

Output:
  <peptides.txt>.proteinMaxLFQ.tsv  (tab-separated, log2 intensities)

Dependencies: pandas, numpy, scipy
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse import csgraph
from scipy.optimize import least_squares


def delayed_normalization(run_matrix, sample_map, max_iter=200):
	"""
	Lightweight delayed-normalization solver.
	- run_matrix: DataFrame rows=features, cols=runs (log2 intensities)
	- sample_map: dict sample_name -> list of run column names
	Returns dict run_name -> additive log2 normalization factor
	"""
	runs = list(run_matrix.columns)
	run_idx = {r: i for i, r in enumerate(runs)}
	sample_names = list(sample_map.keys())
	sample_runs = [sample_map[s] for s in sample_names]
	X = run_matrix.values

	pairs = [(a, b) for a in range(len(sample_names) - 1) for b in range(a + 1, len(sample_names))]

	def residuals(n_log2):
		N = 2.0 ** n_log2
		s_per_sample = []
		for runs_group in sample_runs:
			idxs = [run_idx[r] for r in runs_group]
			vals = np.nansum((2.0 ** X[:, idxs]) * N[np.newaxis, idxs], axis=1)
			s_per_sample.append(vals)

		res = []
		for (a, b) in pairs:
			sa = s_per_sample[a]
			sb = s_per_sample[b]
			valid = (sa > 0) & (sb > 0)
			if not np.any(valid):
				continue
			diff = np.log2(sb[valid]) - np.log2(sa[valid])
			res.extend(diff.tolist())
		return np.array(res)

	x0 = np.zeros(len(runs))
	sol = least_squares(residuals, x0, max_nfev=max_iter)
	return dict(zip(runs, sol.x))


def maxLFQ_do(X, min_pair_count=2, use_large_ratio_stabilization=True):
	"""
	Compute MaxLFQ log2-profile for a connected set of samples.
	X: numpy array (peptides x samples) in log2-space (NaN for missing)
	Returns: log2 intensities (length = n_samples)
	"""
	# handle trivial single-column case
	if X.shape[1] == 1:
		col = X[:, 0]
		if np.all(np.isnan(col)):
			return np.array([np.nan])
		return np.array([np.nanmedian(col)])

	N = X.shape[1]
	AtA = np.zeros((N, N), dtype=float)
	Atb = np.zeros(N, dtype=float)

	counts = np.sum(~np.isnan(X), axis=0).astype(float)
	sums_linear = np.nansum(2.0 ** X, axis=0)

	for i in range(N - 1):
		for j in range(i + 1, N):
			mask = (~np.isnan(X[:, i])) & (~np.isnan(X[:, j]))
			n_shared = int(mask.sum())
			if n_shared < min_pair_count:
				continue

			# median of log2 ratios (j - i)
			try:
				rm = np.nanmedian(X[mask, j] - X[mask, i])
			except IndexError:
				continue

			# large-ratio stabilization
			if use_large_ratio_stabilization:
				N_max = max(int(counts[i]), int(counts[j]))
				x = float(N_max) / float(n_shared) if n_shared > 0 else np.inf

				if x <= 2.5:
					r_ij = rm
				elif x >= 5.0:
					si = sums_linear[i] if not np.isnan(sums_linear[i]) else 0.0
					sj = sums_linear[j] if not np.isnan(sums_linear[j]) else 0.0
					if si == 0 or sj == 0:
						r_ij = rm
					else:
						r_ij = np.log2(sj / si)
				else:
					si = sums_linear[i] if not np.isnan(sums_linear[i]) else 0.0
					sj = sums_linear[j] if not np.isnan(sums_linear[j]) else 0.0
					if si == 0 or sj == 0:
						r_ij = rm
					else:
						rs = np.log2(sj / si)
						w = (x - 2.5) / 2.5
						r_ij = w * rs + (1.0 - w) * rm
			else:
				r_ij = rm

			if np.isnan(r_ij):
				continue

			w = float(n_shared)
			AtA[i, j] -= w
			AtA[j, i] -= w
			AtA[i, i] += w
			AtA[j, j] += w

			Atb[i] -= w * r_ij
			Atb[j] += w * r_ij

	# If no valid pairs were added, fall back to per-column medians
	if not np.any(AtA):
		out = np.full(N, np.nan)
		for k in range(N):
			col = X[:, k]
			if np.all(np.isnan(col)):
				out[k] = np.nan
			else:
				out[k] = np.nanmedian(col)
		return out

	# Augment system with sum-to-zero constraint to make the laplacian invertible
	A = AtA.copy()
	b = Atb.copy()
	A2 = np.zeros((N + 1, N + 1), dtype=float)
	A2[:N, :N] = A
	A2[:N, N] = 1.0
	A2[N, :N] = 1.0
	b2 = np.zeros(N + 1, dtype=float)
	b2[:N] = b

	# solve (use least squares fallback for numerical safety)
	try:
		sol = np.linalg.solve(A2, b2)
	except np.linalg.LinAlgError:
		sol, *_ = np.linalg.lstsq(A2, b2, rcond=None)

	x = sol[:N]

	# rescale so that total linear abundance matches the sum of per-sample observed sums
	total_target = np.nansum(sums_linear)
	total_est = np.nansum(2.0 ** x)
	if total_est > 0 and not np.isnan(total_target) and total_target > 0:
		scale = total_target / total_est
		x = x + np.log2(scale)

	return x


def maxLFQ(sub_df, min_pair_count=2, use_large_ratio_stabilization=True):
	"""
	Compute MaxLFQ profile for a pandas DataFrame `sub_df` of peptides x samples
	where values are already log2-transformed intensities. Returns a dict with
	'estimate' (numpy array, length = n_samples) and 'annotation' (string).
	"""
	if sub_df.shape[1] == 0:
		return {'estimate': np.array([]), 'annotation': 'NA'}

	X = sub_df.values.astype(float)

	mask = ~np.isnan(X)
	adj = (mask.T @ mask) > 0
	n_components, labels = csgraph.connected_components(csr_matrix(adj), directed=False)

	N = X.shape[1]
	out = np.full(N, np.nan)

	for comp in range(n_components):
		ind = np.where(labels == comp)[0]
		if len(ind) == 1:
			col = X[:, ind[0]]
			if np.all(np.isnan(col)):
				out[ind] = np.nan
			else:
				out[ind] = np.nanmedian(col)
		else:
			sub = X[:, ind]
			out[ind] = maxLFQ_do(sub, min_pair_count=min_pair_count,
								 use_large_ratio_stabilization=use_large_ratio_stabilization)

	if np.all(np.isnan(out)):
		return {'estimate': out, 'annotation': 'NA'}

	quantified = ~np.isnan(out)
	q_labels = labels[quantified]
	if len(q_labels) > 0 and np.all(q_labels == q_labels[0]):
		return {'estimate': out, 'annotation': ''}
	else:
		anno = labels.astype(object)
		anno[~quantified] = np.nan
		return {'estimate': out, 'annotation': ';'.join(str(x) for x in anno)}


def main():
	if len(sys.argv) < 2:
		print('Usage: python iq/pepMaxLFQ.py peptides.txt')
		sys.exit(1)

	pathFiles = sys.argv[1]
	df = pd.read_csv(pathFiles, low_memory=False, sep='\t')

	# normalize column names we expect
	if 'Leading razor protein' in df.columns:
		df = df.rename(columns={'Leading razor protein': 'ID'})

	if 'ID' not in df.columns:
		raise ValueError("Expected column 'Leading razor protein' or 'ID' in peptide table")

	# explode protein IDs
	df['IDs'] = df['ID'].astype(str).str.split(';')
	dfE = df.explode('IDs').copy()

	# ensure Sequence exists for a stable index
	if 'Sequence' in dfE.columns:
		dfE.index = dfE['Sequence']

	# replace 0 with NaN for intensities
	dfE = dfE.replace(0, np.nan)

	# locate intensity columns (MaxQuant style: 'Intensity <sample>')
	int_cols = [c for c in dfE.columns if c.startswith('Intensity')]
	if len(int_cols) == 0:
		# fallback: any numeric columns other than known metadata
		numeric = dfE.select_dtypes(include=[np.number]).columns.tolist()
		int_cols = numeric

	# strip 'Intensity ' prefix for nicer sample names
	sample_names = [c.removeprefix('Intensity ').strip() for c in int_cols]

	dfE_int = dfE[int_cols].copy()
	dfE_int.columns = sample_names
	dfE_int = dfE_int.astype(float)

	# log2 transform
	with np.errstate(divide='ignore', invalid='ignore'):
		dfE_log2 = np.log2(dfE_int)

	print('Computing MaxLFQ per protein...')

	# --- Delayed normalization: estimate additive run factors and apply
	# Build a sample->runs map. By default each sample is a single run, but
	# this supports grouping if your intensity columns encode runs as
	# "Sample_Run" and you want to group by Sample.
	sample_map = {s: [s] for s in sample_names}
	try:
		print('Estimating delayed-normalization factors...')
		run_factors = delayed_normalization(dfE_log2, sample_map)
		# apply additive log2 corrections (subtract per-run factor)
		for run, f in run_factors.items():
			if run in dfE_log2.columns:
				dfE_log2[run] = dfE_log2[run] - float(f)
		print('Applied delayed-normalization factors.')
	except Exception as e:
		print('Delayed-normalization failed, proceeding without it:', e)

	protein_ids = pd.Index(dfE['IDs'].unique())
	samples = sample_names

	prot_tab = pd.DataFrame(index=protein_ids, columns=samples, dtype=float)
	prot_anno = {}

	for pid in protein_ids:
		mask = dfE['IDs'] == pid
		sub = dfE_log2[mask]
		sub = sub.dropna(how='all')
		if sub.shape[0] == 0:
			prot_tab.loc[pid] = np.nan
			prot_anno[pid] = 'NA'
			continue

		out = maxLFQ(sub)
		est = out.get('estimate')
		anno = out.get('annotation', '')

		try:
			prot_tab.loc[pid] = est
		except Exception:
			prot_tab.loc[pid] = np.asarray(est)

		prot_anno[pid] = anno

	out_fname = pathFiles + '.proteinMaxLFQ.tsv'
	prot_tab.to_csv(out_fname, sep='\t', na_rep='NA')
	print(f'Wrote MaxLFQ results to {out_fname}')


if __name__ == '__main__':
	main()

