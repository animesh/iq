**pepMaxLFQ — Purpose & Quick Start**
- **Purpose:**: Compute protein-level quantification using the MaxLFQ-style algorithm from peptide intensities (log2 space). The script processes MaxQuant-style peptide output, aggregates peptides by protein, runs the inlined MaxLFQ implementation, and writes per-protein estimates.
- **Usage:**: Run the script from the `iq` directory (default input file is `peptides.txt`):

```bash
python pepMaxLFQ.py peptides.txt
```

**What the script expects (input)**
- **Input file:**: A tab-separated peptide table (example: MaxQuant `peptides.txt`) containing:
  - a column named `Leading razor protein` (the script renames this to `ID`),
  - peptide intensity columns with names matching `Intensity ...` (script finds these via regex),
  - a `Sequence` column (used as peptide index during explosion).
- **Format notes:**: The script explodes multi-protein peptide assignments by splitting the `Leading razor protein` field on semicolons, so one peptide row may be assigned to multiple protein IDs.

**High-level step-by-step behavior**
- **1) Read input file:**: loads the TSV into a pandas DataFrame.
- **2) Normalize protein IDs:**: renames `Leading razor protein` → `ID`, splits it on `;` and `explode()`s rows so each peptide maps to a single protein ID (column `IDs`).
- **3) Aggregate peptide intensities:**: creates several intermediate outputs:
  - `<input>.combinedIntensity.csv` — sum-aggregated intensities per protein (sum of intensities across member peptides).
  - `<input>peptidesSel.sumIntensity.csv` — example subset output saved for the selected protein used in code (keeps original behaviour).
- **4) Replace zeros with NaN:**: zeros are treated as missing (converted to `np.nan`) before log transformation.
- **5) Log2 transform:**: intensity columns are transformed to log2.
- **6) Compute MaxLFQ per protein:**:
  - For each protein, a peptide × sample matrix is built from log2 intensities (rows = peptides assigned to the protein, columns = samples).
  - If the protein has only one quantified peptide the peptide vector is used directly.
  - Otherwise the code runs an inlined MaxLFQ implementation (two functions: `maxLFQ_do` and `maxLFQ`).

    maxLFQ_do (core):
    - For each pair of samples (i, j) it computes the median of (log2 intensity_j - log2 intensity_i) across peptides that are quantified in both samples. This yields pairwise medians r_{i,j}.
    - It accumulates a Laplacian-like matrix `AtA` and right-hand side `Atb` from these pairwise medians.
    - An augmented linear system is constructed to fix the global additive constant, and the system is solved (least squares) to recover sample offsets (relative protein abundances).

    maxLFQ (wrapper):
    - Identifies connected components of the sample graph (samples connected if they share at least one quantified peptide).
    - Runs `maxLFQ_do` on each connected component to compute estimates for that subset of samples.
    - Returns a vector of sample estimates (log2 space) and a small annotation string when components are disconnected.

- **7) Write outputs:**: the script writes
  - `<input>.proteinMaxLFQ.tsv` — the main per-protein MaxLFQ table (rows = proteins, columns = samples, values in log2 space; `NA` for missing),
  - other intermediate CSVs created earlier remain in place (sample-by-sample diffs, medians, etc.).

**Files produced by the run (examples)**
- `peptides.txt.proteinMaxLFQ.tsv` — per-protein log2 MaxLFQ table (primary output).
- `peptides.txt.combinedIntensity.csv` — sum of peptide intensities per protein.
- `peptides.txtpeptidesSel.sumIntensity.csv` — sample of peptide intensity extraction for a selected protein.
- `peptides.txt.sampleBySampleLog2IntensityDiff.csv` — pairwise peptide log2 intensity differences (sample×sample columns).
- `peptides.txt.proteins.merged.sampleBySampleIntensityMedian.csv` — per-protein median of sample-by-sample intensity differences.

**Dependencies**
- **Python 3.8+** (tested in this repo environment)
- **Python packages:**: `pandas`, `numpy`, `scipy`

Install with pip if needed:

```bash
pip install pandas numpy scipy
```

**Notes about algorithm and behavior**
- **Missing data handling:**: zeros are replaced with `NaN` and medians are computed ignoring NaNs. Pairs with no overlapping peptides are skipped (no connection added between those samples).
- **Anchoring the solution:**: The MaxLFQ linear system has an arbitrary additive constant; the implementation anchors the system using the overall mean heuristic (`np.nanmean(X) * N`). This is a valid option but you may prefer a `sum(x)=0` constraint or a median-based anchor.
- **Disconnected components:**: If samples split into disconnected components (no shared peptides between groups), estimates are computed separately; the `annotation` result string indicates component labels.
- **Runtime:**: time is approximately O(N^2 * F) with N samples and F peptides per protein. For many samples or very large proteins this can be slow. Consider numba or precomputed masks if performance is an issue.

**How to verify / compare to MaxQuant outputs**
- If you also have a `proteinGroups.txt` (MaxQuant output), you can compare results: the repository includes a helper `iq/compare_maxlfq_to_proteingroups.py` which
  - matches protein IDs between our MaxLFQ output and `proteinGroups.txt`,
  - converts MaxQuant `LFQ intensity` to log2 space,
  - computes per-sample Pearson correlations and writes a merged CSV `maxlfq_vs_proteingroups_comparison.csv`.

Example run-and-compare:

```bash
# run MaxLFQ from peptides
python pepMaxLFQ.py peptides.txt

# compare with MaxQuant proteinGroups
python compare_maxlfq_to_proteingroups.py
```

**Suggested improvements / customizations**
- Add `argparse` to accept arbitrary input and output paths (currently the script defaults to `peptides.txt`).
- Expose anchoring method as an option (`sum=0`, median-based, or mean-based).
- Provide an option to skip creation of large intermediate files.
- Speedups: implement `maxLFQ_do` pairwise loop in `numba` or precompute peptide presence sets to avoid repeated mask creation.

**Troubleshooting**
- If you see `RuntimeWarning: All-NaN slice encountered` — this occurred previously for single-sample components with no quantified peptides; the implementation contains a guard that assigns `NA` in that case. Re-run after the guard is present.
- If your input uses different column names, update the script to point to the correct `Leading razor protein` or `Intensity` column patterns.

**Where to look in the repository**
- Script: `iq/pepMaxLFQ.py`
- Comparison helper: `iq/compare_maxlfq_to_proteingroups.py`
- ProteinGroups file (example): `iq/proteinGroups.txt`

If you want, I can now:
- Add `argparse` and make input/output paths configurable,
- Add example plots (scatter plots of our log2 vs MaxQuant log2) and save PNGs,
- Or add `numba` acceleration for large datasets.

---
Generated by the repo assistant — let me know which improvement to add next.
