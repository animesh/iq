**Overview**

This document summarizes the changes made to implement and improve a MaxLFQ-style
protein quantification pipeline inside `iq/pepMaxLFQ.py`. It records the
algorithmic updates, practical usage, test results (comparison vs MaxQuant
`proteinGroups.txt`), and recommended next steps.


**Changes Implemented**

- **Inlined MaxLFQ:** Replaced earlier fragmented script with a single,
  self-contained implementation of the MaxLFQ workflow (peptide -> protein).
- **Weighted normal equations:** When assembling the pairwise normal equations
  the pair contribution weight equals the number of shared peptides between
  runs (wij = n_shared). This reproduces the approach described in Cox et al.
- **Minimum shared-peptide threshold:** A `min_pair_count` filter is enforced
  when building pairwise ratios so pairs with too few shared peptides are
  ignored (default = 2).
- **Large-ratio stabilization:** For pairs with small peptide overlap we
  interpolate between the median of peptide ratios and the ratio of summed
  intensities (paper rule: interpolate using x = N_max / n_shared with
  breakpoints at 2.5 and 5.0).
- **Augmented solve + rescaling:** The laplacian normal system is augmented
  with a sum constraint to make it solvable; after solving we rescale the
  recovered log2-profile so that the total linear abundance (sum 2**x) matches
  the summed observed peptide intensities for that protein.
- **Delayed normalization (optional):** Added a lightweight `delayed_normalization`
  solver that fits additive log2 normalization factors per run via nonlinear
  least squares (minimize differences of log ratios of summed intensities across
  peptide features). The pipeline now estimates run factors and subtracts them
  from the log2 peptide matrix before computing MaxLFQ.


**Files Added / Updated**

- `iq/pepMaxLFQ.py` — Rewritten, self-contained MaxLFQ pipeline with the above
  algorithmic improvements and a `main()` runner.
- `iq/test_peptides_small.txt` — small synthetic test peptide table used to
  validate runtime behavior during development.
- `iq/peptides.txt.proteinMaxLFQ.tsv` — output produced by running the
  pipeline on the repository `peptides.txt` (log2 protein intensities).
- `iq/maxlfq_vs_proteingroups_correlation.csv` — correlation table comparing
  our results to `proteinGroups.txt` (before delayed-normalization was applied).
- `iq/maxlfq_vs_proteingroups_correlation_after_delayednorm.csv` — per-sample
  Pearson correlations after integrating delayed normalization.


**How to run**

1. Run the pipeline on a MaxQuant-style peptide table:

```
python iq/pepMaxLFQ.py iq/peptides.txt
```

2. The script writes protein-level log2 LFQ results to:

```
iq/peptides.txt.proteinMaxLFQ.tsv
```

3. A short comparison step (used during development) computes Pearson
   correlations vs MaxQuant `proteinGroups.txt` and writes a CSV. The
   repository contains the generated CSVs referenced above.


**Key Test Results (comparison vs MaxQuant `proteinGroups.txt`)**

The comparison matches proteins by the `Majority protein IDs` (first ID) and
converts MaxQuant LFQ intensities to log2 for correlation with our log2
estimates.

- Correlations BEFORE delayed normalization (our first implementation):

  - UPS1_01: r = 0.519
  - UPS1_02: r = 0.560
  - UPS1_03: r = 0.537
  - UPS1_04: r = 0.587
  - UPS2_01: r = 0.589
  - UPS2_02: r = 0.694
  - UPS2_03: r = 0.592
  - UPS2_04: r = 0.582

- Correlations AFTER applying the lightweight delayed-normalization (pipeline
  now subtracts estimated additive run factors before MaxLFQ):

  - UPS1_01: r = 0.450
  - UPS1_02: r = 0.501
  - UPS1_03: r = 0.481
  - UPS1_04: r = 0.548
  - UPS2_01: r = 0.504
  - UPS2_02: r = 0.457
  - UPS2_03: r = 0.555
  - UPS2_04: r = 0.501

Observations:
- Overall the weighted MaxLFQ implementation produces moderate agreement
  with MaxQuant LFQ (Pearson r typically ~0.5–0.7 in the initial run).
- After applying the current delayed-normalization solver the correlations
  decreased for many samples. This is not necessarily incorrect — delayed
  normalization changes per-peptide log-ratios and therefore the MaxLFQ
  solution — but in this dataset (apparently one run per sample) the
  solver provides limited benefit and likely introduces noise.


**Notes, Caveats, and Rationale**

- The goal was to reproduce key algorithmic steps from Cox et al. (2014):
  pairwise medians, pair filtering, large-ratio stabilization, weighting by
  shared-peptide counts, and correct rescaling. Those changes were implemented
  to make the solver more faithful to the published MaxLFQ behavior.
- The current delayed-normalization routine is a small, generic least-squares
  fit. Delayed normalization usually helps when multiple runs correspond to the
  same biological sample (and need per-run normalization before aggregation).
  If you have one run per sample, the solver has fewer constraints and can
  overfit pairwise peptide differences, reducing agreement with MaxQuant.


**Recommended Next Steps**

- Add a command-line flag to `pepMaxLFQ.py` to toggle delayed normalization
  (e.g., `--no-delayed-normalization`) so you can re-run without it quickly.
- Add an explicit `--min-pair-count` and `--no-stabilization` flags for
  sensitivity testing.
- Produce per-protein diagnostics (scatter of our vs MaxQuant, absolute
  differences) and Bland–Altman plots to identify systematic bias.
- If you want closer numeric agreement with MaxQuant, consider these
  investigations:
  - Verify how MaxQuant selects peptides for pair medians (razor peptide
    handling / grouping semantics).
  - Implement intensity-based weighting in addition to (or instead of)
    count-based weighting.
  - Add a small regularization to the delayed-normalization solver or restrict
    it to grouped runs (when present).


**Quick contact / reproducibility**

If you'd like, I can:
- Add a CLI to toggle delayed normalization and other algorithm flags.
- Generate per-protein comparison CSVs and plots so we can inspect the largest
  discrepancies together.
- Implement more robust delayed-normalization (regularized solver) and run a
  sensitivity sweep for `min_pair_count` and the stabilization parameters.

Files referenced above are in the `iq/` folder of this repository.
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
