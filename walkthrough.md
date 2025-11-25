
# Verification of iq.py Port

I have verified that the Python port (`iq.py`) produces results identical to the original R package (`iq.R`) within floating-point tolerance.

## Verification Steps

1.  **Run R Script**: Executed `run_r_test.R` to generate `r_result.tsv`.
2.  **Run Python Script**: Executed `run_py_test.py` to generate `py_result.tsv`.
3.  **Compare Results**: Created and ran `compare_outputs.py` to compare the two output files.
4.  **Unit Tests**: Ran `test_iq.py` to verify individual functions (`preprocess`, `maxLFQ`) against expected logic.

## Results

### Comparison
The maximum difference between the R and Python outputs for the test dataset is approximately **3.55e-14**, which is attributable to floating-point arithmetic differences between the two languages.

```
Comparing r_result.tsv and py_result.tsv...
R shape: (1, 3)
Py shape: (1, 3)
Max difference: 3.552713678800501e-14
SUCCESS: Results match within tolerance.
```

### Unit Tests
All tests in `test_iq.py` passed.

```
Testing maxLFQ...
maxLFQ Test PASSED

Testing preprocess...
preprocess Test PASSED

Testing full pipeline...
Pipeline Test PASSED
```

## Large Data Verification (Simulated Ground Truth)

I performed an additional verification using a larger simulated dataset where the ground truth protein abundances are known.

### Steps
1.  **Data Generation**: Modified `generate_large_data.py` to generate a dataset with 100 proteins, 20 samples, and known ground truth abundances. The data simulates peptide/fragment ionization efficiencies and biological variation.
2.  **Execution**: Ran both `run_r_large.R` and `run_py_large.py` on the generated `large_test_data.tsv`.
3.  **Verification**: Created `verify_results.py` to:
    *   Compare R and Python outputs.
    *   Compare Python outputs against the Ground Truth (using Pearson correlation and RMSD).

### Results

#### R vs Python
The results match extremely well, with differences negligible (likely due to floating point precision).

```
Comparing r_result_large.tsv and py_result_large.tsv...
Max difference (R vs Py): 1.31e-13
PASS: R and Python results match.
```

#### Python vs Ground Truth
The estimated abundances correlate highly with the ground truth, confirming the algorithm correctly recovers relative protein abundances from fragment intensities.

```
Comparing py_result_large.tsv and ground_truth.tsv...
Average Pearson Correlation: 0.9891
Average RMSD (centered): 0.2807
PASS: High correlation with ground truth.
```

## Conclusion
The `iq.py` script is a faithful port of `iq.R` and accurately recovers protein abundances from simulated data.
