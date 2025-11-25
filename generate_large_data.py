
import pandas as pd
import numpy as np
import random

def generate_large_data(filename, truth_filename, n_proteins=1000, n_samples=20, n_peptides_per_protein=5, n_fragments_per_peptide=3):
    print(f"Generating data with {n_proteins} proteins, {n_samples} samples...")
    
    proteins = [f"P{i}" for i in range(n_proteins)]
    samples = [f"S{i}" for i in range(n_samples)]
    
    rows = []
    truth_rows = []
    
    for p in proteins:
        # Generate peptides
        peptides = [f"{p}_Seq{j}" for j in range(n_peptides_per_protein)]
        
        # True abundances for this protein across samples
        # Let's make it interesting: some proteins change, some don't.
        # But for now, just random is fine.
        # We need to store the true abundance for each sample.
        
        for s in samples:
            # Base intensity for the protein in this sample
            protein_abundance = np.random.normal(20, 2) 
            
            truth_rows.append({
                'Protein': p,
                'Sample': s,
                'TrueAbundance': protein_abundance
            })
            
            for pep in peptides:
                # Generate fragments
                fragments = [f"y{k}" for k in range(1, n_fragments_per_peptide + 1)]
                
                for f in fragments:
                    # Peptide/Fragment ionization efficiency
                    # This should be constant for the same peptide/fragment across samples
                    # But here it is inside the sample loop?
                    # WAIT! frag_efficiency must be constant for the fragment across samples!
                    # The original code had it inside the sample loop?
                    # Let's check the original code.
                    pass

    # The original code:
    # for s in samples:
    #    protein_abundance = ...
    #    for f in fragments:
    #       frag_efficiency = np.random.normal(0, 1)
    
    # This implies frag_efficiency varies per sample! That violates the model.
    # The model assumes I_fs = A_s + E_f.
    # If E_f varies with s, it's just noise.
    # I must fix this in the generation logic too.
    
    rows = []
    truth_rows = []
    
    for p in proteins:
        peptides = [f"{p}_Seq{j}" for j in range(n_peptides_per_protein)]
        
        # Pre-calculate fragment efficiencies
        frag_efficiencies = {}
        for pep in peptides:
            fragments = [f"y{k}" for k in range(1, n_fragments_per_peptide + 1)]
            for f in fragments:
                frag_efficiencies[(pep, f)] = np.random.normal(0, 1)
        
        for s in samples:
            protein_abundance = np.random.normal(20, 2)
            
            truth_rows.append({
                'Protein': p,
                'Sample': s,
                'TrueAbundance': protein_abundance
            })
            
            for pep in peptides:
                fragments = [f"y{k}" for k in range(1, n_fragments_per_peptide + 1)]
                for f in fragments:
                    eff = frag_efficiencies[(pep, f)]
                    
                    log2_intensity = protein_abundance + eff + np.random.normal(0, 0.5)
                    
                    if random.random() < 0.1:
                        intensity = "NA"
                    else:
                        intensity = 2**log2_intensity
                        
                    rows.append({
                        'PG.ProteinGroups': p,
                        'R.Condition': s,
                        'EG.ModifiedSequence': pep,
                        'FG.Charge': 2,
                        'F.FrgIon': f,
                        'F.Charge': 1,
                        'F.PeakArea': intensity
                    })

    df = pd.DataFrame(rows)
    print(f"Generated {len(df)} rows.")
    df.to_csv(filename, sep="\t", index=False)
    print(f"Saved data to {filename}")
    
    truth_df = pd.DataFrame(truth_rows)
    # Pivot truth to match output format (Protein x Sample)
    truth_pivot = truth_df.pivot(index='Protein', columns='Sample', values='TrueAbundance')
    truth_pivot.to_csv(truth_filename, sep="\t")
    print(f"Saved truth to {truth_filename}")

if __name__ == "__main__":
    # 100 proteins * 5 peptides * 3 fragments * 20 samples = 30,000 rows (smaller for speed)
    generate_large_data("large_test_data.tsv", "ground_truth.tsv", n_proteins=100, n_samples=20, n_peptides_per_protein=5, n_fragments_per_peptide=3)
