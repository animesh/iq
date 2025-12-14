# python ./pepMaxLFQ.py peptides.txt
# wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/09/PXD000279/dynamicrangebenchmark.zip
# unzip dynamicrangebenchmark.zip
# %% setup
#install dependencies: pip install pandas scipy
import sys,numpy as np,pandas as pd
from scipy.sparse import csgraph, csr_matrix


def maxLFQ_do(X):
	# X is a numpy array (fragments x samples)
	N = X.shape[1]
    
	AtA = np.zeros((N, N))
	Atb = np.zeros(N)
    
	# Iterate over pairs of samples
	for i in range(N - 1):
		for j in range(i + 1, N):
			# Calculate ratio using only non-NA overlapping fragments
			mask = (~np.isnan(X[:, i])) & (~np.isnan(X[:, j]))
			if mask.sum() == 0:
				continue
			diff = -X[mask, i] + X[mask, j]
			r_i_j = np.nanmedian(diff)
            
			if not np.isnan(r_i_j):
				AtA[i, j] = -1
				AtA[j, i] = -1
				AtA[i, i] += 1
				AtA[j, j] += 1
                
				Atb[i] -= r_i_j
				Atb[j] += r_i_j
                
	# Solve system
	A_final = np.vstack([
		np.hstack([2 * AtA, np.ones((N, 1))]),
		np.hstack([np.ones(N), [0]])
	])
    
	# Anchor using overall mean (heuristic)
	b_final = np.concatenate([2 * Atb, [np.nanmean(X) * N]])
    
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
	mask = ~np.isnan(X_np)
	adj = (mask.T @ mask) > 0
	n_components, labels = csgraph.connected_components(csr_matrix(adj), directed=False)
    
	w = np.full(N, np.nan)
    
	for i in range(n_components):
		ind = np.where(labels == i)[0]
        
		if len(ind) == 1:
			# handle possible all-NaN column explicitly to avoid RuntimeWarning
			col = X_np[:, ind[0]]
			if np.all(np.isnan(col)):
				w[ind] = np.nan
			else:
				w[ind] = np.nanmedian(col)
		else:
			w[ind] = maxLFQ_do(X_np[:, ind])
            
	if np.all(np.isnan(w)):
		return {'estimate': w, 'annotation': "NA"}
	else:
		quantified_samples = ~np.isnan(w)
		q_labels = labels[quantified_samples]
        
		if len(q_labels) > 0 and np.all(q_labels == q_labels[0]):
			return {'estimate': w, 'annotation': ""}
		else:
			anno = labels.astype(object)
			anno[np.isnan(w)] = np.nan
			return {'estimate': w, 'annotation': ";".join(str(x) for x in anno)}
# %% data
pathFiles = sys.argv[1]
#pathFiles = "peptides.txt"
dfS=pd.read_csv(pathFiles,low_memory=False,sep='\t')
print(dfS.columns)
print(dfS.head())
#plt.close('all')
#import matplotlib
#matplotlib.use('Agg')
#plt.plot(df['Intensity'])
#plt.show()
#dfS=df[df['PEP']<0.01]
#print(dfS['Mod. peptide IDs'].value_counts())
#plt.plot(dfS['Intensity'])
#dfS.rename({'Leading razor protein':'ID'},inplace=True,axis='columns')
# %% separate IDs by semicolon
dfS.rename({'Leading razor protein':'ID'},inplace=True,axis='columns')
print(dfS[dfS==0].count())
#dfS.replace(0, np.nan, inplace=True)
#print(dfS[dfS==0].count())
dfS['IDs']=dfS.ID.str.split(';')
dfSE=dfS.explode('IDs')
dfSE.index=dfSE["Sequence"]
# %% groupy IDs using sum of intensities
dfSEG=dfSE.groupby(dfSE['IDs']).aggregate('sum')
dfSEG.to_csv(pathFiles+'.combinedIntensity.csv')#,sep="\")#,rownames=FALSE)
# %% check for a protein
dfA5A614=dfSE[dfSE['IDs'].isin(["A5A614"])]
dfA5A614.index=dfA5A614.Sequence
dfA5A614=dfA5A614.filter(regex='Intensity ',axis=1)
dfA5A614.columns=dfA5A614.columns.str.removeprefix('Intensity ')
#dfA5A614.columns=dfA5A614.columns.str.split('_').str[0]
#dfA5PFJ8=dfSE[dfSE['IDs'].isin(["A5PFJ8"])]
#dfSE[dfSE['IDs'].isin(["A5PFJ8"])]
#dfA5A614.T.plot.line().figure.savefig(pathFiles+'peptidesSel.sumIntensity.png',dpi=100,bbox_inches = "tight")
dfA5A614.to_csv(pathFiles+'peptidesSel.sumIntensity.csv')#,sep="\")#,rownames=FALSE)
# %% extract intensity
dfSE.replace(0,np.nan,inplace=True)
dfSEint=dfSE.filter(regex='Intensity ',axis=1)
dfSEint.columns=dfSEint.columns.str.removeprefix('Intensity ')
print(dfSEint.describe())
# %% tranform with log2
dfSEintLog2=np.log2(dfSEint)
print(dfSEintLog2.describe())
# %% compute MaxLFQ per protein
print("Computing MaxLFQ per protein...")

if 'IDs' not in dfSE.columns:
	raise ValueError("Expected column 'IDs' in exploded peptide table")

protein_ids = pd.Index(dfSE['IDs'].unique())
samples = dfSEintLog2.columns

# Prepare output table
prot_tab = pd.DataFrame(index=protein_ids, columns=samples, dtype=float)
prot_anno = {}

for pid in protein_ids:
	mask = dfSE['IDs'] == pid
	sub = dfSEintLog2[mask]
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
		# If est is a scalar or differently shaped, try to broadcast
		prot_tab.loc[pid] = np.asarray(est)

	prot_anno[pid] = anno

out_fname = pathFiles + '.proteinMaxLFQ.tsv'
prot_tab.to_csv(out_fname, sep='\t', na_rep='NA')
print(f"Wrote MaxLFQ results to {out_fname}")
# %% pairwise substract log2 sample values per peptide
dfSEintRepS=dfSEintLog2[np.repeat(dfSEintLog2.columns.values,dfSEintLog2.shape[1])]
#dfSEintRepS.reset_index(inplace=True)
dfSEintRepS.columns
dfSEintRepB=pd.concat([dfSEintLog2]*(dfSEintLog2.shape[1]),axis=1)
#dfSEintRepB.reset_index(inplace=True)
dfSEintRepB.columns
dfSEintRepS.columns=dfSEintRepS.columns+';'+dfSEintRepB.columns
dfSEintRepB.columns=dfSEintRepS.columns
dfSEintRepSB=dfSEintRepS-dfSEintRepB
dfSEintRepSB.replace(0,np.nan,inplace=True)
dfSEintRepSB.dropna(axis = 1, how = 'all',inplace=True)
print(dfSEintRepSB.describe())
dfSEintRepSB["Protein"]=dfSE["IDs"]
dfSEintRepSB=dfSEintRepSB.sort_values('Protein')
dfSEintRepSB.to_csv(pathFiles+'sampleBySampleLog2IntensityDiff.csv')
print(dfSEintRepSB.columns)
print("Intensity Difference Summary\n",dfSEintRepSB.describe())
# %% group peptides to protein using median of log2differences 
dfSEintRepSBPG=dfSEintRepSB.groupby(dfSEintRepSB['Protein']).aggregate('median')
#print(dfSEintRepSBPG.describe())
#dfSEintRepSBPG.hist()
#dfSEintRepSBPG.replace(np.nan,0,inplace=True)
dfSEintRepSBPG.to_csv(pathFiles+'proteins.merged.sampleBySampleIntensityMedian.csv')
#print(dfSEintRepSBPG.columns)
print("Median Intensity Difference Summary\n",dfSEintRepSBPG.describe())
