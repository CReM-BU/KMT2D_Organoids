```python
# Below is the MultiVelo to assess temporal dynamics of our cerebral organoid multiome data. 
# Multivelo is a computational tool that models gene regulation by integrating scRNA-seq and snATAC-seq data (Li et al., 2023. PMCID: PMC10246490). 
# Our goal is to investigate the dynamic relationship between gene expression and chromatin accessibility in wild-type (WT) and KMT2D-LOF samples.

# Import required packages
import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

# scanpy, scvelo, and multivelo are used for processing, visualizing, and modeling single-cell RNA and ATAC data
# matplotlib.pyplot is used for plotting
# os, numpy, pandas, scipy are general Python libraries for data handling and scientific computing
```


```python
# Set environment settings
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)
```


```python
# Load and clean RNA data

# Read RNA velocity-formatted .loom file containing spliced and unspliced reads 
adata_rna = scv.read("WT_aggregated_corrected.loom", cache=True)

# Reformat cell barcodes to match expected naming (e.g., for integration with ATAC)
adata_rna.obs_names = [x.split(':')[1][:-1] + '-1' for x in adata_rna.obs_names]

# Ensures gene names are unique and computes total counts per cell
adata_rna.var_names_make_unique()
adata_rna.obs['n_counts'] = adata_rna.X.sum(axis=1).A1

adata_rna
```




    AnnData object with n_obs × n_vars = 12142 × 58395
        obs: 'n_counts'
        var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'
        layers: 'ambiguous', 'matrix', 'spliced', 'unspliced'




```python
# Plot total counts per cell to identify low- and high-count outliers
sc.pl.violin(adata_rna, ['n_counts'], log=True)
```


    
![png](output_3_0.png)
    



```python
# Remove cells with too few or too many transcripts (possibly dead or doublets)
sc.pp.filter_cells(adata_rna, min_counts=800)
sc.pp.filter_cells(adata_rna, max_counts=30000)

adata_rna
```




    AnnData object with n_obs × n_vars = 6907 × 58395
        obs: 'n_counts'
        var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'
        layers: 'ambiguous', 'matrix', 'spliced', 'unspliced'




```python
sc.pl.violin(adata_rna, ['n_counts'], log=True)
```


    
![png](output_5_0.png)
    



```python
# Filter genes based on how often they’re shared and normalize expression
scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10)
```

    Filtered out 44184 genes that are detected 10 counts (shared).
    Normalized count data: X, spliced, unspliced.
    Logarithmized X.


    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/preprocessing/utils.py:705: DeprecationWarning: `log1p` is deprecated since scVelo v0.3.0 and will be removed in a future version. Please use `log1p` from `scanpy.pp` instead.
      log1p(adata)



```python
# Ensure unique adata_rna barcodes

# Check and remove duplicate barcodes
is_unique = adata_rna.obs_names.is_unique
print(f"Are barcodes unique? {is_unique}")

duplicates = adata_rna.obs_names[adata_rna.obs_names.duplicated()]
print(f"Duplicate barcodes: {duplicates}")
adata_rna = adata_rna[~adata_rna.obs_names.duplicated()].copy()

is_unique = adata_rna.obs_names.is_unique
print(f"Are barcodes unique? {is_unique}")

adata_rna
```

    Are barcodes unique? False
    Duplicate barcodes: Index(['ATATGCATCCCGCAAA-1', 'CGGCTAATCAATCTAG-1', 'GAGCCACTCAAGGACA-1',
           'TAAGCTATCACTTCAT-1', 'TGCTTGTGTCGCAAAC-1', 'TTATTGCTCTAAGTCA-1',
           'AGTCGCATCTGTGAGT-1', 'CCCGCTTCAACTCGCG-1', 'CGGTTTCTCGGGATTT-1',
           'CTACAACAGGAACCAA-1', 'CTGGCTTTCCCGCATT-1', 'GCCTGCTGTGCAATGC-1',
           'GTAGGATCAACTCGCG-1', 'TGGACCGGTTTATTCG-1', 'AAGGTGCAGTAGGATG-1',
           'CTTTGTCCATACCCGG-1', 'GATTAAGCAACACTTG-1', 'GCTCGATCAATACTGT-1',
           'GGTACTAGTTCCGGGA-1', 'GTATGTGGTAATCCCT-1', 'TAAGTAGCATTAGGTT-1',
           'TAGCGCGGTGATGAAA-1'],
          dtype='object')
    Are barcodes unique? True





    AnnData object with n_obs × n_vars = 6885 × 14211
        obs: 'n_counts', 'initial_size_unspliced', 'initial_size_spliced', 'initial_size'
        var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand', 'gene_count_corr'
        uns: 'log1p'
        layers: 'ambiguous', 'matrix', 'spliced', 'unspliced'




```python
# Match barcodes with external cell type annotations

# First used shell commands to re-format the barcode names
#cut -d '_' -f1,3 cell_annotations.tsv | sed 's/_1//g' | sed 's/_2//g' > cell_annotations_2.tsv
#then extra barcodes from adata_rna

print(adata_rna.obs_names[:10])

# Export barcodes for grep-matching with the external cell annotation file
barcodes_df = pd.DataFrame(adata_rna.obs_names, columns=['Barcode'])
barcodes_df.to_csv("barcodes.csv", index=False)

# the extarcted barcodes are then greped out from the cell_annotaion_2.tsv (run in terminal)
# grep -wFf barcodes.csv cell_annotations_2.tsv | awk '!seen[$1]++' > cell_annotations_updated.tsv
```

    Index(['AACAAGCCACCGGCTA-1', 'AAGACATAGCCATCAG-1', 'AATCAGGAGAGGATAT-1',
           'AACTACTCAATGCGCT-1', 'AAGGCCCTCCTTGCGT-1', 'AAAGCTTGTGTTCCCA-1',
           'AACTAGCTCATTACGA-1', 'AATCTCAAGGCGCTTA-1', 'AAAGCGGGTGCAATGC-1',
           'AAAGCTTGTTAATCGG-1'],
          dtype='object')



```python
# Load matched cell-type annotations and merge them into adata_rna
cell_annot = pd.read_csv('cell_annotations_updated.tsv', sep='\t', index_col=0)
adata_rna = adata_rna[cell_annot.index,:]
adata_rna.obs['celltype'] = cell_annot['celltype']
adata_rna

#above is the finl RNA object and we have 6882 cells with 14211 genes
```


```python
# Save the cleaned and annotated RNA data
adata_rna.write('adata_rna.h5ad')
```


```python
# Load and process ATAC data

# Load 10x Genomics ATAC matrix, select only peak features
adata_atac = sc.read_10x_mtx('/restricted/projectnb/crem-bioinfo/project_workspace/22_05_13_angie/MultiVelo/cellranger/WT/WTaggr/outs/filtered_feature_bc_matrix', var_names='gene_symbols', cache=True, gex_only=False)
adata_atac = adata_atac[:,adata_atac.var['feature_types'] == "Peaks"]
```


```python
# We aggregate peaks around each gene as well as those that have high correlations with promoter peak or gene expression.
# Peak annotation contains the metadata for all peaks.
# Feature linkage contains pairs of correlated genomic features.

adata_atac = mv.aggregate_peaks_10x(adata_atac,'atac_peak_annotation.tsv','feature_linkage.bedpe')
adata_atac
```

    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/ipywidgets/widgets/widget.py:503: DeprecationWarning: The `ipykernel.comm.Comm` class has been deprecated. Please use the `comm` module instead.For creating comms, use the function `from comm import create_comm`.
      self.comm = Comm(**args)



      0%|          | 0/21086 [00:00<?, ?it/s]





    View of AnnData object with n_obs × n_vars = 12142 × 21086




```python
# Filter and normalize ATAC
# Examine the total count distribution and remove outliers
plt.hist(adata_atac.X.sum(1), bins=100, range=(0, 100000));
```


    
![png](output_13_0.png)
    



```python
# Filter cells based on peak count distribution
sc.pp.filter_cells(adata_atac, min_counts=2000)
sc.pp.filter_cells(adata_atac, max_counts=60000)
plt.hist(adata_atac.X.sum(1), bins=100, range=(0, 100000));
```


    
![png](output_14_0.png)
    



```python
# Normalize aggregated peaks with TF-IDF, which adjusts for sequencing depth and emphasizes informative peaks
mv.tfidf_norm(adata_atac)
```


```python
print(adata_atac.obs_names[:10])
```

    Index(['AAACAGCCAATACTGT-4', 'AAACAGCCAATCTCTC-3', 'AAACAGCCAATTAACC-3',
           'AAACAGCCACAGAACG-4', 'AAACAGCCAGAATGAC-2', 'AAACAGCCAGGAACCA-1',
           'AAACAGCCATTAGCCA-4', 'AAACAGCCATTCCTGT-4', 'AAACAGCCATTGACAT-2',
           'AAACATGCAATCATGT-3'],
          dtype='object')



```python
# Match barcode naming convention to RNA data (-1 suffix)
# Since barcode of atac has -2,-3 and-4 in the end, the overlap with rna is low

# Update obs_names to ensure all barcodes end with "-1"
adata_atac.obs_names = [barcode.split('-')[0] + '-1' for barcode in adata_atac.obs_names]
print(adata_atac.obs_names[:10])
```

    Index(['AAACAGCCAATACTGT-1', 'AAACAGCCAATCTCTC-1', 'AAACAGCCAATTAACC-1',
           'AAACAGCCACAGAACG-1', 'AAACAGCCAGAATGAC-1', 'AAACAGCCAGGAACCA-1',
           'AAACAGCCATTAGCCA-1', 'AAACAGCCATTCCTGT-1', 'AAACAGCCATTGACAT-1',
           'AAACATGCAATCATGT-1'],
          dtype='object')



```python
#first make adata_atac unique
is_unique = adata_atac.obs_names.is_unique
print(f"Are barcodes unique? {is_unique}")

duplicates = adata_atac.obs_names[adata_atac.obs_names.duplicated()]
print(f"Duplicate barcodes: {duplicates}")

adata_atac = adata_atac[~adata_atac.obs_names.duplicated()].copy()

is_unique = adata_atac.obs_names.is_unique
print(f"Are barcodes unique? {is_unique}")

adata_atac
```

    Are barcodes unique? False
    Duplicate barcodes: Index(['AACTAGCTCAATGAGG-1', 'AAGCATGAGGCCAATT-1', 'AAGCCTGTCCTAGTAA-1',
           'AAGGATCCAACTGGCT-1', 'AAGGTGCAGTAGGATG-1', 'AATTAGCGTCACAGCG-1',
           'ACTTATCTCACCGGTA-1', 'AGGTCAAAGTGACCTG-1', 'AGTCGCATCTGTGAGT-1',
           'ATATGCATCCCGCAAA-1', 'ATCACCCTCGCTATAA-1', 'ATCGGCCAGCTATTGA-1',
           'ATCTATGAGTCTTGAA-1', 'ATGGCTAGTCGTAAAT-1', 'ATTGCACAGTAACCCG-1',
           'CAATGAACAGGCTGTT-1', 'CAGGTCCAGAGCCGCT-1', 'CCAGACTCATCTAGCA-1',
           'CCATAATCACCATATG-1', 'CCATCACTCTTGACCC-1', 'CCCAATTGTGTTTGTC-1',
           'CGATTTGCAATCATGT-1', 'CGGTTTCTCGGGATTT-1', 'CTACAACAGGAACCAA-1',
           'CTAGTAATCAAGACTC-1', 'CTCTATGTCCTCACAC-1', 'CTGGTTTGTAATGGCC-1',
           'CTTTGTCCATACCCGG-1', 'GACCGAACAAGGTAAC-1', 'GACCTAAGTAGACAAA-1',
           'GAGAGGCGTAAAGCGG-1', 'GAGCCACTCAAGGACA-1', 'GAGCCTTCACGTGCTG-1',
           'GATGGCTGTGATTTGG-1', 'GATTAAGCAACACTTG-1', 'GATTATGTCCTTCTAG-1',
           'GCAAGTGCAGGCCAAA-1', 'GCTCGATCAATACTGT-1', 'GGAAGTATCCGTAAAC-1',
           'GGTACTAGTTCCGGGA-1', 'GTAGGATCAACTCGCG-1', 'GTTCGCGCAGATAGAC-1',
           'GTTCTTGTCCTTCTAG-1', 'GTTTACCGTTAGTTGG-1', 'TAAGCTATCACTTCAT-1',
           'TAAGTAGCATTAGGTT-1', 'TATAGCCAGGATCACT-1', 'TCATACTTCCAAGTGT-1',
           'TCCTTGCAGTAAAGGT-1', 'TGCTTAAAGCTATTAG-1', 'TGCTTGTGTAGCTGGT-1',
           'TGCTTGTGTCGCAAAC-1', 'TGGACCGGTTTATTCG-1', 'TGGTCCTTCGTAATCA-1',
           'TTATTGCTCTAAGTCA-1', 'TTGCAATCATAGTCAT-1'],
          dtype='object')
    Are barcodes unique? True





    AnnData object with n_obs × n_vars = 9841 × 21086
        obs: 'n_counts'




```python
#Finding shared barcodes and features between RNA and ATAC

# Keep only cells and genes found in both RNA and ATAC datasets
shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
len(shared_cells), len(shared_genes)
```




    (5602, 12283)




```python
# Keep only cells and genes found in both RNA and ATAC datasets
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]
```


```python
#after selecting only shared cell and shared genes, we are left with 5602 cells and 12283 genes with peaks
adata_rna
adata_atac
```




    View of AnnData object with n_obs × n_vars = 5602 × 12283
        obs: 'n_counts'




```python
# Dimensionality reduction (for UMAP + smoothing)
scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
sc.pp.neighbors(adata_rna, n_pcs=30, n_neighbors=50)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)
```

    WARNING: Did not normalize X as it looks processed already. To enforce normalization, set `enforce=True`.
    WARNING: Did not normalize spliced as it looks processed already. To enforce normalization, set `enforce=True`.
    WARNING: Did not normalize unspliced as it looks processed already. To enforce normalization, set `enforce=True`.
    WARNING: adata.X seems to be already log-transformed.


    /scratch/2038529.1.crem/ipykernel_1189482/1384757388.py:3: DeprecationWarning: `log1p` is deprecated since scVelo v0.3.0 and will be removed in a future version. Please use `log1p` from `scanpy.pp` instead.
      scv.pp.log1p(adata_rna)


    computing moments based on connectivities
        finished (0:00:08) --> added 
        'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)



```python
scv.tl.umap(adata_rna)
scv.pl.umap(adata_rna, color='celltype')
```

    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])
    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])
    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])
    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])
    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])
    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])
    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])
    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/scvelo/plotting/utils.py:63: DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, pd.CategoricalDtype) instead
      return isinstance(c, str) and c in data.obs.keys() and cat(data.obs[c])



    
![png](output_23_1.png)
    



```python
# Neighbor smoothing with Seurat WNN results

# Write out filtered cells and prepare to run Seurat WNN --> Refer to the R script on Github

# Export final filtered cells for use in Seurat's weighted nearest neighbors (WNN) in R
adata_rna.obs_names.to_frame().to_csv('filtered_cells.txt', header=False, index=False)
```


```python
# The data above are the final cells. We have generted Seurat WNN based on the script shared. For it we have used the aggr WT object. See the RMD file in the folder.

# Load neighbor index and distance from Seurat
nn_idx = np.loadtxt("nn_idx_update.txt", delimiter=',')
nn_dist = np.loadtxt("nn_dist_update.txt", delimiter=',')
nn_cells = pd.Index(pd.read_csv("nn_cells_update.txt", header=None)[0])

adata_atac = adata_atac[nn_cells]

# Make sure cell names match
np.all(nn_cells == adata_atac.obs_names)
```




    True




```python
# Apply KNN smoothing to ATAC data to reduce sparsity
mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)
```


```python
adata_atac
```




    AnnData object with n_obs × n_vars = 5602 × 12283
        obs: 'n_counts'
        layers: 'Mc'
        obsp: 'connectivities'




```python
#save both RNA and ATAC data individually
adata_atac.write('adata_atac.h5ad')
adata_rna.write('adata_rna.h5ad')
```


```python
# Run MultiVelo: model transcriptional & chromatin kinetics
# This will take a while. Parallelization is high recommended.
# Uses "invert" initialization mode (based on a backwards model guess)
# n_anchors=500 controls how many representative cells are used for speed

adata_result = mv.recover_dynamics_chrom(adata_rna,
                                         adata_atac,
                                         max_iter=5,
                                         init_mode="invert",
                                         parallel=True,
                                         save_plot=False,
                                         rna_only=False,
                                         fit=True,
                                         n_anchors=500,
                                         extra_color_key='celltype'
                                        )
```

    /projectnb/crem-bioinfo/bpushpin/.conda/envs/velotest/lib/python3.10/site-packages/ipywidgets/widgets/widget.py:503: DeprecationWarning: The `ipykernel.comm.Comm` class has been deprecated. Please use the `comm` module instead.For creating comms, use the function `from comm import create_comm`.
      self.comm = Comm(**args)



      0%|          | 0/12283 [00:00<?, ?it/s]



```python
# Save the result for use later on
adata_result.write("multivelo_result.h5ad")
```


```python
# Load final results
adata_result = sc.read_h5ad("multivelo_newWT.h5ad")
adata_result2 = sc.read_h5ad("multivelo_newKO.h5ad")

```


```python
# Plot summary of model genes in pie chart

#WT
mv.pie_summary(adata_result)

#KMT2D-LOF
mv.pie_summary(adata_result2)
```


    
![png](output_32_0.png)
    



    
![png](output_32_1.png)
    



```python
#Components in the Pie Chart
#Induction:Represents genes that show increased transcription activity, meaning that their unspliced RNA fraction is higher relative to spliced RNA.
#These genes are being actively transcribed.
#Repression:Represents genes undergoing decreased transcription activity, meaning their spliced RNA fraction is higher relative to unspliced RNA.
#These genes are slowing down their expression.
#Model 1 & Model 2:These refer to different mathematical models used to estimate RNA velocity based on splicing kinetics.
#They help categorize genes based on their transcriptional state transitions.
#The exact difference between Model 1 and Model 2 depends on the scVelo settings but generally reflects different assumptions about the balance between splicing and degradation.

#If Induction is dominant → many genes are being actively transcribed.
#If Repression is dominant → many genes are shutting down transcription.
#If Model 1 and Model 2 are well-represented → different groups of genes follow distinct velocity patterns
```


```python
# Summarize the distributions of switch times for each model gene class

#WT 
mv.switch_time_summary(adata_result)

#KMT2D-LOF
mv.switch_time_summary(adata_result2)
```


    
![png](output_34_0.png)
    



    
![png](output_34_1.png)
    



```python
#Prime:Represents genes that are in a "primed" state for activation. These genes show a strong pre-transcriptional state, often marked by an accumulation of unspliced RNA.
#It suggests that these genes are prepared for transcription but haven't fully started expressing mature RNA yet.
#Coupled-on:These are genes that are actively transcribing and show a coupling between spliced and unspliced RNA.
#The unspliced and spliced RNA fractions are both high, indicating active transcription and processing. The term "on" suggests that these genes are engaged in transcriptional activity.
#Decoupled:Genes in this state show disconnection between splicing and transcription. There might be a mismatch between the unspliced and spliced RNA dynamics.
#For example, genes that are actively transcribed but where splicing doesn't proceed as efficiently, leading to an accumulation of unspliced RNA.
#Coupled-off:These genes are repressed or off in terms of transcription, where splicing has been completed, and unspliced RNA is minimal.
#Essentially, these genes have turned off transcription, and the mature (spliced) RNA is the dominant form.

#Each box represents the distribution of switch times for genes that fall into one of these categories.
#Prime: Likely to have early switch times, as genes are prepped for activation.
#Coupled-on: Might show mid-phase switch times, reflecting ongoing transcription.
#Decoupled: Could show diverse switch times, indicating a mixture of active transcription and delayed splicing.
#Coupled-off: Genes in this state would have later switch times, indicating that transcription has been turned off and splicing is complete.
```


```python
# plot likelihood and variable distributions

#WT
mv.likelihood_plot(adata_result)

#KMT2D-LOF
mv.likelihood_plot(adata_result2)
```


    
![png](output_36_0.png)
    



    
![png](output_36_1.png)
    



```python
# Computing velocity stream and latent time

# Calculates the velocity graph by evaluating cosine similarity between velocity vectors and potential future cell states. 
# This measures how closely the observed direction of gene expression change aligns with the predicted trajectory indicated by the RNA velocity.

#WT
mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

#KMT2D-LOF
mv.velocity_graph(adata_result2)
mv.latent_time(adata_result2)

scv.pl.scatter(adata_result, color='latent_time', color_map='gnuplot', size=80, save='WT_latenttime.png', dpi=300)
scv.pl.scatter(adata_result2, color='latent_time', color_map='gnuplot', size=80, save='KO_latenttime.png', dpi=300)
```

    computing velocity graph (using 1/20 cores)



      0%|          | 0/5602 [00:00<?, ?cells/s]


        finished (0:04:44) --> added 
        'velo_s_norm_graph', sparse matrix with cosine correlations (adata.uns)
    computing terminal states
    WARNING: Uncertain or fuzzy root cell identification. Please verify.
        identified 1 region of root cells and 1 region of end points .
        finished (0:00:00) --> added
        'root_cells', root cells of Markov diffusion process (adata.obs)
        'end_points', end points of Markov diffusion process (adata.obs)
    computing latent time using root_cells as prior
        finished (0:00:03) --> added 
        'latent_time', shared time (adata.obs)



```python
# Embed stream plots onto UMAPs

#WT
mv.velocity_embedding_stream(adata_result, basis='umap', color='celltype', legend_loc='right', save='velocity_WT.png', dpi=300)
mv.velocity_embedding_stream(adata_result, basis='umap', color='anno_7', legend_loc='right', save='velocity_WT_anno7.png', dpi=300)

#KMT2D-LOF
mv.velocity_embedding_stream(adata_result2, basis='umap', color='celltype', legend_loc='right', save='velocity_KO.png', dpi=300)
mv.velocity_embedding_stream(adata_result2, basis='umap', color='anno_7', legend_loc='right', save='velocity_KO_anno7.png', dpi=300)
```


    
![png](output_38_0.png)
    



    
![png](output_38_1.png)
    



    
![png](output_38_2.png)
    



    
![png](output_38_3.png)
    



```python
# Embed velocities onto UMAPs as grids

#WT
scv.pl.velocity_embedding_grid(
    adata_result,
    basis='umap',
    color='celltype', 
    vkey='velo_s_norm',  # matches the graph key in .uns
    show=True,
    arrow_size='25', save='WT_Grid.png', dpi=300)

scv.pl.velocity_embedding_grid(
    adata_result,
    basis='umap',
    color='anno_7', 
    vkey='velo_s_norm',  # matches the graph key in .uns
    show=True,
    arrow_size='25', save='WT_Grid_Anno7.png', dpi=300)

#KMT2D-LOF
scv.pl.velocity_embedding_grid(
    adata_result2,
    basis='umap',
    color='celltype', 
    vkey='velo_s_norm',  # matches the graph key in .uns
    show=True,
    arrow_size='25', save='LOF_Grid.png', dpi=300)

scv.pl.velocity_embedding_grid(
    adata_result2,
    basis='umap',
    color='anno_7',
    vkey='velo_s_norm',  # matches the graph key in .uns
    show=True,
    arrow_size='25' , save='LOF_Grid_Anno7.png', dpi=300)
```


    
![png](output_39_0.png)
    



    
![png](output_39_1.png)
    



    
![png](output_39_2.png)
    



    
![png](output_39_3.png)
    



```python
# Extract the genes driving the velocity vectors
# This is done using a differential expression analysis (Welch’s t-test with conservatively inflated variance) on the velocity expression matrix to detect genes that exhibit distinct transcriptional dynamics in one cluster relative to others, such as being actively induced in one group but stable elsewhere.
# If predefined clusters are not provided, velocity-based clustering is first performed using Louvain modularity on the velocity expression data.

#WT
scv.tl.rank_velocity_genes(adata_result, groupby='celltype', min_corr=.3, vkey='velo_s_norm', n_genes=10000)
df1 = pd.DataFrame(adata_result.uns['rank_velocity_genes']['names'])
df1.to_csv("WT_RankedVelocityGenes.csv", index=False)

scv.tl.rank_velocity_genes(adata_result, groupby='anno_7', min_corr=.3, vkey='velo_s_norm', n_genes=10000)
df2 = pd.DataFrame(adata_result.uns['rank_velocity_genes']['names'])
df2.to_csv("WT_RankedVelocityGenesAnno7.csv", index=False)

#KMT2D-LOF
scv.tl.rank_velocity_genes(adata_result2, groupby='celltype', min_corr=.3, vkey='velo_s_norm', n_genes=10000)
df4 = pd.DataFrame(adata_result2.uns['rank_velocity_genes']['names'])
df4.to_csv("KO_RankedVelocityGenes.csv", index=False)

scv.tl.rank_velocity_genes(adata_result2, groupby='anno_7', min_corr=.3, vkey='velo_s_norm', n_genes=10000)
df5 = pd.DataFrame(adata_result2.uns['rank_velocity_genes']['names'])
df5.to_csv("KO_RankedVelocityGenesAnno7.csv", index=False)
```

    ranking velocity genes
        finished (0:00:00) --> added 
        'rank_velocity_genes', sorted scores by group ids (adata.uns) 
        'spearmans_score', spearmans correlation scores (adata.var)
    ranking velocity genes
        finished (0:00:00) --> added 
        'rank_velocity_genes', sorted scores by group ids (adata.uns) 
        'spearmans_score', spearmans correlation scores (adata.var)
    ranking velocity genes
        finished (0:00:01) --> added 
        'rank_velocity_genes', sorted scores by group ids (adata.uns) 
        'spearmans_score', spearmans correlation scores (adata.var)
    ranking velocity genes
        finished (0:00:00) --> added 
        'rank_velocity_genes', sorted scores by group ids (adata.uns) 
        'spearmans_score', spearmans correlation scores (adata.var)
    ranking velocity genes
        finished (0:00:00) --> added 
        'rank_velocity_genes', sorted scores by group ids (adata.uns) 
        'spearmans_score', spearmans correlation scores (adata.var)
    ranking velocity genes
        finished (0:00:01) --> added 
        'rank_velocity_genes', sorted scores by group ids (adata.uns) 
        'spearmans_score', spearmans correlation scores (adata.var)



```python
gene_list = ['PAX6', 'NEUROD4', 'OLIG1', 'RBPJ', 'CDK6']

# Accessibility/expression by gene time, colored by the four potential states.
# The solid black curve indicates anchors.
fig=mv.dynamic_plot(adata_result, gene_list, color_by='state', axis_on=True, frame_on=False)
plt.savefig("phase_WTmay.pdf", dpi=300, bbox_inches='tight', transparent=True) 

fig=mv.dynamic_plot(adata_result2, gene_list, color_by='state', axis_on=True, frame_on=False)
plt.savefig("phase_KOmay.pdf", dpi=300, bbox_inches='tight', transparent=True) 
```


    
![png](output_41_0.png)
    



    
![png](output_41_1.png)
    



```python
gene_list = ['ELAVL3', 'RBFOX3', 'MAP2', 'SLIT2', 'DCX', 'KMT2C', 'KAT6A', 'NCOR1', 'EZH2']

# Accessibility/expression by gene time, colored by the four potential states.
# The solid black curve indicates anchors.
fig=mv.dynamic_plot(adata_result, gene_list, color_by='state', axis_on=True, frame_on=False)
plt.savefig("phasesup_WTmay.pdf", dpi=300, bbox_inches='tight', transparent=True) 

fig=mv.dynamic_plot(adata_result2, gene_list, color_by='state', axis_on=True, frame_on=False)
plt.savefig("phasesup_KOmay.pdf", dpi=300, bbox_inches='tight', transparent=True) 
```


    
![png](output_42_0.png)
    



    
![png](output_42_1.png)
    



```python

```
