# KMT2D Cerebral Organoid Multiome Analysis

**Lead Lab:** Serrano Lab, Center for Regenerative Medicine (CReM), Boston University 
**Lead Analysis:** Center for Regenerative Medicine Biopinformatics Core. Director: Pushpinder Bawa, PhD. 
**Contact:** [maserr@bu.edu]  
**License:** See [LICENSE](LICENSE) file

## Related Publication

**KMT2D-deficiency destabilizes lineage progression in immature neural progenitors**  
Carly S. Golden<sup>1</sup>, Pushpinder Bawa<sup>1,2</sup>, Feiya Wang<sup>1,2</sup>, Olaf Bodamer<sup>3</sup>, Joseph H. Yost<sup>4</sup>, Maria A. Serrano<sup>1*</sup>  
<sup>1</sup> Center for Regenerative Medicine, Department of Medicine, Boston University Chobanian & Avedisian School of Medicine, Boston, USA  
<sup>2</sup> Boston Medical Center, Boston, Massachusetts, USA  
<sup>3</sup> Department of Pediatrics, Division of Genetics and Genomics, Boston Children’s Hospital  
<sup>4</sup> The University of Utah, Salt Lake City, Utah and The Catholic University of America, Washington DC, USA  
*Corresponding Author: maserr@bu.edu*

## Overview

This repository hosts all code, documentation and processed results for our multi‑omic profiling of human cerebral organoids over early neurodevelopmental time‑points (Day 51 & Day 74).  
It includes fully reproducible pipelines in **R (R Markdown)** and **Python (Jupyter/Markdown)** for:

| Modality | Main notebook | Purpose |
|----------|---------------|---------|
| **WNN integration** | `10x_ATAC_RNA_WNN_CombiningPeaks.Rmd` | Joint ATAC + RNA preprocessing, quality control, weighted‑nearest‑neighbour (WNN) integration and clustering (Seurat v4 + Signac v1.6). |
| **Regulatory network inference** | `figR.Rmd` | Identification of Peak–Gene links, DORCs and TF‑DORC regulatory circuits (FigR). |
| **Chromatin–transcription dynamics** | `MultiVeloAnalysis_CerebralOrganoids.md` | Velocity‑based modelling of coupled RNA/ATAC dynamics (MultiVelo + scVelo). |

**Key objectives:**
- Map chromatin accessibility and transcriptional states in KMT2D-LOF and wild-type (WT) organoids.
- Identify regulatory elements (DORCs) and transcription factors driving cell fate decisions.
- Model dynamic chromatin-transcription coupling using MultiVelo.

## Folder Structure


```
.
├── 10x\_ATAC\_RNA\_WNN\_CombingPeaks.Rmd     # R analysis pipeline (Seurat/Signac WNN)
├── figR.Rmd                                  # R analysis pipeline for FigR (regulatory networks)
├── MultiVeloAnalysis\_CerebralOrganoids.md   # MultiVelo pipeline (Python; gene regulation dynamics)
├── LICENSE
├── README.md
└── .gitignore

```
Raw data, processed objects, and intermediate results are stored on secure institutional servers due to size and privacy.
Data is also deposited in the Gene Expression Omnibus (GEO) data repository. *Reviewer token can be provided upon request.*


## Experimental Methods

### Organoid Dissociation & Multiome Library Preparation

- **Dissociation:** Organoids were pooled by size, dissociated to single-cell suspensions using a papain-based enzymatic protocol with careful pH adjustment, temperature control, and debris removal.
- **Nuclei Isolation:** Adapted 10x Genomics protocol with sequential filtration, hypotonic lysis, and triple wash steps.
- **Library Preparation:** ~5,000 nuclei/sample loaded on Chromium Next GEM Single Cell Multiome (10x Genomics) for simultaneous snATAC-seq and scRNA-seq. Libraries sequenced at >200M reads/sample (each modality).

### Data Preprocessing

- **Alignment & Quantification:** FASTQ generated via Cell Ranger ARC v2.0.0, aligned to GRCh38.
- **Initial Filtering:** Cells with <800 detected genes or >20% mitochondrial reads were excluded.
- **Quality Control:** Used Seurat/Signac to assess metrics (nucleosome signal, TSS enrichment, reads in peaks, blacklist ratio). 
- **Dimensionality Reduction:** PCA for RNA, LSI for ATAC. UMAPs generated for joint and separate modalities. 

## Analysis Pipelines

### 1. Weighted Nearest Neighbor (WNN) Analysis

- **File:** `10x_ATAC_RNA_WNN_CombingPeaks.Rmd`
- **Key Steps:**
    - Merge multiome datasets from multiple time points/genotypes.
    - Construct union peak sets for ATAC.
    - Perform SCTransform (RNA), TF-IDF (ATAC), PCA/LSI, WNN integration, and clustering.
    - Visualize cell states and perform marker gene analysis.
    - Differential expression via MAST, cluster annotation, and UMAP plotting.

### 2. Regulatory Network Inference (FigR)

- **File:** `figR.Rmd`
- **Key Steps:**
    - Subset neural progenitor clusters.
    - Run `runGenePeakcorr()` for peak-gene correlations (±10kb TSS, n=100 permutations).
    - Identify DORC genes (≥7 linked peaks).
    - Smooth DORC scores across KNN graphs.
    - Infer TF-DORC regulatory networks via `runFigRGRN()` and visualize with custom heatmaps.

### 3. Chromatin-Transcription Dynamics (MultiVelo)

- **File:** `MultiVeloAnalysis_CerebralOrganoids.md` (with Python scripts)
- **Key Steps:**
    - Preprocess RNA and ATAC data (Scanpy/scVelo).
    - Intersect barcodes/genes; normalize and smooth data.
    - Model gene regulatory dynamics (induction, repression, state transitions) across pseudotime.
    - Visualize latent time, velocity streams, gene state assignments, and key regulatory genes.
    - All scripts are compatible with conda environments and contain detailed stepwise comments.


## Reproducibility

- **Software:**  
    - R 4.2+, Seurat v4, Signac v1.6, FigR, chromVAR, MAST, patchwork  
    - Python 3.9+, scanpy, scvelo, multivelo, numpy, pandas, matplotlib  
- **Code/Parameter Details:**  
    - See Rmd/MD scripts for exact filtering, clustering, and plotting parameters.
    - All steps are heavily commented for reproducibility and adaptation.


## Data & Results

- Processed count matrices, intermediate RDS and H5AD objects are available on request (pending data sharing policies).
- Plots, DEG tables, and summary statistics are generated automatically to `/plots` and `/rds` directories (paths specified in code).
- Figures illustrate chromatin accessibility, cell states, velocity dynamics, DORC regulatory modules, and TF activity in both genotypes.


## Contributing & Community Guidelines

### How to Use This Repository

1. **Clone the repository:**  
2. **Install required R and Python packages** (see instructions above).
3. **Run R or Python notebooks/scripts** for your analyses of interest.
4. **Use feature branches:**  
   - Create a new branch for any edits or new analyses (`git checkout -b my-feature`).
   - Never commit directly to `main`—this branch serves as the stable, published version for reproducibility and publication.
   - When your changes are ready and fully reproducible, open a pull request (PR) for code review before merging into `main`.

### Contributing Guidelines

- **Branch first:**  
  Always start by creating a new branch (e.g., `git checkout -b feature/my-analysis`).  
  Avoid direct commits to `main` to maintain its stability as the published version.
- **Keep raw data out of the repository:**  
  Use symbolic links or follow `.gitignore` rules to prevent large files from entering version control.
- **Atomic commits:**  
  Make small, logical commits with clear, descriptive messages.
- **Pull requests and review:**  
  Open a PR only when your analysis runs end-to-end and all notebooks/scripts compile without errors.  
  At least one reviewer must approve the PR before merging to `main`.
- **Documentation:**  
  Update README and code comments for any new workflows or significant changes to existing analysis.

### Code of Conduct

We are committed to fostering an open, collaborative, and respectful research environment. All contributors, users, and collaborators of this repository are expected to uphold the highest standards of academic collegiality and integrity.

- **Scientific Rigor and Openness:**  
  All code, data, and analyses should be shared transparently, with clear documentation and appropriate attribution. Constructive feedback, critical review, and reproducibility are essential to our work.
- **Collaboration and Communication:**  
  We encourage open communication and welcome questions, suggestions, and constructive criticism. Disagreements should be handled respectfully and with the shared goal of advancing science.
- **Attribution and Acknowledgement:**  
  Please properly credit original work and contributions, including those of all collaborators and prior publications. Refer to our in-repo acknowledgments and cite the relevant bioRxiv preprint/final publication when using this resource.
- **Ethics and Data Use:**  
  Ensure all data and materials are used in accordance with institutional, legal, and ethical guidelines, including proper handling of sensitive data and privacy considerations.
- **Reporting Issues:**  
  If you experience or witness violations of this Code of Conduct, or have concerns regarding collegiality or data use, please contact the corresponding author ([maserr@bu.edu](mailto:maserr@bu.edu)) or repository maintainers.

We thank all contributors for upholding these values and helping create a supportive and productive research community.

## References

- Buenrostro, J.D., et al. (2018). Integrated Single-Cell Analysis Maps the Continuous Regulatory Landscape of Human Hematopoietic Differentiation. *Cell* 173, 1535–1548.
- Li, X., et al. (2023). MultiVelo: Integration of chromatin accessibility and gene expression dynamics at single-cell resolution. *Nat Methods* 20, 1411–1423. [PMCID: PMC10246490]
- Stuart, T., et al. (2019). Comprehensive Integration of Single-Cell Data. *Cell* 177, 1888–1902.
- For full protocol details, see in-repo method descriptions and [10x Genomics Documentation](https://www.10xgenomics.com/).
 Additional References can be found in the manuscript: [**KMT2D-deficiency destabilizes lineage progression in immature neural progenitors**](LINK)

## Acknowledgements

- Peer review acknowledgement
- Key contributors acknowledgement 


