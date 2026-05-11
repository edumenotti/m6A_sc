# scRNA-seq Processing Pipeline — Charles Dataset Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible Nextflow pipeline that takes the aggregated CellRanger H5 matrix from the Charles dataset (24 samples, mouse hematopoiesis, DMSO vs STM, LSK/LK/Input populations) and produces QC-filtered, batch-integrated, clustered, and annotated AnnData objects with cell type labels validated against canonical markers.

**Architecture:** Six sequential Python scripts (one per analysis stage) are orchestrated by a Nextflow DSL2 pipeline. Each script reads a `.h5ad` from the previous step and writes the next. The pipeline can run locally (pixi env) or on a SLURM cluster (Yale HPC/Grace).

**Tech Stack:** Python 3.11 · scanpy · scvi-tools/scANVI · scArches-style query mapping · popV · Nextflow DSL2 · pixi

---

## State-of-the-Art Report: Tool Selection Rationale

### Integration: Why scVI?
The **Luecken et al. 2022** benchmark (Nature Methods, "Benchmarking atlas-level data integration in single-cell genomics") evaluated 68 method/preprocessing combinations on >1.2 million cells. **scVI and scANVI consistently ranked among the top performers**, especially on complex multi-batch tasks. Key finding: highly variable gene selection before integration improves results; scaling hurts biological variation. This dataset has 6 pools → scVI is the correct choice.

### Cell Type Annotation: Strategy for Specialized Populations

This is the most important decision point. The data contains only hematopoietic progenitor populations (LSK = HSC/MPP, LK = progenitors, Input = whole BM). General databases like the full Human Cell Atlas contain 500+ types, making direct application noisy.

**The practical solution for specialized populations:**

**Option A — Use a domain-specific reference (recommended here):**
The **Nestorowa et al. 2016** (Blood) dataset contains 1,656 mouse HSPCs with detailed labels (LT-HSC, ST-HSC, MPP, LMPP, CMP, GMP, MEP, CLP, MkP, etc.). It is the canonical mouse HSPC reference, specifically built from LSK and LK sorted populations — identical to our experimental gates. Available via Bioconductor `scRNAseq` package. **No subsetting needed** because the reference is already specialized.

**Option B — Subset a large atlas:**
If using the **Meta-Analytic Mouse Bone Marrow Atlas** (Gillis lab, CSHL; 300,000+ cells across 12 datasets), you pre-filter to hematopoietic cell types: `adata_ref = adata_ref[adata_ref.obs['cell_type'].isin(hspc_types)]`. This is the "subset" workflow the user asked about.

**Option C — popV ensemble (2024, Nature Genetics):**
popV runs 8 annotation algorithms in parallel (KNN, scVI, SVM, logistic regression, etc.) and takes a consensus vote weighted by Cell Ontology hierarchy. It provides uncertainty scores per cell. This is the most robust approach and handles cases where your query cells don't match any reference cell type (they get flagged as uncertain). **We will use popV as the primary annotator.**

**Option D — CellTypist:**
CellTypist is disabled for this dataset. The human `Immune_All_Low.pkl` model matched only 7 query genes in the local test, so its labels are not reliable here.

**Recommended strategy (this plan):**
1. **popV** with Nestorowa 2016 reference → primary labels + uncertainty scores
2. **scANVI/scArches-style reference mapping** → secondary label transfer
3. **Manual marker dotplot and marker-score consistency tables** → ground truth validation

### Doublet Detection: DoubletFinder via Python-driven R bridge
Use **DoubletFinder** for doublet detection because this project should match the Seurat/R ecosystem used in the original processing logic. The pipeline entry point remains Python (`pipeline/scripts/02_doublets.py`) so the Nextflow workflow stays Python-oriented, but the script exports each sample's sparse raw counts to temporary Matrix Market files and calls `Rscript` to run Seurat + DoubletFinder. This avoids brittle in-memory AnnData→Seurat conversion while still satisfying the requirement to run the R package from a Python-controlled stage. Run **per-sample** (`sample_id`) to avoid synthetic cross-sample doublets that cannot exist biologically.

### Ambient RNA: Skipped (data limitation)
CellBender requires the **raw (unfiltered)** feature-barcode matrix to model ambient RNA. This dataset only provides the `filtered_feature_bc_matrix.h5`. SoupX can work with the filtered matrix but requires the cluster topology — it would be circular to run QC before clustering. **Decision: skip ambient RNA correction.** Note this in Methods.

---

## File Structure

```
Charles/
├── pipeline/
│   ├── main.nf                  # Nextflow entry point (DSL2)
│   ├── nextflow.config          # Executor config (local + SLURM profiles)
│   ├── params.yaml              # All tunable parameters
│   ├── conf/
│   │   ├── local.config         # pixi local profile
│   │   └── slurm.config         # Yale HPC SLURM profile
│   ├── modules/
│   │   ├── qc.nf
│   │   ├── doublets.nf
│   │   ├── normalize.nf
│   │   ├── integrate.nf
│   │   ├── cluster.nf
│   │   ├── annotate.nf
│   │   └── markers.nf
│   └── scripts/
│       ├── 01_qc.py
│       ├── 02_doublets.py
│       ├── 03_normalize.py
│       ├── 04_integrate.py
│       ├── 05_cluster.py
│       ├── 06_annotate.py
│       └── 07_markers.py
├── data/
│   └── count/filtered_feature_bc_matrix.h5  (existing)
├── results/
│   ├── 00_qc_scan/
│   ├── 01_qc/
│   ├── 02_doublets/
│   ├── 03_hvg_scan/
│   ├── 03_normalize/
│   ├── 04_integrate/
│   ├── 05_cluster/
│   ├── 06_annotate/
│   └── 07_markers/
└── docs/
    └── superpowers/plans/  (this file)
```

---

## Task 0: QC Parameter Scan (`scripts/00_qc_parameter_scan.py`)

**Files:**
- Create: `pipeline/scripts/00_qc_parameter_scan.py`
- Writes: `results/00_qc_scan/`

### What this script does
- Loads the raw aggregated H5 and annotates all cells by sample.
- Computes `% mitochondrial`, `% ribosomal` (`Rpl*`, `Rps*`), and `% hemoglobin` (`Hba*`, `Hbb*`) metrics.
- Produces per-sample QC summaries and plots.
- Scans combinations of `count_nmads`, `mito_nmads`, optional absolute mitochondrial caps, and optional ribosomal caps.
- Writes `qc_filter_parameter_scan.csv` so the filtering choice can be justified before integration.

- [x] **Step 0.1: Run QC scan**

```bash
cd /home/edu-pc/Yale/Charles
pixi run python pipeline/scripts/00_qc_parameter_scan.py \
  --h5 data/count/filtered_feature_bc_matrix.h5 \
  --out results/00_qc_scan/
```

Observed in the first run:
- 13 mitochondrial genes
- 101 ribosomal genes
- 8 hemoglobin genes

### Final QC decision
- Ambient RNA correction: skipped because only filtered CellRanger matrices are available locally.
- Cell filtering: `count_nmads=6`, `mito_nmads=4`, no absolute ribosomal cutoff.
- Rationale: compared with the stricter `count_nmads=5`, `mito_nmads=3` setting, 6/4 retains more cells while keeping post-filter mitochondrial burden nearly unchanged.
- Final rerun output: 53,494 cells before QC, 8,239 removed, 45,255 retained.

---

## Task 1: QC Script (`scripts/01_qc.py`)

**Files:**
- Create: `pipeline/scripts/01_qc.py`
- Create: `results/01_qc/` (runtime)

### What this script does
- Loads the H5 file
- Separates GEX from HTO layers
- Annotates cells with sample metadata
- Computes QC metrics (n_genes, n_counts, pct_mito)
- Computes ribosomal and hemoglobin QC metrics for reporting
- Generates QC violin/scatter plots
- Filters cells using MAD-based thresholds (not fixed cutoffs)
- Saves filtered AnnData

### Why MAD-based thresholds?
Fixed thresholds (e.g., "remove cells with < 200 genes") are arbitrary and dataset-dependent. MAD (median absolute deviation) thresholds adapt to each sample's distribution. After scanning candidate settings, use a moderately permissive final setting: flag cells that are >6 MADs below the median for n_genes or n_counts, or >4 MADs above the median for pct_mito. This retained substantially more cells than 5/3 while keeping the post-filter mitochondrial burden nearly unchanged.

- [ ] **Step 1.1: Create the QC script**

```python
# pipeline/scripts/01_qc.py
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def mad_filter(values: np.ndarray, nmads: float = 5, direction: str = "both") -> np.ndarray:
    """Return boolean mask of cells to KEEP (True = keep)."""
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    if direction in ("lower", "both"):
        low = med - nmads * mad
    else:
        low = -np.inf
    if direction in ("upper", "both"):
        high = med + nmads * mad
    else:
        high = np.inf
    return (values >= low) & (values <= high)

SAMPLE_MAP = {
    "1":  ("D1", "DMSO", "LSK", "2"), "2":  ("D1", "DMSO", "LSK", "1"),
    "3":  ("D1", "DMSO", "LK",  "2"), "4":  ("D1", "DMSO", "LK",  "1"),
    "5":  ("D1", "STM",  "LSK", "2"), "6":  ("D1", "STM",  "LSK", "1"),
    "7":  ("D1", "STM",  "LK",  "2"), "8":  ("D1", "STM",  "LK",  "1"),
    "9":  ("D1", "STM",  "I",   "2"), "10": ("D1", "STM",  "I",   "1"),
    "11": ("D1", "DMSO", "I",   "2"), "12": ("D1", "DMSO", "I",   "1"),
    "13": ("D2", "STM",  "LSK", "2"), "14": ("D2", "STM",  "LSK", "1"),
    "15": ("D2", "DMSO", "LSK", "2"), "16": ("D2", "DMSO", "LSK", "1"),
    "17": ("D2", "STM",  "LK",  "2"), "18": ("D2", "STM",  "LK",  "1"),
    "19": ("D2", "DMSO", "LK",  "2"), "20": ("D2", "DMSO", "LK",  "1"),
    "21": ("D2", "STM",  "I",   "2"), "22": ("D2", "STM",  "I",   "1"),
    "23": ("D2", "DMSO", "I",   "2"), "24": ("D2", "DMSO", "I",   "1"),
}

def main(h5_path: str, out_dir: str, mito_nmads: float = 3.0, count_nmads: float = 5.0):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Load
    adata = sc.read_10x_h5(h5_path, gex_only=False)
    adata.var_names_make_unique()

    # Separate GEX and HTO
    is_hto = adata.var["feature_types"] == "Multiplexing Capture"
    hto = adata[:, is_hto].copy()
    adata = adata[:, ~is_hto].copy()

    # Store HTO in obsm
    adata.obsm["HTO"] = hto.X.toarray()
    adata.uns["HTO_names"] = hto.var_names.tolist()

    # Add sample metadata
    suffixes = [b.split("-")[-1] for b in adata.obs_names]
    meta = pd.DataFrame(
        [SAMPLE_MAP.get(s, ("unknown",) * 4) for s in suffixes],
        columns=["donor", "treatment", "population", "replicate"],
        index=adata.obs_names,
    )
    adata.obs = meta
    adata.obs["sample_id"] = [
        f"{r.donor}_{r.treatment}_{r.population}_45_{r.replicate}"
        for _, r in adata.obs.iterrows()
    ]
    adata.obs["pool"] = adata.obs["sample_id"].map(lambda x: _pool_from_sample(x))

    # QC metrics
    mito_genes = adata.var_names.str.startswith("mt-")
    adata.var["mt"] = mito_genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # MAD filters (per-sample)
    keep = np.ones(adata.n_obs, dtype=bool)
    for sid in adata.obs["sample_id"].unique():
        mask = adata.obs["sample_id"] == sid
        idx = np.where(mask)[0]
        keep[idx] &= mad_filter(adata.obs.loc[mask, "log1p_n_genes_by_counts"].values, count_nmads)
        keep[idx] &= mad_filter(adata.obs.loc[mask, "log1p_total_counts"].values, count_nmads)
        keep[idx] &= mad_filter(adata.obs.loc[mask, "pct_counts_mt"].values, mito_nmads, "upper")

    # QC plots before filtering
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                 groupby="sample_id", rotation=90, show=False)
    plt.savefig(out / "qc_violin_before.png", dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Cells before QC: {adata.n_obs}")
    print(f"Cells removed: {(~keep).sum()}")
    adata = adata[keep].copy()
    print(f"Cells after QC: {adata.n_obs}")

    # QC plots after filtering
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                 groupby="sample_id", rotation=90, show=False)
    plt.savefig(out / "qc_violin_after.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Save cells-per-sample summary
    summary = adata.obs.groupby("sample_id").size().rename("n_cells").reset_index()
    summary.to_csv(out / "cells_per_sample.csv", index=False)

    adata.write_h5ad(out / "adata_qc.h5ad")
    print(f"Saved: {out / 'adata_qc.h5ad'}")

POOL_MAP = {
    "D1_DMSO_LSK": "Pool_A", "D1_DMSO_LK": "Pool_A",
    "D1_STM_LSK":  "Pool_B", "D1_STM_LK":  "Pool_B",
    "D1_STM_I":    "Pool_C", "D1_DMSO_I":  "Pool_C",
    "D2_STM_LSK":  "Pool_D", "D2_DMSO_LSK":"Pool_D",
    "D2_STM_LK":   "Pool_E", "D2_DMSO_LK": "Pool_E",
    "D2_STM_I":    "Pool_F", "D2_DMSO_I":  "Pool_F",
}

def _pool_from_sample(sample_id: str) -> str:
    key = "_".join(sample_id.split("_")[:3])
    return POOL_MAP.get(key, "unknown")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--mito_nmads", type=float, default=4.0)
    parser.add_argument("--count_nmads", type=float, default=6.0)
    args = parser.parse_args()
    main(args.h5, args.out, args.mito_nmads, args.count_nmads)
```

- [ ] **Step 1.2: Run to verify it works**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/01_qc.py \
  --h5 data/count/filtered_feature_bc_matrix.h5 \
  --out results/01_qc/ \
  --count_nmads 6 \
  --mito_nmads 4
```

Expected output:
```
Cells before QC: 53494
Cells removed: ~8000
Cells after QC: ~45000
Saved: results/01_qc/adata_qc.h5ad
```

- [ ] **Step 1.3: Commit**

```bash
git add pipeline/scripts/01_qc.py
git commit -m "feat: QC script with MAD-based per-sample filtering"
```

---

## Task 2: Doublet Detection with DoubletFinder (`scripts/02_doublets.py`)

**Files:**
- Create: `pipeline/scripts/02_doublets.py`
- Reads: `results/01_qc/adata_qc.h5ad`
- Writes: `results/02_doublets/adata_no_doublets.h5ad`

### Why per-sample DoubletFinder?
DoubletFinder creates artificial doublets inside a Seurat object and scores cells against that synthetic neighborhood. If run on the full aggregated matrix, it can create impossible cross-sample doublets because each cell barcode suffix corresponds to one sample. Running per `sample_id` produces realistic synthetic doublets and keeps the expected doublet rate interpretable.

- [ ] **Step 2.1: Create the doublet script**

```python
# pipeline/scripts/02_doublets.py
"""Run per-sample DoubletFinder from a Python pipeline stage.

Python owns AnnData I/O and writes one sparse Matrix Market count matrix per
sample. Rscript owns Seurat + DoubletFinder and writes back a CSV with barcode,
score, class, pK, and expected doublet counts. This keeps the pipeline runnable
from Python/Nextflow without relying on fragile in-memory rpy2 conversion.
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import subprocess
import tempfile
from pathlib import Path

# Full implementation writes temporary sparse matrices and invokes:
# Rscript pipeline/scripts/run_doubletfinder.R --counts ... --genes ... --cells ...
# The R helper loads Seurat + DoubletFinder, optionally installs DoubletFinder
# with remotes::install_github("chris-mcginnis-ucsf/DoubletFinder"), runs
# NormalizeData, FindVariableFeatures, ScaleData, RunPCA, FindNeighbors,
# FindClusters, modelHomotypic, and doubletFinder/doubletFinder_v3.

def main(in_path: str, out_dir: str, expected_rate: float = 0.075,
         rscript: str = "Rscript", install_missing: bool = False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    # The executable script in the repository contains the full bridge.
    # It writes `doubletfinder_calls.csv`, `doubletfinder_summary.csv`,
    # `doublet_scores.png`, and `adata_no_doublets.h5ad`.
    ...
```

- [ ] **Step 2.2: Add R bridge dependencies to pixi**

Add to `pixi.toml` under `[dependencies]`:
```toml
r-base = ">=4.5,<4.6"
r-seurat = "*"
r-seuratobject = "*"
r-matrix = "*"
r-remotes = "*"
rpy2 = "*"
```

Then run: `cd /home/edu-pc/Yale/Charles && pixi install`

- [ ] **Step 2.3: Install DoubletFinder if missing**

DoubletFinder is not available as `r-doubletfinder` in conda-forge. Install it from GitHub into the pixi R environment:

```bash
cd /home/edu-pc/Yale/Charles
pixi run Rscript -e 'if (!requireNamespace("DoubletFinder", quietly=TRUE)) remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade="never")'
```

- [ ] **Step 2.4: Run**

```bash
cd /home/edu-pc/Yale/Charles
pixi run python pipeline/scripts/02_doublets.py \
  --input results/01_qc/adata_qc.h5ad \
  --out results/02_doublets/ \
  --expected-rate 0.075 \
  --pk 0.09 \
  --install-missing
```

Expected: per-sample DoubletFinder scores and calls in `doubletfinder_calls.csv`, a per-sample `doubletfinder_summary.csv`, a score distribution plot, and ~5-10% doublets removed depending on recovered cells/sample.

For a slower tuning run, add `--auto-pk` to perform DoubletFinder's pK sweep per sample. The default local execution uses `--pk 0.09` to keep the first complete pipeline run tractable.

- [ ] **Step 2.5: Commit**

```bash
git add pipeline/scripts/02_doublets.py pipeline/scripts/run_doubletfinder.R pixi.toml pixi.lock
git commit -m "feat: per-sample doublet detection with DoubletFinder"
```

---

## Task 3: Normalization & Feature Selection (`scripts/03_normalize.py`)

**Files:**
- Create: `pipeline/scripts/03_normalize.py`
- Reads: `results/02_doublets/adata_no_doublets.h5ad`
- Writes: `results/03_normalize/adata_normalized.h5ad`

### Design decisions
- Store raw counts in `adata.layers["counts"]` before any normalization (required by scVI)
- Use `sc.pp.normalize_total` + `sc.pp.log1p` for visualization layers only
- Select 3000 highly variable genes (HVGs) using the `seurat_v3` method (variance-stabilizing, works on counts directly)
- Exclude mitochondrial, ribosomal, and hemoglobin genes from HVG list (they are QC/technical-dominant signals here, not desired integration features)

### HVG parameter scan
Before finalizing `n_hvgs`, run:

```bash
cd /home/edu-pc/Yale/Charles
pixi run python pipeline/scripts/03_hvg_parameter_scan.py \
  --input results/02_doublets/adata_no_doublets.h5ad \
  --out results/03_hvg_scan/ \
  --n-hvgs-values 1000,2000,3000,4000,5000
```

Observed in the first run:
- 3,000 requested HVGs produced 2,915 effective HVGs after excluding 10 mitochondrial, 71 ribosomal, and 4 hemoglobin genes.
- 2,000 requested HVGs produced 1,918 effective HVGs; 4,000 requested HVGs produced 3,910 effective HVGs.

Final rerun after QC 6/4 and DoubletFinder retained 2,909 effective HVGs after excluding 10 mitochondrial, 77 ribosomal, and 4 hemoglobin genes.

- [ ] **Step 3.1: Create the normalization script**

```python
# pipeline/scripts/03_normalize.py
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def main(in_path: str, out_dir: str, n_hvgs: int = 3000):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    # Keep raw counts as a layer (required by scVI)
    adata.layers["counts"] = adata.X.copy()

    # Normalize for visualization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    # HVG selection on counts (seurat_v3 uses raw counts via layer)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvgs,
        flavor="seurat_v3",
        layer="counts",
        batch_key="sample_id",  # select HVGs that are variable across samples
        span=0.3,
    )

    # Exclude technical/QC-dominant genes from HVGs
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))
    adata.var["ribo"] = adata.var_names.str.match(r"^(Rpl|Rps)")
    adata.var["hb"] = adata.var_names.str.match(r"^(Hba|Hbb)")
    adata.var.loc[adata.var["mt"] | adata.var["ribo"] | adata.var["hb"], "highly_variable"] = False

    n_hvg = adata.var["highly_variable"].sum()
    print(f"HVGs selected: {n_hvg}")

    # Plot HVG dispersion
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(out / "hvg_dispersion.png", dpi=150, bbox_inches="tight")
    plt.close()

    # PCA for visualization only (not used for integration — scVI handles that)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
    plt.savefig(out / "pca_variance_ratio.png", dpi=150)
    plt.close()

    adata.write_h5ad(out / "adata_normalized.h5ad")
    print(f"Saved: {out / 'adata_normalized.h5ad'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n_hvgs", type=int, default=3000)
    args = parser.parse_args()
    main(args.input, args.out, args.n_hvgs)
```

- [ ] **Step 3.2: Run**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/03_normalize.py \
  --input Charles/results/02_doublets/adata_no_doublets.h5ad \
  --out Charles/results/03_normalize/
```

Expected: `HVGs selected: 3000`

- [ ] **Step 3.3: Commit**

```bash
git add pipeline/scripts/03_normalize.py
git commit -m "feat: normalization, log1p, HVG selection (seurat_v3 + batch_key)"
```

---

## Task 4: Batch Integration with scVI (`scripts/04_integrate.py`)

**Files:**
- Create: `pipeline/scripts/04_integrate.py`
- Reads: `results/03_normalize/adata_normalized.h5ad`
- Writes: `results/04_integrate/adata_integrated.h5ad`

### Design decisions
- **Batch key**: `sample_id` (24 samples) — the most granular. Using `pool` (6 pools) also valid; `sample_id` is more conservative.
- **Categorical covariates**: `donor`, `treatment`, `population` — tell scVI about biological covariates to preserve
- **n_latent**: 30 (default 10 is too low for complex datasets; 30 captures more biology)
- **n_layers**: 2 (deeper encoder/decoder)
- GPU training if available (pixi.toml has CUDA 12.8)

- [x] **Step 4.1: Create the integration script**

```python
# pipeline/scripts/04_integrate.py
import scanpy as sc
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def main(in_path: str, out_dir: str,
         n_latent: int = 30, n_layers: int = 2, max_epochs: int = 400):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    # scVI setup: use counts layer, HVGs only
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="sample_id",
        categorical_covariate_keys=["donor", "treatment", "population"],
    )

    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood="nb",  # negative binomial for UMI counts
    )

    model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        early_stopping_patience=20,
        batch_size=256,
    )

    # Save training history
    fig, ax = plt.subplots()
    train_loss = model.history["elbo_train"]
    val_loss = model.history["elbo_validation"]
    ax.plot(train_loss, label="train ELBO")
    ax.plot(val_loss, label="validation ELBO")
    ax.set_xlabel("Epoch")
    ax.legend()
    plt.savefig(out / "training_loss.png", dpi=150)
    plt.close()

    # Extract latent representation
    adata.obsm["X_scVI"] = model.get_latent_representation()

    # UMAP on integrated space
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30, metric="euclidean")
    sc.tl.umap(adata, min_dist=0.3)

    # Plot UMAP colored by batch variables
    for color in ["sample_id", "pool", "donor", "treatment", "population"]:
        sc.pl.umap(adata, color=color, show=False)
        plt.savefig(out / f"umap_{color}.png", dpi=150, bbox_inches="tight")
        plt.close()

    # Save model
    model.save(str(out / "scvi_model"), overwrite=True)

    adata.write_h5ad(out / "adata_integrated.h5ad")
    print(f"Saved: {out / 'adata_integrated.h5ad'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n_latent", type=int, default=30)
    parser.add_argument("--n_layers", type=int, default=2)
    parser.add_argument("--max_epochs", type=int, default=400)
    args = parser.parse_args()
    main(args.input, args.out, args.n_latent, args.n_layers, args.max_epochs)
```

- [x] **Step 4.2: Run (GPU recommended, ~15-30 min)**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/04_integrate.py \
  --input Charles/results/03_normalize/adata_normalized.h5ad \
  --out Charles/results/04_integrate/ \
  --n_latent 30 --max_epochs 400
```

Expected: ELBO converging, UMAP without obvious pool-based separation.

- [x] **Step 4.3: Check integration quality**

Inspect `results/04_integrate/umap_pool.png`. If cells cluster strongly by pool (not biology), increase `max_epochs` or adjust `n_latent`.

Run result: GPU training used CUDA on the NVIDIA GeForce RTX 5060 Laptop GPU and stopped by early stopping at epoch 314/400 with best validation ELBO 6732.945. The integrated AnnData has 42,390 cells, `X_scVI` shape `(42390, 30)`, `X_umap` shape `(42390, 2)`, and a neighbors graph. Visual QC note: `pool` is partially confounded with population/donor/treatment, so do not interpret `umap_pool.png` alone as pure batch separation; donor and treatment are broadly mixed while the dominant structure follows population.

- [ ] **Step 4.4: Commit**

```bash
git add pipeline/scripts/04_integrate.py
git commit -m "feat: scVI batch integration with 30 latent dims, per-sample batch key"
```

---

## Task 5: Clustering (`scripts/05_cluster.py`)

**Files:**
- Create: `pipeline/scripts/05_cluster.py`
- Reads: `results/04_integrate/adata_integrated.h5ad`
- Writes: `results/05_cluster/adata_clustered.h5ad`

- [x] **Step 5.1: Create the clustering script**

```python
# pipeline/scripts/05_cluster.py
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def main(in_path: str, out_dir: str, resolution: float = 0.5):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    # Leiden clustering (Leiden > Louvain — more stable, same complexity)
    sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_r{resolution}")

    # Also run at multiple resolutions for comparison
    for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_r{res}")

    # UMAP plots per resolution
    for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
        sc.pl.umap(adata, color=f"leiden_r{res}", legend_loc="on data", show=False)
        plt.savefig(out / f"umap_leiden_r{res}.png", dpi=150, bbox_inches="tight")
        plt.close()

    # UMAP with population (ground truth gate) overlaid
    sc.pl.umap(adata, color=["population", "treatment", "donor"], show=False)
    plt.savefig(out / "umap_metadata.png", dpi=150, bbox_inches="tight")
    plt.close()

    adata.write_h5ad(out / "adata_clustered.h5ad")
    n_clusters = adata.obs[f"leiden_r{resolution}"].nunique()
    print(f"Leiden r={resolution}: {n_clusters} clusters")
    print(f"Saved: {out / 'adata_clustered.h5ad'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--resolution", type=float, default=0.5)
    args = parser.parse_args()
    main(args.input, args.out, args.resolution)
```

- [x] **Step 5.2: Run**

```bash
cd /home/edu-pc/Yale/Charles
pixi run python pipeline/scripts/05_cluster.py \
  --input results/04_integrate/adata_integrated.h5ad \
  --out results/05_cluster/ \
  --resolution 0.5 \
  --resolutions 0.3,0.5,0.8,1.0,1.5,2.0,2.5,3.0
```

Run result:
- r=0.3: 10 clusters
- r=0.5: 14 clusters
- r=0.8: 20 clusters
- r=1.0: 23 clusters
- r=1.5: 27 clusters
- r=2.0: 36 clusters
- r=2.5: 40 clusters
- r=3.0: 46 clusters

The script writes `cluster_resolution_summary.csv`, `cluster_dominance_summary.csv`, `cluster_composition_long.csv`, and per-resolution count/proportion tables for `pool`, `sample_id`, `population`, `donor`, and `treatment`. At r=0.5, clusters 10 and 12 are >80% dominated by `Pool_C` / `D1_DMSO_I_45_2`; both are mostly `Input` and should be checked during marker validation.

Decision: use `leiden_r0.5` as the primary initial clustering resolution. It gives 14 clusters, avoids the overfragmentation seen at r=2.0-3.0, and still separates interpretable major compartments. Keep r=0.8/r=1.0 as secondary resolutions if marker validation shows mixed populations inside a r=0.5 cluster.

- [ ] **Step 5.3: Commit**

```bash
git add pipeline/scripts/05_cluster.py
git commit -m "feat: Leiden clustering at multiple resolutions"
```

---

## Task 6: Cell Type Annotation (`scripts/06_annotate.py`)

**Files:**
- Create: `pipeline/scripts/06_annotate.py`
- Reads: `results/05_cluster/adata_clustered.h5ad`
- Writes: `results/06_annotate/adata_annotated.h5ad`

### Strategy: popV with mouse hematopoietic reference

**Step-by-step annotation workflow:**

1. Download Nestorowa 2016 reference via `scrnaseq` R package → convert to AnnData (or use the pre-built version from scvi-hub)
2. Run popV with the reference
3. Run scANVI/scArches-style label transfer as the secondary method
4. Validate labels manually against an audited marker table

**About the reference subset question (practical explanation):**

The Nestorowa 2016 reference contains ONLY mouse HSPCs (HSC, MPP, LMPP, CMP, GMP, MEP, CLP, MkP, ErP). There are no other cell types. So for our LSK/LK/Input populations this reference is ideal — no subsetting needed for the HSPC populations.

For the Input (whole BM) fraction, there will be mature immune cells (T cells, B cells, NK cells, neutrophils) that are not fully covered by Nestorowa. For these, popV/scANVI may return an "unknown" or the closest HSPC match with low agreement score. These cells must be resolved by marker validation, not by CellTypist.

CellTypist is disabled in this pipeline. In this dataset, `Immune_All_Low.pkl` matched only 7 of 6,639 model genes, making its labels unusable.

**Practical subsetting if you need it (for other datasets):**
```python
# Example: subset the meta-analytic atlas to only HSPC types
hspc_types = ["HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "MkP", "ErP", "DC_progenitor"]
ref_sub = ref[ref.obs["cell_type"].isin(hspc_types)].copy()
# Then use ref_sub as the reference for popV/scArches
```

- [ ] **Step 6.1: Prepare Nestorowa reference (run once)**

```python
cd /home/edu-pc/Yale/Charles
pixi run Rscript pipeline/scripts/prepare_nestorowa_reference.R \
  --out data/references/nestorowa_2016_hspc.h5ad \
  --label-col auto \
  --install-missing
```

- [x] **Step 6.2: Create the annotation script**

```python
`pipeline/scripts/06_annotate.py` now implements:
- popV primary annotation, using Nestorowa labels and excluding CellTypist from the default popV method list.
- scANVI secondary label transfer on the same reference/query gene intersection.
- per-cluster summaries for `popv_prediction`, `popv_majority_vote_prediction`, and `scanvi_label`.
See the script for executable details.
```

- [ ] **Step 6.3: Run popV + scANVI annotation**

```bash
cd /home/edu-pc/Yale/Charles
pixi run -e popv python pipeline/scripts/06_annotate.py \
  --input results/05_cluster/adata_clustered.h5ad \
  --ref data/references/nestorowa_2016_hspc.h5ad \
  --out results/06_annotate/ \
  --ref-labels-key cell_type \
  --ref-batch-key reference_batch \
  --query-batch-key sample_id \
  --cluster-key leiden_r0.5
```

Previous Paul15 KNN output is retained only as a provisional fallback and should not be treated as the final annotation.

- [ ] **Step 6.5: Commit**

```bash
git add pipeline/scripts/06_annotate.py pipeline/scripts/prepare_nestorowa_reference.R pipeline/config/cell_type_markers.tsv
git commit -m "feat: popV primary annotation with scANVI secondary transfer"
```

---

## Task 7: Marker Validation (`scripts/07_markers.py`)

**Files:**
- Create: `pipeline/scripts/07_markers.py`
- Reads: `results/06_annotate/adata_annotated.h5ad`
- Writes: `results/07_markers/` (plots + CSVs)

### Canonical mouse HSPC markers for manual validation

| Population | Markers (positive) | Markers (negative) |
|-----------|-------------------|-------------------|
| LT-HSC | Slamf1 (CD150), Kit, Ly6a (Sca1) | Cd48, Cd34 |
| ST-HSC / MPP | Slamf1 (low), Kit, Ly6a | — |
| CMP | Cd34, Kit | Fcgr3 (CD16/32) |
| GMP | Cd34, Kit, Fcgr3 | — |
| MEP | Cd34 (low), Kit | Fcgr3 |
| CLP | Il7r (CD127) | Kit (low) |
| Cycling/proliferating | Mki67, Top2a, Pcna | — |

- [x] **Step 7.1: Create the marker validation script**

```python
# pipeline/scripts/07_markers.py
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

HSPC_MARKERS = {
    "HSC/MPP": ["Kit", "Ly6a", "Slamf1", "Cd150", "Cd48"],
    "LT-HSC": ["Slamf1", "Kit", "Ly6a"],
    "CMP_GMP_MEP": ["Cd34", "Fcgr3", "Fcgr2b"],
    "CLP": ["Il7r", "Dntt", "Rag1"],
    "Erythroid": ["Gypa", "Hba-a1", "Klf1", "Gata1"],
    "Megakaryocyte": ["Pf4", "Itga2b", "Vwf"],
    "Granulocyte": ["Elane", "Mpo", "Cebpa"],
    "Monocyte": ["Csf1r", "Cx3cr1", "Ly6c2"],
    "Proliferating": ["Mki67", "Top2a", "Pcna"],
    "Basophil": ["Mcpt8", "Il4"],
    "Dendritic": ["Siglech", "Bst2", "Flt3"],
    "STM_response": ["Cebpb", "Fos", "Jun", "Atf3"],  # generic stress/drug response
}

def main(in_path: str, out_dir: str, cluster_key: str = "leiden_r0.5"):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    # Filter to genes present in data
    markers_flat = [g for genes in HSPC_MARKERS.values() for g in genes]
    present = [g for g in markers_flat if g in adata.var_names]
    missing = [g for g in markers_flat if g not in adata.var_names]
    if missing:
        print(f"Missing markers (not in dataset): {missing}")

    # Dot plot: clusters × markers
    sc.pl.dotplot(
        adata,
        var_names=HSPC_MARKERS,
        groupby=cluster_key,
        use_raw=False,
        layer="lognorm",
        standard_scale="var",
        show=False,
    )
    plt.savefig(out / "dotplot_clusters_markers.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Dot plot: paul15 labels × markers
    sc.pl.dotplot(
        adata,
        var_names=HSPC_MARKERS,
        groupby="paul15_label",
        use_raw=False,
        layer="lognorm",
        standard_scale="var",
        show=False,
    )
    plt.savefig(out / "dotplot_paul15_markers.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Violin plots for key markers per population gate
    sc.pl.violin(
        adata,
        keys=["Kit", "Ly6a", "Slamf1", "Cd48"],
        groupby="population",
        show=False,
    )
    plt.savefig(out / "violin_gate_markers.png", dpi=150, bbox_inches="tight")
    plt.close()

    # UMAP colored by individual marker expression
    for gene in ["Kit", "Ly6a", "Slamf1", "Mki67", "Il7r", "Gata1", "Mpo"]:
        if gene in adata.var_names:
            sc.pl.umap(adata, color=gene, layer="lognorm", show=False, vmax="p99")
            plt.savefig(out / f"umap_expr_{gene}.png", dpi=150, bbox_inches="tight")
            plt.close()

    # Differential expression per cluster (top 10 markers per cluster)
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon",
                            layer="lognorm", use_raw=False)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
    plt.savefig(out / "dotplot_top_de_genes.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Save DE results
    import pandas as pd
    de_results = sc.get.rank_genes_groups_df(adata, group=None)
    de_results.to_csv(out / "de_genes_per_cluster.csv", index=False)

    print(f"Saved marker validation plots to {out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cluster_key", default="leiden_r0.5")
    args = parser.parse_args()
    main(args.input, args.out, args.cluster_key)
```

- [x] **Step 7.2: Run**

```bash
cd /home/edu-pc/Yale/Charles
pixi run python pipeline/scripts/07_markers.py \
  --input results/05_cluster/adata_clustered.h5ad \
  --out results/07_markers/ \
  --cluster_key leiden_r0.5 \
  --watchlist-clusters 10,12
```

Run result: marker/QC validation was run before automated reference annotation so the `Pool_C`/sample-dominated clusters could be assessed directly. Cluster 10 is B-like (`Cd79a`, `Ighm`, `Cd79b`, `Ebf1`, `Pax5`) and cluster 12 is T-like (`Cd3d`, `Cd3g`, `Skap1`, `Ms4a4b`). Their median mitochondrial percentages are low (~2.65% and ~2.47%), so they are not obvious low-quality clusters. Treat them as possible mature lymphoid/input contamination or real lineage cells pending reference annotation.

- [ ] **Step 7.3: Commit**

```bash
git add pipeline/scripts/07_markers.py
git commit -m "feat: marker validation dotplots, violin plots, and DE per cluster"
```

---

## Task 8: Nextflow Pipeline (`pipeline/main.nf`)

**Files:**
- Create: `pipeline/main.nf`
- Create: `pipeline/nextflow.config`
- Create: `pipeline/params.yaml`
- Create: `pipeline/conf/local.config`
- Create: `pipeline/conf/slurm.config`
- Create: `pipeline/modules/qc.nf` through `markers.nf`

### Why Nextflow DSL2?
- Automatic resume (`-resume`) — if annotation fails, reruns only annotation, not QC
- Process-level resource management (GPU for scVI, CPU for QC)
- Portable: same pipeline runs locally and on Yale HPC/Grace cluster

- [ ] **Step 8.1: Create modules**

```groovy
// pipeline/modules/qc.nf
process QC {
    tag "qc"
    publishDir "${params.outdir}/01_qc", mode: 'copy'
    memory '16 GB'
    cpus 4

    input:
    path h5_file

    output:
    path "adata_qc.h5ad", emit: h5ad
    path "*.png"
    path "cells_per_sample.csv"

    script:
    """
    python ${projectDir}/scripts/01_qc.py \
        --h5 ${h5_file} \
        --out . \
        --mito_nmads ${params.mito_nmads} \
        --count_nmads ${params.count_nmads}
    """
}
```

```groovy
// pipeline/modules/doublets.nf
process DOUBLETS {
    tag "doublets"
    publishDir "${params.outdir}/02_doublets", mode: 'copy'
    memory '16 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_no_doublets.h5ad", emit: h5ad
    path "*.png"

    script:
    """
    python ${projectDir}/scripts/02_doublets.py \
        --input ${h5ad} \
        --out .
    """
}
```

```groovy
// pipeline/modules/normalize.nf
process NORMALIZE {
    tag "normalize"
    publishDir "${params.outdir}/03_normalize", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_normalized.h5ad", emit: h5ad
    path "*.png"

    script:
    """
    python ${projectDir}/scripts/03_normalize.py \
        --input ${h5ad} \
        --out . \
        --n_hvgs ${params.n_hvgs}
    """
}
```

```groovy
// pipeline/modules/integrate.nf
process INTEGRATE {
    tag "integrate"
    publishDir "${params.outdir}/04_integrate", mode: 'copy'
    memory '64 GB'
    cpus 8
    accelerator 1  // request 1 GPU on SLURM

    input:
    path h5ad

    output:
    path "adata_integrated.h5ad", emit: h5ad
    path "scvi_model/", emit: model
    path "*.png"

    script:
    """
    python ${projectDir}/scripts/04_integrate.py \
        --input ${h5ad} \
        --out . \
        --n_latent ${params.n_latent} \
        --max_epochs ${params.max_epochs}
    """
}
```

```groovy
// pipeline/modules/cluster.nf
process CLUSTER {
    tag "cluster"
    publishDir "${params.outdir}/05_cluster", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "adata_clustered.h5ad", emit: h5ad
    path "*.png"

    script:
    """
    python ${projectDir}/scripts/05_cluster.py \
        --input ${h5ad} \
        --out . \
        --resolution ${params.leiden_resolution}
    """
}
```

```groovy
// pipeline/modules/annotate.nf
process ANNOTATE {
    tag "annotate"
    publishDir "${params.outdir}/06_annotate", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad
    path ref_h5ad

    output:
    path "adata_annotated.h5ad", emit: h5ad
    path "*.png"
    path "*.csv"

    script:
    """
    python ${projectDir}/scripts/06_annotate.py \
        --input ${h5ad} \
        --ref ${ref_h5ad} \
        --out .
    """
}
```

```groovy
// pipeline/modules/markers.nf
process MARKERS {
    tag "markers"
    publishDir "${params.outdir}/07_markers", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    path h5ad

    output:
    path "*.png"
    path "*.csv"

    script:
    """
    python ${projectDir}/scripts/07_markers.py \
        --input ${h5ad} \
        --out . \
        --cluster_key leiden_r${params.leiden_resolution}
    """
}
```

- [ ] **Step 8.2: Create main.nf**

```groovy
// pipeline/main.nf
nextflow.enable.dsl = 2

include { QC }        from './modules/qc'
include { DOUBLETS }  from './modules/doublets'
include { NORMALIZE } from './modules/normalize'
include { INTEGRATE } from './modules/integrate'
include { CLUSTER }   from './modules/cluster'
include { ANNOTATE }  from './modules/annotate'
include { MARKERS }   from './modules/markers'

workflow {
    h5_ch   = Channel.fromPath(params.h5_input)
    ref_ch  = Channel.fromPath(params.ref_h5ad)

    QC(h5_ch)
    DOUBLETS(QC.out.h5ad)
    NORMALIZE(DOUBLETS.out.h5ad)
    INTEGRATE(NORMALIZE.out.h5ad)
    CLUSTER(INTEGRATE.out.h5ad)
    ANNOTATE(CLUSTER.out.h5ad, ref_ch)
    MARKERS(ANNOTATE.out.h5ad)
}
```

- [ ] **Step 8.3: Create nextflow.config**

```groovy
// pipeline/nextflow.config
params {
    h5_input          = "${projectDir}/../data/count/filtered_feature_bc_matrix.h5"
    ref_h5ad          = "${projectDir}/../data/paul15_reference.h5ad"
    outdir            = "${projectDir}/../results"
    mito_nmads        = 4.0
    count_nmads       = 6.0
    n_hvgs            = 3000
    n_latent          = 30
    max_epochs        = 400
    leiden_resolution = 0.5
}

profiles {
    local {
        includeConfig 'conf/local.config'
    }
    slurm {
        includeConfig 'conf/slurm.config'
    }
}

process.shell = ['/bin/bash', '-euo', 'pipefail']
```

- [ ] **Step 8.4: Create local and SLURM configs**

```groovy
// pipeline/conf/local.config
process {
    executor = 'local'
    withName: 'INTEGRATE' {
        // local GPU run
        conda = null
    }
}
executor {
    cpus   = 8
    memory = '64 GB'
}
```

```groovy
// pipeline/conf/slurm.config
process {
    executor    = 'slurm'
    queue       = 'gpu'
    clusterOptions = '--partition=gpu --gres=gpu:1'

    withName: 'QC|DOUBLETS|NORMALIZE|CLUSTER|ANNOTATE|MARKERS' {
        queue          = 'general'
        clusterOptions = '--partition=general'
    }
    withName: 'INTEGRATE' {
        queue          = 'gpu'
        clusterOptions = '--partition=gpu --gres=gpu:1 --time=04:00:00'
    }
}
```

- [ ] **Step 8.5: Create params.yaml for alternative runs**

```yaml
# pipeline/params.yaml
h5_input: "../data/count/filtered_feature_bc_matrix.h5"
ref_h5ad: "../data/paul15_reference.h5ad"
outdir: "../results"
mito_nmads: 4.0
count_nmads: 6.0
n_hvgs: 3000
n_latent: 30
max_epochs: 400
leiden_resolution: 0.5
```

- [ ] **Step 8.6: Test local run**

```bash
cd /home/edu-pc/Yale/Charles/pipeline
nextflow run main.nf -profile local -resume
```

If Nextflow is not installed:
```bash
curl -s https://get.nextflow.io | bash
mv nextflow /home/edu-pc/.local/bin/
```

- [ ] **Step 8.7: Commit pipeline**

```bash
git add pipeline/
git commit -m "feat: Nextflow DSL2 pipeline wrapping all 7 analysis scripts"
```

---

## Task 9: GitHub Setup

**Files:**
- Create: `.gitignore`
- Create: `README.md`

- [ ] **Step 9.1: Initialize git repo**

```bash
cd /home/edu-pc/Yale/Charles
git init
```

- [ ] **Step 9.2: Create .gitignore**

```gitignore
# .gitignore
results/
data/count/
*.h5ad
*.h5
*.pyc
__pycache__/
.pixi/
work/           # Nextflow work directory
.nextflow/
.nextflow.log*
*.cloupe
```

- [ ] **Step 9.3: Create README.md**

```markdown
# Charles scRNA-seq Pipeline

Single-cell RNA-seq processing pipeline for mouse hematopoiesis (DMSO vs STM treatment).

## Data
- 24 samples: 2 donors × 2 treatments (DMSO/STM) × 3 populations (LSK/LK/Input) × 2 replicates
- Mouse (mm10), 10X Genomics CellPlex multiplexing
- 53,494 cells, 32,285 genes

## Pipeline
QC → Doublet removal → Normalization → scVI integration → Leiden clustering → Cell type annotation → Marker validation

## Run

```bash
# Local
cd pipeline
nextflow run main.nf -profile local -resume

# Yale HPC (Grace)
nextflow run main.nf -profile slurm -resume
```

## Environment
Uses pixi (see `~/Yale/pixi.toml`).
```

- [ ] **Step 9.4: Create GitHub repo and push**

```bash
cd /home/edu-pc/Yale/Charles
git add pipeline/ docs/ README.md .gitignore
git commit -m "feat: initial pipeline commit — scRNA-seq QC to annotation"
gh repo create charles-scrna-pipeline --private --push --source .
```

(Requires `gh` CLI authenticated with `gh auth login`)

---

## Self-Review Checklist

- [x] QC: MAD-based per-sample filtering, HTO separated to obsm
- [x] Doublets: DoubletFinder per-sample through Python-driven R bridge (not full matrix)
- [x] Normalization: raw counts stored before normalization (required by scVI)
- [x] HVG: `batch_key=sample_id` for cross-sample HVG selection
- [x] Integration: scVI with `categorical_covariate_keys` preserving biology
- [x] Clustering: Leiden at multiple resolutions
- [ ] Annotation: popV primary with Nestorowa reference + scANVI/scArches secondary
- [x] Markers: dotplot, violin, UMAP expression, Wilcoxon DE
- [x] Nextflow: all 7 scripts wrapped, local + SLURM profiles
- [x] GitHub: .gitignore excludes data/results, README explains experiment

**Known gaps to discuss with PI:**
- D1_DMSO_LSK_45_2 (134 cells) and D1_STM_LSK_45_2 (290 cells) — decide whether to exclude before integration
- Ambient RNA: skipped because raw matrix unavailable; note in Methods
- Nestorowa reference still needs to be generated at `data/references/nestorowa_2016_hspc.h5ad`
