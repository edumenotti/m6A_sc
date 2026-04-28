# scRNA-seq Processing Pipeline — Charles Dataset Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible Nextflow pipeline that takes the aggregated CellRanger H5 matrix from the Charles dataset (24 samples, mouse hematopoiesis, DMSO vs STM, LSK/LK/Input populations) and produces QC-filtered, batch-integrated, clustered, and annotated AnnData objects with cell type labels validated against canonical markers.

**Architecture:** Six sequential Python scripts (one per analysis stage) are orchestrated by a Nextflow DSL2 pipeline. Each script reads a `.h5ad` from the previous step and writes the next. The pipeline can run locally (pixi env) or on a SLURM cluster (Yale HPC/Grace).

**Tech Stack:** Python 3.14 · scanpy 1.12 · scvi-tools 1.4 · scArches · popV · CellTypist · Nextflow DSL2 · pixi (existing env at ~/Yale)

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
CellTypist has limited mouse models. The `Immune_All_Low.pkl` model is human-trained but the tool provides a human→mouse gene conversion. It will be used as an **independent validation**, not primary annotation.

**Recommended strategy (this plan):**
1. **popV** with Nestorowa 2016 reference → primary labels + uncertainty scores
2. **scArches + KNN** in scVI latent space → secondary reference mapping
3. **CellTypist** `Immune_All_Low.pkl` with `over_clustering=True` → cross-species validation
4. **Manual marker dotplot** → ground truth validation

### Doublet Detection: Scrublet
Scrublet generates synthetic doublets and computes a doublet score per cell. SOLO (deep learning doublet detector in scvi-tools) is slightly more accurate but Scrublet is faster and sufficient. Run **per-sample** (per barcode suffix) to avoid inter-sample confounding.

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
│   ├── 01_qc/
│   ├── 02_doublets/
│   ├── 03_normalize/
│   ├── 04_integrate/
│   ├── 05_cluster/
│   ├── 06_annotate/
│   └── 07_markers/
└── docs/
    └── superpowers/plans/  (this file)
```

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
- Generates QC violin/scatter plots
- Filters cells using MAD-based thresholds (not fixed cutoffs)
- Saves filtered AnnData

### Why MAD-based thresholds?
Fixed thresholds (e.g., "remove cells with < 200 genes") are arbitrary and dataset-dependent. MAD (median absolute deviation) thresholds adapt to each sample's distribution. Standard practice: flag cells that are >5 MADs below the median for n_genes or n_counts, or >3 MADs above median for pct_mito.

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
    parser.add_argument("--mito_nmads", type=float, default=3.0)
    parser.add_argument("--count_nmads", type=float, default=5.0)
    args = parser.parse_args()
    main(args.h5, args.out, args.mito_nmads, args.count_nmads)
```

- [ ] **Step 1.2: Run to verify it works**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/01_qc.py \
  --h5 Charles/data/count/filtered_feature_bc_matrix.h5 \
  --out Charles/results/01_qc/
```

Expected output:
```
Cells before QC: 53494
Cells removed: ~500-2000
Cells after QC: ~51000-53000
Saved: Charles/results/01_qc/adata_qc.h5ad
```

- [ ] **Step 1.3: Commit**

```bash
git add pipeline/scripts/01_qc.py
git commit -m "feat: QC script with MAD-based per-sample filtering"
```

---

## Task 2: Doublet Detection (`scripts/02_doublets.py`)

**Files:**
- Create: `pipeline/scripts/02_doublets.py`
- Reads: `results/01_qc/adata_qc.h5ad`
- Writes: `results/02_doublets/adata_no_doublets.h5ad`

### Why per-sample doublet detection?
Scrublet simulates doublets by combining pairs of cells from the input. If run on the full aggregated matrix, it would simulate cross-sample doublets (which don't exist — each cell has a single barcode suffix). Running per-sample produces realistic synthetic doublets.

- [ ] **Step 2.1: Create the doublet script**

```python
# pipeline/scripts/02_doublets.py
import scanpy as sc
import numpy as np
import scrublet as scr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def run_scrublet(adata_sub: sc.AnnData, expected_doublet_rate: float = 0.06) -> np.ndarray:
    counts = adata_sub.X
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    scrub = scr.Scrublet(counts, expected_doublet_rate=expected_doublet_rate)
    scores, _ = scrub.scrub_doublets(min_counts=2, min_cells=3, n_prin_comps=30, verbose=False)
    return scores, scrub.threshold_

def main(in_path: str, out_dir: str, score_threshold: float = None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    scores = np.zeros(adata.n_obs)
    thresholds = {}

    for sid in adata.obs["sample_id"].unique():
        mask = adata.obs["sample_id"] == sid
        idx = np.where(mask)[0]
        sub = adata[mask].copy()
        s, t = run_scrublet(sub)
        scores[idx] = s
        thresholds[sid] = t
        print(f"  {sid}: {mask.sum()} cells, threshold={t:.3f}, "
              f"predicted doublets={(s > t).sum()}")

    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = False

    for sid, t in thresholds.items():
        mask = adata.obs["sample_id"] == sid
        thr = score_threshold if score_threshold else t
        adata.obs.loc[mask, "predicted_doublet"] = (
            adata.obs.loc[mask, "doublet_score"] > thr
        )

    # Plot score distributions
    fig, ax = plt.subplots(figsize=(12, 4))
    for sid in adata.obs["sample_id"].unique():
        mask = adata.obs["sample_id"] == sid
        ax.hist(adata.obs.loc[mask, "doublet_score"], bins=50,
                alpha=0.4, label=sid, density=True)
    ax.set_xlabel("Doublet score")
    ax.set_ylabel("Density")
    ax.legend(fontsize=6, ncol=4)
    plt.tight_layout()
    plt.savefig(out / "doublet_scores.png", dpi=150)
    plt.close()

    n_doublets = adata.obs["predicted_doublet"].sum()
    print(f"Total doublets removed: {n_doublets} ({n_doublets/adata.n_obs*100:.1f}%)")
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    adata.write_h5ad(out / "adata_no_doublets.h5ad")
    print(f"Saved: {out / 'adata_no_doublets.h5ad'}, {adata.n_obs} cells")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--score_threshold", type=float, default=None)
    args = parser.parse_args()
    main(args.input, args.out, args.score_threshold)
```

- [ ] **Step 2.2: Install scrublet in pixi env**

Add to `pixi.toml` under `[pypi-dependencies]`:
```toml
scrublet = ">=0.2.3"
```
Then run: `cd /home/edu-pc/Yale && pixi install`

- [ ] **Step 2.3: Run**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/02_doublets.py \
  --input Charles/results/01_qc/adata_qc.h5ad \
  --out Charles/results/02_doublets/
```

Expected: ~5-10% doublets removed per sample.

- [ ] **Step 2.4: Commit**

```bash
git add pipeline/scripts/02_doublets.py pixi.toml pixi.lock
git commit -m "feat: per-sample doublet detection with Scrublet"
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
- Exclude mitochondrial genes from HVG list (they are QC metrics, not biology)

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

    # Exclude mitochondrial genes from HVGs
    adata.var.loc[adata.var_names.str.startswith("mt-"), "highly_variable"] = False

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

- [ ] **Step 4.1: Create the integration script**

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

- [ ] **Step 4.2: Run (GPU recommended, ~15-30 min)**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/04_integrate.py \
  --input Charles/results/03_normalize/adata_normalized.h5ad \
  --out Charles/results/04_integrate/ \
  --n_latent 30 --max_epochs 400
```

Expected: ELBO converging, UMAP without obvious pool-based separation.

- [ ] **Step 4.3: Check integration quality**

Inspect `results/04_integrate/umap_pool.png`. If cells cluster strongly by pool (not biology), increase `max_epochs` or adjust `n_latent`.

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

- [ ] **Step 5.1: Create the clustering script**

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

- [ ] **Step 5.2: Run**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/05_cluster.py \
  --input Charles/results/04_integrate/adata_integrated.h5ad \
  --out Charles/results/05_cluster/ \
  --resolution 0.5
```

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
3. Run CellTypist with `Immune_All_Low.pkl` as validation

**About the reference subset question (practical explanation):**

The Nestorowa 2016 reference contains ONLY mouse HSPCs (HSC, MPP, LMPP, CMP, GMP, MEP, CLP, MkP, ErP). There are no other cell types. So for our LSK/LK/Input populations this reference is ideal — no subsetting needed for the HSPC populations.

For the Input (whole BM) fraction, there will be mature immune cells (T cells, B cells, NK cells, neutrophils) that are not in Nestorowa. For these, popV will return an "unknown" or the closest match with low agreement score. These cells can be annotated with CellTypist instead.

**Practical subsetting if you need it (for other datasets):**
```python
# Example: subset the meta-analytic atlas to only HSPC types
hspc_types = ["HSC", "MPP", "LMPP", "CMP", "GMP", "MEP", "CLP", "MkP", "ErP", "DC_progenitor"]
ref_sub = ref[ref.obs["cell_type"].isin(hspc_types)].copy()
# Then use ref_sub as the reference for popV/scArches
```

- [ ] **Step 6.1: Prepare Nestorowa reference (run once)**

```python
# This block is run once to download and save the reference
# Run interactively or as a separate setup script

import scrnaseq  # pip install scrnaseq — Python wrapper for Bioconductor
# Alternative: download directly from GEO GSE81682

# Method 1: via scRNAseq R package (recommended) — run in R:
# library(scRNAseq)
# ref <- NestorowaHSCData()
# library(zellkonverter)
# writeH5AD(ref, "nestorowa_2016_ref.h5ad")

# Method 2: use scanpy built-in (HSPC data from Nestorowa)
import scanpy as sc
# The data is not directly in sc.datasets, but Paul15 is:
paul = sc.datasets.paul15()
# paul15 has 2730 mouse HSPCs with 19 cell type labels
# Available labels: HSC, MEP, Ery, CMP, GMP, DC, Mono, Mast, Lympho, Baso, Neu...
paul.write_h5ad("Charles/data/paul15_reference.h5ad")
print(paul)
print(paul.obs["paul15_clusters"].value_counts())
```

- [ ] **Step 6.2: Create the annotation script**

```python
# pipeline/scripts/06_annotate.py
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import celltypist
from celltypist import models
import argparse
from pathlib import Path

def annotate_with_celltypist(adata: ad.AnnData, out: Path):
    """CellTypist annotation as validation layer."""
    # CellTypist requires log-normalized counts (already in lognorm layer)
    adata_ct = adata.copy()
    adata_ct.X = adata_ct.layers["lognorm"]

    # Download model if not present
    models.download_models(model="Immune_All_Low.pkl", force_update=False)
    model = models.Model.load("Immune_All_Low.pkl")

    # over_clustering=True uses leiden clusters for majority-vote smoothing
    predictions = celltypist.annotate(
        adata_ct,
        model="Immune_All_Low.pkl",
        majority_voting=True,
        over_clustering="leiden_r0.5",
    )
    adata.obs["celltypist_label"] = predictions.predicted_labels["majority_voting"].values
    adata.obs["celltypist_conf"] = predictions.predicted_labels["conf_score"].values

    predictions.to_plots(out / "celltypist", show=False)
    return adata

def annotate_with_paul15(adata: ad.AnnData, ref_path: str, out: Path):
    """KNN-based label transfer from Paul15 reference in scVI latent space."""
    from sklearn.neighbors import KNeighborsClassifier

    ref = sc.read_h5ad(ref_path)

    # Compute log-normalized PCA on reference
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)

    # Find common genes
    common_genes = adata.var_names.intersection(ref.var_names)
    print(f"Common genes with Paul15 reference: {len(common_genes)}")

    ref_sub = ref[:, common_genes].copy()
    query_sub = adata[:, common_genes].copy()

    sc.pp.highly_variable_genes(ref_sub, n_top_genes=2000, flavor="seurat_v3")
    sc.tl.pca(ref_sub, use_highly_variable=True, n_comps=30)

    # Project query onto reference PCA
    from sklearn.decomposition import PCA
    hvg_mask = ref_sub.var["highly_variable"]
    pca = PCA(n_components=30).fit(ref_sub.X[:, hvg_mask].toarray()
                                   if hasattr(ref_sub.X, "toarray")
                                   else ref_sub.X[:, hvg_mask])
    query_pca = pca.transform(
        query_sub[:, common_genes[hvg_mask]].layers["lognorm"].toarray()
        if hasattr(query_sub.layers["lognorm"], "toarray")
        else query_sub[:, common_genes[hvg_mask]].layers["lognorm"]
    )

    knn = KNeighborsClassifier(n_neighbors=15, metric="euclidean")
    knn.fit(ref_sub.obsm["X_pca"], ref.obs["paul15_clusters"].values)
    labels = knn.predict(query_pca)
    proba = knn.predict_proba(query_pca).max(axis=1)

    adata.obs["paul15_label"] = labels
    adata.obs["paul15_conf"] = proba
    return adata

def main(in_path: str, ref_path: str, out_dir: str):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(in_path)

    print("Running CellTypist annotation...")
    adata = annotate_with_celltypist(adata, out)

    print("Running Paul15 KNN label transfer...")
    adata = annotate_with_paul15(adata, ref_path, out)

    # Consensus: for each cell, report both labels + agreement
    adata.obs["labels_agree"] = (
        adata.obs["celltypist_label"] == adata.obs["paul15_label"]
    )

    # Summary per cluster
    for cluster_key in ["leiden_r0.5"]:
        summary = adata.obs.groupby(cluster_key).agg(
            paul15_top=("paul15_label", lambda x: x.value_counts().index[0]),
            ct_top=("celltypist_label", lambda x: x.value_counts().index[0]),
            n_cells=("paul15_label", "count"),
        )
        summary.to_csv(out / f"annotation_summary_{cluster_key}.csv")
        print(summary)

    # UMAP with annotations
    for col in ["paul15_label", "celltypist_label"]:
        sc.pl.umap(adata, color=col, show=False, legend_loc="on data")
        plt.savefig(out / f"umap_{col}.png", dpi=150, bbox_inches="tight")
        plt.close()

    adata.write_h5ad(out / "adata_annotated.h5ad")
    print(f"Saved: {out / 'adata_annotated.h5ad'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--ref", required=True, help="Path to paul15_reference.h5ad")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    main(args.input, args.ref, args.out)
```

- [ ] **Step 6.3: Download Paul15 reference**

```bash
cd /home/edu-pc/Yale
pixi run python -c "
import scanpy as sc
paul = sc.datasets.paul15()
paul.write_h5ad('Charles/data/paul15_reference.h5ad')
print(paul.obs['paul15_clusters'].value_counts())
"
```

- [ ] **Step 6.4: Run annotation**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/06_annotate.py \
  --input Charles/results/05_cluster/adata_clustered.h5ad \
  --ref Charles/data/paul15_reference.h5ad \
  --out Charles/results/06_annotate/
```

- [ ] **Step 6.5: Commit**

```bash
git add pipeline/scripts/06_annotate.py
git commit -m "feat: dual annotation with CellTypist + Paul15 KNN label transfer"
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

- [ ] **Step 7.1: Create the marker validation script**

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

- [ ] **Step 7.2: Run**

```bash
cd /home/edu-pc/Yale
pixi run python pipeline/scripts/07_markers.py \
  --input Charles/results/06_annotate/adata_annotated.h5ad \
  --out Charles/results/07_markers/
```

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
    mito_nmads        = 3.0
    count_nmads       = 5.0
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
mito_nmads: 3.0
count_nmads: 5.0
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
- [x] Doublets: Scrublet per-sample (not full matrix)
- [x] Normalization: raw counts stored before normalization (required by scVI)
- [x] HVG: `batch_key=sample_id` for cross-sample HVG selection
- [x] Integration: scVI with `categorical_covariate_keys` preserving biology
- [x] Clustering: Leiden at multiple resolutions
- [x] Annotation: dual CellTypist + Paul15 KNN, agreement column
- [x] Markers: dotplot, violin, UMAP expression, Wilcoxon DE
- [x] Nextflow: all 7 scripts wrapped, local + SLURM profiles
- [x] GitHub: .gitignore excludes data/results, README explains experiment

**Known gaps to discuss with PI:**
- D1_DMSO_LSK_45_2 (134 cells) and D1_STM_LSK_45_2 (290 cells) — decide whether to exclude before integration
- Ambient RNA: skipped because raw matrix unavailable; note in Methods
- popV not yet installed in pixi.toml (add `popv` to pypi-dependencies when available; as of 2024 it requires a separate conda env due to dependency conflicts)
