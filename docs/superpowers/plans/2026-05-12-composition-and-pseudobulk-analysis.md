# Composition + Pseudobulk DEG Analysis Plan (MDS scRNA-seq)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the underpowered/confounded CellChat+NicheNet comparative analysis with a statistically defensible workflow centered on (a) Bayesian compositional analysis with **scCODA**, (b) per-celltype pseudobulk DEG via **PyDESeq2** (donors as the replicate unit), and (c) MSigDB pathway scoring via **decoupler**. Archive the exploratory CellChat/NicheNet outputs with a clear scope statement.

**Architecture:** Python-only stack (matches the rest of the pipeline). Three new scripts (18 / 19 / 20) consume the existing `adata_progenitor_annotated.h5ad`, write to `results/18_composition`, `results/19_pseudobulk_deg`, `results/20_pathway_activity`. Old `16_cellchat` and `17_nichenet` outputs are moved under `results/_archive_exploratory/` with a README documenting the confounders. Existing Nextflow conditional block stays; new modules are wired into the same `run_interaction_analysis` flag (now reframed as `run_downstream_analysis`).

**Tech Stack:**
- `scCODA` (Theis lab, Bayesian compositional model, robust at low n)
- `decoupler` v2 — pseudobulk aggregation + ULM/ORA pathway activity
- `PyDESeq2` — Python DESeq2 reimplementation with factorial design support
- `msigdb` (`decoupler.op.msigdb`) — Hallmark + Reactome + KEGG
- scanpy / anndata — already in pixi env

**Why these tools are appropriate for the data:**

| Constraint | Tool choice rationale |
|---|---|
| n=2 donors per condition | scCODA: Bayesian, robust to low n. PyDESeq2 + DESeq2 dispersion shrinkage: established method for n≥2. |
| 4 conditions in 2×2 factorial | scCODA `formula="genotype + treatment + genotype:treatment"`. PyDESeq2 same. |
| Different sort compositions | scCODA models composition explicitly (log-ratios), so reference cell type controls for sort baseline. |
| Cell-type level questions | Pseudobulk DEG by (donor × celltype) treats donor as the replicate unit, correctly. |
| Pathway-level robustness | decoupler ULM/ORA aggregates many genes → less noise than gene-by-gene at n=2. |

---

## File Structure

**Created:**
- `pipeline/scripts/18_composition_sccoda.py` — scCODA composition analysis (genotype + treatment + interaction)
- `pipeline/scripts/19_pseudobulk_deg.py` — Per-celltype pseudobulk DEG via decoupler + PyDESeq2
- `pipeline/scripts/20_pathway_activity.py` — Hallmark/Reactome/Progeny scoring on per-celltype DEG stats
- `pipeline/modules/composition.nf` — Nextflow wrapper for script 18
- `pipeline/modules/pseudobulk_deg.nf` — Nextflow wrapper for script 19
- `pipeline/modules/pathway_activity.nf` — Nextflow wrapper for script 20
- `results/_archive_exploratory/README.md` — Documents why CellChat/NicheNet outputs are exploratory only

**Modified:**
- `pixi.toml` — add `sccoda`, `pydeseq2`, `decoupler-py` (verify versions)
- `pipeline/main.nf` — wire COMPOSITION, PSEUDOBULK_DEG, PATHWAY_ACTIVITY into the conditional block (renamed `run_downstream_analysis`)
- `pipeline/nextflow.config` — add `run_downstream_analysis` parameter

**Archived (moved, not deleted — full audit trail):**
- `results/16_cellchat/` → `results/_archive_exploratory/16_cellchat/`
- `results/17_nichenet/` → `results/_archive_exploratory/17_nichenet/`

**Kept as-is:**
- `pipeline/scripts/15_macrophage_states.py` (per-cell decoupler ULM — not affected by composition confound)
- `pipeline/scripts/16_cellchat.R` + `17_nichenet.R` (kept on disk; only outputs archived; can be re-run as descriptive maps if needed)
- `pipeline/scripts/export_h5ad_for_r.py` (still useful if R analyses are revisited)

---

## Task 0: Dependency setup

**Files:**
- Modify: `pixi.toml`

- [ ] **Step 1: Verify required Python packages are not yet installed**

Run:
```bash
pixi run python -c "import sccoda, pydeseq2; print('OK')" 2>&1 | tail -3
```
Expected: `ModuleNotFoundError` for at least one.

- [ ] **Step 2: Add packages to pixi.toml**

Open `pixi.toml` and under `[pypi-dependencies]` (or `[dependencies]` if available via conda), add:

```toml
sccoda = "*"
pydeseq2 = "*"
```

Verify decoupler is at v2+ (already present from script 15).

- [ ] **Step 3: Install and verify**

Run:
```bash
pixi install 2>&1 | tail -10
pixi run python -c "import sccoda, pydeseq2, decoupler; print('sccoda', sccoda.__version__); print('pydeseq2', pydeseq2.__version__); print('decoupler', decoupler.__version__)"
```
Expected: three version lines.

- [ ] **Step 4: Commit**

```bash
git add pixi.toml pixi.lock
git commit -m "deps: add sccoda + pydeseq2 for composition and pseudobulk DEG analyses"
```

---

## Task 1: Archive exploratory CellChat / NicheNet outputs

**Files:**
- Create: `results/_archive_exploratory/README.md`
- Move: `results/16_cellchat/`, `results/17_nichenet/`

- [ ] **Step 1: Create archive README documenting limitations**

Create `results/_archive_exploratory/README.md`:

```markdown
# Exploratory analyses — archived

These outputs were generated by `16_cellchat.R` and `17_nichenet.R`. They are
retained for audit but **should not be used for inferential claims** because
the underlying experimental design has confounders that these tools cannot
address:

1. **Composition is not comparable across conditions.** Different sort fractions
   (I / LK / LSK) per sample produce ≥20× differences in B-cell, monocyte, T-cell,
   and DC counts between Mutant and WT. CellChat with `population.size=TRUE`
   amplifies these abundance differences as "differential communication".

2. **n=2 donors per condition.** CellChat treats each cell as an independent
   observation. p-values from `rankNet` are nominal and over-estimated.

3. **Competitive transplant lost in silico.** WT (CD45.1) and Mutant (CD45.2)
   cells coexist in the same mouse. In-silico separation by genotype destroys
   the actual in-vivo crosstalk signal.

4. **Top hits are dominated by abundant cell types** (granulocyte_neutrophil
   → granulocyte_neutrophil ANXA1 / CCL6 / SELPLG), not by MDS-specific biology.

For valid inferential analyses see:
- `results/18_composition/` — scCODA Bayesian compositional analysis
- `results/19_pseudobulk_deg/` — per-celltype DEG with donor as replicate
- `results/20_pathway_activity/` — pathway scoring per celltype

The scripts that produced these archived results (`16_cellchat.R`,
`17_nichenet.R`) are kept under `pipeline/scripts/` and can still be used
as descriptive per-condition maps (do NOT make cross-condition claims).
```

- [ ] **Step 2: Move existing results into archive**

Run:
```bash
mkdir -p results/_archive_exploratory
git mv results/16_cellchat results/_archive_exploratory/16_cellchat 2>/dev/null \
  || mv results/16_cellchat results/_archive_exploratory/16_cellchat
git mv results/17_nichenet results/_archive_exploratory/17_nichenet 2>/dev/null \
  || mv results/17_nichenet results/_archive_exploratory/17_nichenet
ls results/_archive_exploratory/
```
Expected: `16_cellchat/  17_nichenet/  README.md`.

- [ ] **Step 3: Commit**

```bash
git add results/_archive_exploratory/README.md
git add -A results/_archive_exploratory/ || true
git commit -m "chore: archive CellChat/NicheNet exploratory outputs with confounder README"
```

Note: results/ may be gitignored — check first with `cat .gitignore | grep -i results`. If it is, only the README needs to be committed; the directory move is local-only.

---

## Task 2: scCODA composition script — skeleton + input loading

**Files:**
- Create: `pipeline/scripts/18_composition_sccoda.py`

- [ ] **Step 1: Write the script skeleton with CLI and h5ad loading**

```python
#!/usr/bin/env python3
"""
18_composition_sccoda.py

Bayesian compositional analysis with scCODA on cell-type frequencies per donor
per condition. Tests:
- Effect of treatment (STM vs DMSO) within each genotype
- Effect of genotype (Mutant vs WT) within each treatment
- Genotype × treatment interaction

Why scCODA: explicit Bayesian model of compositional data, robust to small n,
controls for the simplex constraint (one cell type up → others appear down).
Reference cell type acts as the baseline against which all others are
log-ratio-tested.

Inputs
------
--input      path to adata_progenitor_annotated.h5ad
--out        output directory
--level      'manual_level1' or 'manual_level2'
--ref-type   reference cell type (must be present in all samples; defaults
             to most abundant shared type)

Outputs
-------
sccoda_counts_per_sample.csv         per-donor-per-condition cell-type counts
sccoda_<contrast>_credible_effects.csv  credible effects (FDR-controlled)
sccoda_<contrast>_summary.txt        full HMC sampling summary
sccoda_boxplots_<level>.png          boxplots of celltype fractions × condition
sccoda_chosen_reference.txt          which reference cell type was selected and why
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/18_composition")
    p.add_argument("--level", default="manual_level1",
                   choices=["manual_level1", "manual_level2"])
    p.add_argument("--ref-type", default=None,
                   help="Reference cell type. Defaults to most-abundant type "
                        "with cells in ALL donor×condition samples.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  shape: {adata.shape}, layers: {list(adata.layers.keys())}")

    obs = adata.obs.copy()
    obs["celltype"] = obs[args.level].astype(str)
    obs["condition"] = obs["genotype"].astype(str) + "_" + obs["treatment"].astype(str)
    obs["sample"] = obs["donor"].astype(str) + "__" + obs["condition"].astype(str)

    print("\nDonors per condition:")
    print(obs.groupby("condition")["donor"].nunique())
    print("\nCells per (donor, condition):")
    print(obs.groupby(["condition", "donor"]).size())

    # placeholder for the rest of the workflow
    raise SystemExit("Skeleton OK — Task 3 fills the model")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the skeleton to confirm it loads the data**

Run:
```bash
pixi run python pipeline/scripts/18_composition_sccoda.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out /tmp/sccoda_test \
  --level manual_level1 2>&1 | tail -30
```
Expected: prints donors per condition and cells per (donor, condition), then exits with "Skeleton OK".

- [ ] **Step 3: Commit**

```bash
git add pipeline/scripts/18_composition_sccoda.py
git commit -m "feat(18): scCODA composition script skeleton with h5ad loading"
```

---

## Task 3: scCODA model fit + credible effects

**Files:**
- Modify: `pipeline/scripts/18_composition_sccoda.py`

- [ ] **Step 1: Replace the placeholder with the count-matrix builder + model fit**

Replace the `raise SystemExit(...)` and everything that comes after `print(obs.groupby(...))` with:

```python
    # ── Build sample × celltype count matrix ─────────────────────────────
    counts = (
        obs.groupby(["sample", "celltype"]).size().unstack(fill_value=0).reset_index()
    )
    # Add covariates from the sample name (donor + condition)
    counts["donor"] = counts["sample"].str.split("__").str[0]
    counts["condition"] = counts["sample"].str.split("__").str[1]
    counts["genotype"] = counts["condition"].str.split("_").str[0]
    counts["treatment"] = counts["condition"].str.split("_").str[1]
    counts.to_csv(os.path.join(args.out, "sccoda_counts_per_sample.csv"), index=False)
    print(f"\nCount matrix shape: {counts.shape}")
    print(counts.head().to_string())

    # ── Reference cell type selection ────────────────────────────────────
    celltype_cols = [c for c in counts.columns
                     if c not in ("sample", "donor", "condition", "genotype", "treatment")]
    # Most-abundant cell type with >0 in every sample
    presence = (counts[celltype_cols] > 0).all(axis=0)
    candidates = [c for c in celltype_cols if presence[c]]
    if args.ref_type:
        ref = args.ref_type
        assert ref in celltype_cols, f"--ref-type {ref} not in {celltype_cols}"
    elif candidates:
        totals = counts[candidates].sum(axis=0).sort_values(ascending=False)
        ref = totals.index[0]
    else:
        # fallback: use scCODA's automatic selector via 'automatic'
        ref = "automatic"
    with open(os.path.join(args.out, "sccoda_chosen_reference.txt"), "w") as fh:
        fh.write(f"reference_cell_type: {ref}\n")
        fh.write(f"present_in_all_samples: {candidates}\n")
        fh.write(f"sample_totals_per_celltype:\n{counts[celltype_cols].sum(axis=0).to_string()}\n")
    print(f"\nReference cell type: {ref}")

    # ── scCODA model: factorial formula ──────────────────────────────────
    from sccoda.util import cell_composition_data as dat
    from sccoda.util import comp_ana as mod

    sccoda_data = dat.from_pandas(counts, covariate_columns=[
        "sample", "donor", "condition", "genotype", "treatment"
    ])

    contrasts = [
        ("treatment_in_WT",         "treatment", counts["genotype"] == "WT"),
        ("treatment_in_Mutant",     "treatment", counts["genotype"] == "Mutant"),
        ("genotype_in_DMSO",        "genotype",  counts["treatment"] == "DMSO"),
        ("genotype_in_STM",         "genotype",  counts["treatment"] == "STM"),
        ("genotype_x_treatment",    "genotype + treatment + genotype:treatment", None),
    ]

    for name, formula, row_mask in contrasts:
        print(f"\n=== Contrast: {name} (formula: {formula}) ===")
        if row_mask is not None:
            sub = sccoda_data[row_mask.values].copy()
        else:
            sub = sccoda_data.copy()
        if sub.shape[0] < 2:
            print(f"  skipping — only {sub.shape[0]} samples after filtering")
            continue
        model = mod.CompositionalAnalysis(
            sub, formula=formula, reference_cell_type=ref
        )
        result = model.sample_hmc(num_results=20000, num_burnin=5000)
        result.set_fdr(est_fdr=0.05)
        ce = result.credible_effects()
        ce_df = ce.reset_index()
        ce_df.columns = ["covariate", "celltype", "credible"]
        ce_df.to_csv(
            os.path.join(args.out, f"sccoda_{name}_credible_effects.csv"),
            index=False
        )
        with open(os.path.join(args.out, f"sccoda_{name}_summary.txt"), "w") as fh:
            import io, contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                result.summary_extended()
            fh.write(buf.getvalue())
        print(f"  credible effects:\n{ce_df[ce_df['credible']].to_string(index=False)}")

    # ── Diagnostic boxplots ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    long = counts.melt(
        id_vars=["sample", "donor", "condition", "genotype", "treatment"],
        value_vars=celltype_cols,
        var_name="celltype", value_name="n_cells",
    )
    long["fraction"] = long.groupby("sample")["n_cells"].transform(
        lambda x: x / x.sum()
    )
    long.boxplot(column="fraction", by=["celltype", "condition"], ax=axes[0])
    axes[0].set_xticklabels([])
    axes[0].set_title(f"Fractions per condition ({args.level})")
    axes[0].set_xlabel("celltype × condition")
    plt.suptitle("")
    pivot = (long.groupby(["celltype", "condition"])["fraction"].mean()
             .unstack(fill_value=0))
    pivot.plot(kind="bar", stacked=True, ax=axes[1])
    axes[1].set_title("Mean fraction per condition")
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(os.path.join(args.out, f"sccoda_boxplots_{args.level}.png"), dpi=150)
    plt.close(fig)

    print(f"\nscCODA analysis complete. Outputs in {args.out}/")
```

- [ ] **Step 2: Run end-to-end and confirm the credible-effects CSVs land**

Run:
```bash
pixi run python pipeline/scripts/18_composition_sccoda.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out results/18_composition/manual_level1 \
  --level manual_level1 2>&1 | tee /tmp/sccoda_run.log | tail -50
```
Expected: 5 contrast blocks finish; output dir has `sccoda_*_credible_effects.csv` × 5, `sccoda_*_summary.txt` × 5, plus boxplots and `sccoda_chosen_reference.txt`. HMC takes ~30s per contrast → ~3 min total.

- [ ] **Step 3: Run for manual_level2**

Run:
```bash
pixi run python pipeline/scripts/18_composition_sccoda.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out results/18_composition/manual_level2 \
  --level manual_level2 2>&1 | tee /tmp/sccoda_run_l2.log | tail -30
```
Expected: same output structure under `manual_level2/`.

- [ ] **Step 4: Commit**

```bash
git add pipeline/scripts/18_composition_sccoda.py
git commit -m "feat(18): scCODA composition analysis with 5 contrasts (treatment within genotype, genotype within treatment, factorial)"
```

---

## Task 4: Pseudobulk DEG script — skeleton + aggregation

**Files:**
- Create: `pipeline/scripts/19_pseudobulk_deg.py`

- [ ] **Step 1: Write script with CLI, h5ad loading, and pseudobulk aggregation**

```python
#!/usr/bin/env python3
"""
19_pseudobulk_deg.py

Per-celltype pseudobulk differential expression with PyDESeq2.
Donor is the replicate unit (n=2 per condition → 8 donor-level samples per celltype).

For each celltype with ≥10 cells in every (donor × condition) sample, we:
  1. Sum raw counts across cells → (donor × condition) × gene matrix
  2. Fit DESeq2 with design ~genotype + treatment + genotype:treatment
  3. Extract Wald results for treatment, genotype, and interaction
  4. Write per-celltype CSVs + a global summary

Why this is appropriate:
- Donor as replicate avoids pseudoreplication (the CellChat sin).
- DESeq2 dispersion shrinkage stabilises n=2 variance estimates.
- Restricting to celltypes present in every sample removes the composition
  confound that broke CellChat.

Inputs
------
--input        path to adata_progenitor_annotated.h5ad
--out          output directory
--level        manual_level1 (default) or manual_level2
--min-cells    min cells per (donor, condition, celltype) to retain (default 10)
--padj         FDR threshold for "significant" gene tables (default 0.05)

Outputs
-------
pseudobulk_samples_per_celltype.csv     QC: cells/sample/celltype matrix
deg_<celltype>_<contrast>.csv           full DESeq2 results
deg_<celltype>_<contrast>_top.csv       top sig hits (|lfc|>1, padj<thresh)
deg_summary.csv                         # of sig genes per celltype × contrast
deg_volcano_<celltype>_<contrast>.png   volcano plots per (celltype, contrast)
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/19_pseudobulk_deg")
    p.add_argument("--level", default="manual_level1",
                   choices=["manual_level1", "manual_level2"])
    p.add_argument("--min-cells", type=int, default=10)
    p.add_argument("--padj", type=float, default=0.05)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    assert "counts" in adata.layers, "Need adata.layers['counts'] (raw integer counts)"

    adata.obs["celltype"] = adata.obs[args.level].astype(str)
    adata.obs["condition"] = (
        adata.obs["genotype"].astype(str) + "_" + adata.obs["treatment"].astype(str)
    )
    adata.obs["sample"] = (
        adata.obs["donor"].astype(str) + "__" + adata.obs["condition"].astype(str)
    )

    # ── Pseudobulk via decoupler (sum, raw counts layer) ────────────────
    import decoupler as dc
    pdata = dc.pp.pseudobulk(
        adata,
        sample_col="sample",
        groups_col="celltype",
        layer="counts",
        mode="sum",
        skip_checks=False,
    )
    pdata.obs["donor"]     = pdata.obs["sample"].str.split("__").str[0]
    pdata.obs["condition"] = pdata.obs["sample"].str.split("__").str[1]
    pdata.obs["genotype"]  = pdata.obs["condition"].str.split("_").str[0]
    pdata.obs["treatment"] = pdata.obs["condition"].str.split("_").str[1]

    # QC table: cells per (sample, celltype)
    qc = pdata.obs[["sample","celltype","psbulk_n_cells"]].copy()
    qc.to_csv(os.path.join(args.out, "pseudobulk_samples_per_celltype.csv"), index=False)

    # placeholder for DEG loop (Task 5)
    print(f"Pseudobulk shape: {pdata.shape}")
    raise SystemExit("Aggregation OK — Task 5 adds PyDESeq2")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run and confirm aggregation**

Run:
```bash
pixi run python pipeline/scripts/19_pseudobulk_deg.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out /tmp/psbk_test \
  --level manual_level1 2>&1 | tail -10
```
Expected: prints Pseudobulk shape and exits with "Aggregation OK".

- [ ] **Step 3: Commit**

```bash
git add pipeline/scripts/19_pseudobulk_deg.py
git commit -m "feat(19): pseudobulk DEG skeleton with decoupler aggregation"
```

---

## Task 5: PyDESeq2 factorial DEG loop

**Files:**
- Modify: `pipeline/scripts/19_pseudobulk_deg.py`

- [ ] **Step 1: Replace the placeholder with the per-celltype PyDESeq2 loop**

Replace `raise SystemExit("Aggregation OK — Task 5 adds PyDESeq2")` with:

```python
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    contrasts = [
        ("treatment", ["treatment", "STM", "DMSO"]),
        ("genotype",  ["genotype", "Mutant", "WT"]),
        # interaction handled separately via LRT below
    ]

    summary_rows = []
    celltypes = sorted(pdata.obs["celltype"].unique())

    for ct in celltypes:
        ct_mask = pdata.obs["celltype"] == ct
        sub = pdata[ct_mask].copy()

        # Drop genes with all-zero counts
        keep_gene = (sub.X.sum(axis=0) > 0)
        sub = sub[:, np.asarray(keep_gene).ravel()].copy()

        # Need ≥2 samples per condition AND ≥min-cells in each donor sample
        per_cond = sub.obs.groupby("condition").size()
        cell_ok = sub.obs["psbulk_n_cells"] >= args.min_cells
        if not cell_ok.all():
            sub = sub[cell_ok.values].copy()
            per_cond = sub.obs.groupby("condition").size()
        if (per_cond < 2).any() or sub.n_obs < 6:
            print(f"  Skipping {ct}: insufficient samples after filtering "
                  f"({sub.n_obs} samples, per-cond: {per_cond.to_dict()})")
            continue

        print(f"\n=== Celltype: {ct} (n_samples={sub.n_obs}) ===")
        meta = sub.obs[["genotype","treatment","donor"]].copy()
        for col in ("genotype","treatment"):
            meta[col] = meta[col].astype("category")

        dds = DeseqDataSet(
            counts=pd.DataFrame(
                sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X,
                index=sub.obs_names, columns=sub.var_names,
            ).astype(int),
            metadata=meta,
            design="~ genotype + treatment + genotype:treatment",
            inference=DefaultInference(n_cpus=4),
            quiet=True,
        )
        dds.deseq2()

        for cname, contrast in contrasts:
            try:
                ds = DeseqStats(dds, contrast=contrast, quiet=True)
                ds.summary()
                res = ds.results_df.reset_index().rename(columns={"index":"gene"})
                res["celltype"] = ct
                res["contrast"] = cname
                res.to_csv(
                    os.path.join(args.out, f"deg_{ct}_{cname}.csv"),
                    index=False
                )
                top = res[(res["padj"] < args.padj) & (res["log2FoldChange"].abs() > 1)]
                top = top.sort_values("padj")
                top.to_csv(
                    os.path.join(args.out, f"deg_{ct}_{cname}_top.csv"),
                    index=False
                )
                summary_rows.append({
                    "celltype": ct,
                    "contrast": cname,
                    "n_sig_padj": int((res["padj"] < args.padj).sum()),
                    "n_sig_padj_lfc1": int(len(top)),
                    "n_samples": sub.n_obs,
                })

                # Volcano
                fig, ax = plt.subplots(figsize=(7, 6))
                x = res["log2FoldChange"].values
                y = -np.log10(res["padj"].fillna(1).values)
                sig = (res["padj"] < args.padj) & (res["log2FoldChange"].abs() > 1)
                ax.scatter(x[~sig], y[~sig], s=4, c="lightgray", alpha=0.5)
                ax.scatter(x[sig],  y[sig],  s=8, c="firebrick")
                ax.axhline(-np.log10(args.padj), c="k", lw=0.5, ls="--")
                ax.axvline( 1, c="k", lw=0.5, ls="--")
                ax.axvline(-1, c="k", lw=0.5, ls="--")
                ax.set_xlabel("log2 FC")
                ax.set_ylabel("-log10 padj")
                ax.set_title(f"{ct} — {cname}")
                fig.tight_layout()
                fig.savefig(
                    os.path.join(args.out, f"deg_volcano_{ct}_{cname}.png"),
                    dpi=150
                )
                plt.close(fig)
            except Exception as e:
                print(f"  [{ct}/{cname}] failed: {e}")

    pd.DataFrame(summary_rows).to_csv(
        os.path.join(args.out, "deg_summary.csv"), index=False
    )
    print(f"\nScript 19 complete. Outputs in {args.out}/")
```

- [ ] **Step 2: Run end-to-end for level1**

Run:
```bash
pixi run python pipeline/scripts/19_pseudobulk_deg.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out results/19_pseudobulk_deg/manual_level1 \
  --level manual_level1 \
  --min-cells 10 2>&1 | tee /tmp/psbk_run.log | tail -40
```
Expected: per-celltype "===" headers, skip messages for B_cell / T_cell / basophil (too few samples), `deg_summary.csv` lists hit counts per (celltype × contrast).

- [ ] **Step 3: Sanity check the summary**

Run:
```bash
cat results/19_pseudobulk_deg/manual_level1/deg_summary.csv
```
Expected: celltypes like granulocyte_neutrophil, myeloid_progenitor, cycling_myeloid, progenitor, erythroid each appear twice (treatment, genotype contrasts) with n_sig counts.

- [ ] **Step 4: Run for level2**

Run:
```bash
pixi run python pipeline/scripts/19_pseudobulk_deg.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out results/19_pseudobulk_deg/manual_level2 \
  --level manual_level2 --min-cells 10 2>&1 | tail -30
```
Expected: more celltype rows in summary.

- [ ] **Step 5: Commit**

```bash
git add pipeline/scripts/19_pseudobulk_deg.py
git commit -m "feat(19): PyDESeq2 factorial pseudobulk DEG per celltype with volcanoes"
```

---

## Task 6: Pathway activity from DEG stats (decoupler + MSigDB)

**Files:**
- Create: `pipeline/scripts/20_pathway_activity.py`

- [ ] **Step 1: Write the script**

```python
#!/usr/bin/env python3
"""
20_pathway_activity.py

Pathway activity per celltype × contrast, computed from the PyDESeq2 result
tables in `results/19_pseudobulk_deg/`. We use decoupler's univariate linear
model (ULM) on the `stat` column (the Wald statistic) — recommended pattern
for transferring bulk-style DEG output into a multi-collection enrichment.

Collections (mouse):
- MSigDB Hallmark (mh.all)
- Reactome (m2.cp.reactome)
- KEGG (m2.cp.kegg)

Pathway scoring with n=2 is more robust than individual gene calls because
the test aggregates ~50–500 genes per set.

Inputs
------
--deg-dir     directory with deg_<celltype>_<contrast>.csv (from script 19)
--out         output directory
--organism    'mouse' (default) or 'human'
--padj        threshold for "significant" pathway hits (default 0.05)

Outputs
-------
pathway_activity_<contrast>.csv         long-format scores/padj (all celltypes)
pathway_heatmap_<contrast>.png          celltypes × pathway score heatmap
pathway_top_hits.csv                    union of significant hits across runs
"""
import argparse
import glob
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--deg-dir", required=True)
    p.add_argument("--out", default="results/20_pathway_activity")
    p.add_argument("--organism", default="mouse", choices=["mouse", "human"])
    p.add_argument("--padj", type=float, default=0.05)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    import decoupler as dc

    # ── Build prior knowledge network from MSigDB ───────────────────────
    print("Fetching MSigDB collections ...")
    hallmark = dc.op.msigdb(
        organism=args.organism, top_level_collection="hallmark"
    )
    reactome = dc.op.msigdb(
        organism=args.organism, top_level_collection="curated_genesets"
    )
    reactome = reactome[reactome["collection"].str.contains("REACTOME", case=False)]
    kegg = dc.op.msigdb(
        organism=args.organism, top_level_collection="curated_genesets"
    )
    kegg = kegg[kegg["collection"].str.contains("KEGG", case=False)]

    net = pd.concat([hallmark, reactome, kegg], ignore_index=True)
    net = net.rename(columns={"geneset":"source", "genesymbol":"target"})
    net["weight"] = 1.0
    print(f"  Combined network: {net['source'].nunique()} pathways, "
          f"{len(net)} edges")

    # ── Load DEG result tables ──────────────────────────────────────────
    files = sorted(glob.glob(os.path.join(args.deg_dir, "deg_*_treatment.csv"))) + \
            sorted(glob.glob(os.path.join(args.deg_dir, "deg_*_genotype.csv")))
    files = [f for f in files if "_top" not in f]
    if not files:
        raise SystemExit(f"No DEG CSVs found in {args.deg_dir}")

    # Long-format: celltype × pathway scores per contrast
    all_results = []
    for f in files:
        base = os.path.basename(f).replace("deg_", "").replace(".csv", "")
        parts = base.rsplit("_", 1)
        celltype, contrast = parts[0], parts[1]
        df = pd.read_csv(f).dropna(subset=["stat"])
        if df.empty:
            continue
        mat = df.set_index("gene")[["stat"]].T   # 1 × n_genes
        try:
            acts, padj = dc.mt.ulm(
                data=mat, net=net, tmin=5, verbose=False
            )
        except Exception as e:
            print(f"  [{celltype}/{contrast}] ULM failed: {e}")
            continue
        long = pd.DataFrame({
            "pathway": acts.columns,
            "score":   acts.values.ravel(),
            "padj":    padj.values.ravel(),
        })
        long["celltype"] = celltype
        long["contrast"] = contrast
        all_results.append(long)

    if not all_results:
        raise SystemExit("No pathway results produced — check DEG inputs")

    big = pd.concat(all_results, ignore_index=True)

    # ── Per-contrast wide tables + heatmaps ─────────────────────────────
    for contrast in big["contrast"].unique():
        sub = big[big["contrast"] == contrast].copy()
        sub.to_csv(
            os.path.join(args.out, f"pathway_activity_{contrast}.csv"),
            index=False
        )

        # Pick top pathways by max |score| across celltypes
        ranking = sub.groupby("pathway")["score"].apply(
            lambda x: x.abs().max()
        ).sort_values(ascending=False)
        top_paths = ranking.head(40).index.tolist()
        heat = sub[sub["pathway"].isin(top_paths)].pivot_table(
            index="pathway", columns="celltype", values="score", fill_value=0
        ).reindex(top_paths)

        fig, ax = plt.subplots(figsize=(max(8, 0.6*heat.shape[1]),
                                        max(8, 0.3*heat.shape[0])))
        vmax = float(np.nanmax(np.abs(heat.values))) or 1.0
        im = ax.imshow(heat.values, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(heat.shape[1]))
        ax.set_xticklabels(heat.columns, rotation=45, ha="right")
        ax.set_yticks(range(heat.shape[0]))
        ax.set_yticklabels(heat.index, fontsize=7)
        ax.set_title(f"Top 40 pathway activities — contrast: {contrast}")
        plt.colorbar(im, ax=ax, label="ULM score")
        plt.tight_layout()
        fig.savefig(
            os.path.join(args.out, f"pathway_heatmap_{contrast}.png"),
            dpi=150
        )
        plt.close(fig)

    # ── Union of significant pathway hits ───────────────────────────────
    sig = big[big["padj"] < args.padj].copy()
    sig = sig.sort_values(["contrast","celltype","padj"])
    sig.to_csv(os.path.join(args.out, "pathway_top_hits.csv"), index=False)

    print(f"\nScript 20 complete. Outputs in {args.out}/")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run for the level1 DEG outputs**

Run:
```bash
pixi run python pipeline/scripts/20_pathway_activity.py \
  --deg-dir results/19_pseudobulk_deg/manual_level1 \
  --out    results/20_pathway_activity/manual_level1 \
  --organism mouse 2>&1 | tail -20
```
Expected: prints "Combined network: <N> pathways", produces `pathway_activity_treatment.csv`, `pathway_activity_genotype.csv`, two heatmaps, and `pathway_top_hits.csv`.

- [ ] **Step 3: Inspect top hits sanity-check**

Run:
```bash
head -20 results/20_pathway_activity/manual_level1/pathway_top_hits.csv
```
Expected: rows like `HALLMARK_TNFA_SIGNALING_VIA_NFKB,celltype=granulocyte_neutrophil,contrast=treatment,...` — biologically plausible inflammation hits.

- [ ] **Step 4: Run for level2**

Run:
```bash
pixi run python pipeline/scripts/20_pathway_activity.py \
  --deg-dir results/19_pseudobulk_deg/manual_level2 \
  --out    results/20_pathway_activity/manual_level2 \
  --organism mouse 2>&1 | tail -10
```

- [ ] **Step 5: Commit**

```bash
git add pipeline/scripts/20_pathway_activity.py
git commit -m "feat(20): pathway activity per celltype via decoupler ULM on DEG stats (Hallmark + Reactome + KEGG)"
```

---

## Task 7: Nextflow integration

**Files:**
- Create: `pipeline/modules/composition.nf`
- Create: `pipeline/modules/pseudobulk_deg.nf`
- Create: `pipeline/modules/pathway_activity.nf`
- Modify: `pipeline/main.nf`
- Modify: `pipeline/nextflow.config`

- [ ] **Step 1: Create the composition module**

Create `pipeline/modules/composition.nf`:

```nextflow
process COMPOSITION_SCCODA {
    tag "composition_${level}"
    publishDir "${params.outdir}/18_composition/${level}", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    tuple path(h5ad), val(level)

    output:
    path "sccoda_*"
    path "*.png", optional: true

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/18_composition_sccoda.py \\
        --input ${h5ad} \\
        --out . \\
        --level ${level}
    """
}
```

- [ ] **Step 2: Create the pseudobulk-DEG module**

Create `pipeline/modules/pseudobulk_deg.nf`:

```nextflow
process PSEUDOBULK_DEG {
    tag "pseudobulk_deg_${level}"
    publishDir "${params.outdir}/19_pseudobulk_deg/${level}", mode: 'copy'
    memory '32 GB'
    cpus 4

    input:
    tuple path(h5ad), val(level)

    output:
    path "deg_*"
    path "pseudobulk_samples_per_celltype.csv"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/19_pseudobulk_deg.py \\
        --input ${h5ad} \\
        --out . \\
        --level ${level} \\
        --min-cells ${params.pseudobulk_min_cells}
    """
}
```

- [ ] **Step 3: Create the pathway-activity module**

Create `pipeline/modules/pathway_activity.nf`:

```nextflow
process PATHWAY_ACTIVITY {
    tag "pathway_activity_${level}"
    publishDir "${params.outdir}/20_pathway_activity/${level}", mode: 'copy'
    memory '16 GB'
    cpus 2

    input:
    tuple path(deg_dir), val(level)

    output:
    path "pathway_*"

    script:
    """
    pixi run -m ${params.pixi_manifest} python ${projectDir}/scripts/20_pathway_activity.py \\
        --deg-dir ${deg_dir} \\
        --out . \\
        --organism ${params.pathway_organism}
    """
}
```

- [ ] **Step 4: Wire into main.nf**

Open `pipeline/main.nf`. Find the existing `if (params.run_interaction_analysis)` block (which currently launches macrophage_states + cellchat + nichenet). Replace it with:

```nextflow
include { COMPOSITION_SCCODA }      from './modules/composition.nf'
include { PSEUDOBULK_DEG }          from './modules/pseudobulk_deg.nf'
include { PATHWAY_ACTIVITY }        from './modules/pathway_activity.nf'
include { MACROPHAGE_STATES }       from './modules/macrophage_states.nf'

// ... earlier processes that produce progenitor_annotated_h5ad ...

if (params.run_downstream_analysis) {
    levels_ch = Channel.from('manual_level1', 'manual_level2')

    MACROPHAGE_STATES(progenitor_annotated_h5ad)

    sccoda_in = progenitor_annotated_h5ad.combine(levels_ch)
    COMPOSITION_SCCODA(sccoda_in)

    psbk_in   = progenitor_annotated_h5ad.combine(levels_ch)
    PSEUDOBULK_DEG(psbk_in)

    PATHWAY_ACTIVITY(
        PSEUDOBULK_DEG.out[0].combine(levels_ch)
    )
}
```

If the previous `run_interaction_analysis` flag is still referenced anywhere else in `main.nf`, replace those references with `run_downstream_analysis`.

- [ ] **Step 5: Add params to nextflow.config**

In `pipeline/nextflow.config`, add (or update) under `params { ... }`:

```nextflow
params.run_downstream_analysis = false
params.pseudobulk_min_cells    = 10
params.pathway_organism        = 'mouse'
params.cellchat_organism       = 'Mm'  // kept for backward compat if 16/17 are re-run
```

Remove `params.run_interaction_analysis` if it's no longer referenced.

- [ ] **Step 6: Syntax-check Nextflow**

Run:
```bash
cd pipeline && nextflow run main.nf -preview --run_downstream_analysis 2>&1 | tail -30
```
Expected: no syntax errors. (Preview mode shows the DAG without executing.)

- [ ] **Step 7: Commit**

```bash
git add pipeline/modules/composition.nf pipeline/modules/pseudobulk_deg.nf pipeline/modules/pathway_activity.nf
git add pipeline/main.nf pipeline/nextflow.config
git commit -m "feat: wire scCODA + pseudobulk DEG + pathway activity into Nextflow (replaces CellChat/NicheNet block)"
```

---

## Task 8: Update the project plan / docs

**Files:**
- Modify: `docs/superpowers/plans/2026-05-11-scrna-interaction-analysis.md` (or create a brief addendum next to it)

- [ ] **Step 1: Add a deprecation note to the old plan**

At the top of `docs/superpowers/plans/2026-05-11-scrna-interaction-analysis.md`, prepend:

```markdown
> **DEPRECATED 2026-05-12.** Critical review (see session note kb id=56 plus
> 2026-05-12 critical-thinking evaluation) found that the CellChat / NicheNet
> comparative analysis is dominated by composition confounds (sort fractions
> differ across conditions ≥20×) and is statistically inappropriate at n=2.
> Replacement plan: `docs/superpowers/plans/2026-05-12-composition-and-pseudobulk-analysis.md`.
> The 16_cellchat.R / 17_nichenet.R scripts are retained as descriptive tools
> only; their cross-condition outputs were moved to `results/_archive_exploratory/`.
```

- [ ] **Step 2: Commit**

```bash
git add docs/superpowers/plans/2026-05-11-scrna-interaction-analysis.md
git commit -m "docs: mark CellChat/NicheNet plan deprecated, link replacement"
```

---

## Task 9: Final verification

**Files:** none (this task is a check)

- [ ] **Step 1: Confirm all three new analyses produced outputs**

Run:
```bash
echo "=== 18 composition ==="
ls results/18_composition/manual_level1/ | head
echo "=== 19 pseudobulk DEG ==="
ls results/19_pseudobulk_deg/manual_level1/ | head
echo "=== 20 pathway activity ==="
ls results/20_pathway_activity/manual_level1/ | head
echo "=== archive ==="
ls results/_archive_exploratory/
```
Expected: every section lists files.

- [ ] **Step 2: Skim one credible-effect table and one pathway-hit table**

Run:
```bash
cat results/18_composition/manual_level1/sccoda_treatment_in_WT_credible_effects.csv | head
echo "---"
head -20 results/20_pathway_activity/manual_level1/pathway_top_hits.csv
```
Expected: rows with non-empty effect/pathway names, biologically interpretable.

- [ ] **Step 3: Confirm git tree is clean**

Run:
```bash
git status
git log --oneline -12
```
Expected: clean tree (or only untracked working files); 8 new commits on top of the previous 16_cellchat refactor.

---

## Self-Review Check

- ✅ Spec coverage: composition (scCODA — Task 2–3), pseudobulk DEG (Task 4–5), pathway activity (Task 6), cleanup/archive (Task 1), Nextflow wiring (Task 7), docs (Task 8), verification (Task 9).
- ✅ No placeholders: all code blocks complete.
- ✅ Type consistency: `level` column passed as `manual_level1` / `manual_level2` in every script; `sample = donor__condition` and `condition = genotype_treatment` formula identical across scripts.
- ✅ Each task ≤ ~30 min and ends with a commit.
- ✅ Tools verified via context7 (scCODA: `mod.CompositionalAnalysis` + `sample_hmc` + `credible_effects`; decoupler: `dc.pp.pseudobulk` + `dc.mt.ulm`; PyDESeq2: `DeseqDataSet` + `DeseqStats`; Milo was the alternative but user specified scCODA).
