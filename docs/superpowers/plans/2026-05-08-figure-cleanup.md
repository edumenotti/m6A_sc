# Figure Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Comment out non-central diagnostic plot blocks in the pipeline scripts, keeping only figures that are part of the final annotation story.

**Architecture:** Pure code edit — no logic changes, only wrapping plot-generating blocks with `# DIAGNOSTIC (commented out)` markers. Each block gets a one-line comment explaining why it was removed so future readers know it was intentional.

**Tech Stack:** Python, matplotlib, scanpy

---

## What stays vs. what goes

**KEEP (central to processing/annotation story):**
- `01_qc.py`: QC violin — required for any methods figure
- `03_normalize.py`: all plots (elbow/HVG scatter)
- `04_integrate.py`: training loss + batch correction UMAPs — proof of integration
- `05_cluster.py`: per-resolution leiden UMAPs — core output of clustering step
- `07_markers.py`: `heatmap_marker_scores_by_cluster.png`, `dotplot_clusters_markers.png`, `umap_watchlist_clusters.png`, `dotplot_top_de_genes.png`
- `09_manual_marker_review.py`: all scoring dotplots (used during annotation)
- `11_final_figures.py`: all (these ARE the final figures)
- `13_progenitor_recluster.py`: per-resolution UMAPs + dotplots, `prog_original_umap_overview.png`
- `14_apply_progenitor_annotation.py`: all

**COMMENT OUT:**
- `05_cluster.py:126-129`: multi-panel metadata UMAP (redundant — batch UMAPs are in 04)
- `07_markers.py:310-322`: per-prediction-column dotplot loop (intermediate ML predictions, not final annotations)
- `07_markers.py:326-330`: per-gene expression UMAP loop (exploratory; 10 files per run)
- `09_manual_marker_review.py:382-387`: per-gene expression UMAP loop (exploratory; 23 files per run)
- `13_progenitor_recluster.py:298-312`: QC overlay UMAPs + lineage marker UMAPs (exploratory; dozens of files per run)

---

## File Map

| File | Change |
|------|--------|
| `pipeline/scripts/05_cluster.py` | Comment lines 126-129 |
| `pipeline/scripts/07_markers.py` | Comment lines 310-322 and 326-330 |
| `pipeline/scripts/09_manual_marker_review.py` | Comment lines 382-387 |
| `pipeline/scripts/13_progenitor_recluster.py` | Comment lines 298-312 |

---

### Task 1: Clean 05_cluster.py — remove redundant metadata UMAP

**Files:**
- Modify: `pipeline/scripts/05_cluster.py:126-129`

- [ ] **Step 1: Verify current content at lines 126-129**

Run:
```bash
sed -n '124,132p' pipeline/scripts/05_cluster.py
```
Expected output:
```
        plt.close()

    metadata_colors = [key for key in ["population", "treatment", "donor", "pool"] if key in adata.obs]
    sc.pl.umap(adata, color=metadata_colors, show=False)
    plt.savefig(out / "umap_metadata.png", dpi=150, bbox_inches="tight")
    plt.close()

    write_cluster_composition(adata, out, cluster_keys, COMPOSITION_KEYS)
```

- [ ] **Step 2: Comment out the metadata UMAP block**

Replace lines 126-129 with:
```python
    # Batch/metadata UMAP panel — redundant with per-metadata UMAPs saved in 04_integrate.
    # metadata_colors = [key for key in ["population", "treatment", "donor", "pool"] if key in adata.obs]
    # sc.pl.umap(adata, color=metadata_colors, show=False)
    # plt.savefig(out / "umap_metadata.png", dpi=150, bbox_inches="tight")
    # plt.close()
```

- [ ] **Step 3: Verify the file parses without errors**

Run:
```bash
python -c "import ast; ast.parse(open('pipeline/scripts/05_cluster.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add pipeline/scripts/05_cluster.py
git commit -m "chore: comment out redundant metadata UMAP in 05_cluster"
```

---

### Task 2: Clean 07_markers.py — remove intermediate-prediction dotplots and per-gene UMAPs

**Files:**
- Modify: `pipeline/scripts/07_markers.py:310-330`

- [ ] **Step 1: Verify current content at lines 310-330**

Run:
```bash
sed -n '309,332p' pipeline/scripts/07_markers.py
```
Expected output includes the `for label_col in PREDICTION_COLUMNS:` block (lines 310-322) and the `for gene in ["Kit", "Ly6a", ...]` block (lines 326-330).

- [ ] **Step 2: Comment out the per-prediction-column dotplot loop**

Replace lines 310-322:
```python
    # Per-prediction-column dotplots — intermediate ML outputs (popv, scanvi, etc.), not final annotations.
    # for label_col in PREDICTION_COLUMNS:
    #     if label_col in adata.obs:
    #         sc.pl.dotplot(
    #             adata,
    #             var_names=markers,
    #             groupby=label_col,
    #             use_raw=False,
    #             layer=layer if layer in adata.layers else None,
    #             standard_scale="var",
    #             show=False,
    #         )
    #         plt.savefig(out / f"dotplot_{label_col}_markers.png", dpi=150, bbox_inches="tight")
    #         plt.close()
```

- [ ] **Step 3: Comment out the per-gene expression UMAP loop**

Replace lines 326-330:
```python
    # Per-gene expression UMAPs — exploratory during annotation; generates 10 files per run.
    # for gene in ["Kit", "Ly6a", "Hlf", "Mecom", "Mki67", "Cd79a", "Cd3d", "Gata1", "Mpo", "Csf3r"]:
    #     if gene in adata.var_names:
    #         sc.pl.umap(adata, color=gene, layer=layer if layer in adata.layers else None, show=False, vmax="p99")
    #         plt.savefig(out / f"umap_expr_{gene}.png", dpi=150, bbox_inches="tight")
    #         plt.close()
```

- [ ] **Step 4: Verify the file parses without errors**

Run:
```bash
python -c "import ast; ast.parse(open('pipeline/scripts/07_markers.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 5: Commit**

```bash
git add pipeline/scripts/07_markers.py
git commit -m "chore: comment out intermediate-prediction dotplots and per-gene UMAPs in 07_markers"
```

---

### Task 3: Clean 09_manual_marker_review.py — remove per-gene expression UMAPs

**Files:**
- Modify: `pipeline/scripts/09_manual_marker_review.py:382-387`

- [ ] **Step 1: Verify current content at lines 382-387**

Run:
```bash
sed -n '380,390p' pipeline/scripts/09_manual_marker_review.py
```
Expected output includes the `for gene in args.umap_genes.split(","):` loop.

- [ ] **Step 2: Comment out the per-gene UMAP loop**

Replace lines 382-387:
```python
    # Per-gene expression UMAPs — exploratory during manual review; generates ~23 files per run.
    # for gene in args.umap_genes.split(","):
    #     gene = gene.strip()
    #     if gene and gene in adata.var_names:
    #         sc.pl.umap(adata, color=gene, layer=args.layer if args.layer in adata.layers else None, show=False, vmax="p99")
    #         plt.savefig(out / f"umap_expr_{gene}.png", dpi=150, bbox_inches="tight")
    #         plt.close()
```

- [ ] **Step 3: Verify the file parses without errors**

Run:
```bash
python -c "import ast; ast.parse(open('pipeline/scripts/09_manual_marker_review.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add pipeline/scripts/09_manual_marker_review.py
git commit -m "chore: comment out per-gene expression UMAPs in 09_manual_marker_review"
```

---

### Task 4: Clean 13_progenitor_recluster.py — remove QC overlays and lineage marker UMAPs

**Files:**
- Modify: `pipeline/scripts/13_progenitor_recluster.py:298-312`

- [ ] **Step 1: Verify current content at lines 298-312**

Run:
```bash
sed -n '296,315p' pipeline/scripts/13_progenitor_recluster.py
```
Expected output includes the QC overlays loop and the lineage marker UMAPs loop.

- [ ] **Step 2: Comment out the QC overlays and lineage marker UMAP blocks**

Replace lines 298-312:
```python
    # QC overlays on re-embedded UMAP — exploratory during QC review; not needed in final output.
    # for col in ["pct_counts_mt", "n_genes_by_counts", "n_counts"]:
    #     if col in adata_prog.obs.columns:
    #         save_umap(adata_prog, col, out, "qc")
    # for col in ["score_Stress_IEG", "score_Cycling"]:
    #     if col in adata_prog.obs.columns:
    #         save_umap(adata_prog, col, out, "state")

    # Lineage marker UMAPs — exploratory; dozens of files per run. Dotplots cover the same info.
    # for lineage, genes in PROGENITOR_MARKERS.items():
    #     score_col = f"score_{lineage}"
    #     if score_col in adata_prog.obs.columns:
    #         save_umap(adata_prog, score_col, out, f"score_{lineage}")
    #     for g in present(adata_prog, genes[:2]):
    #         save_umap(adata_prog, g, out, f"gene_{lineage}")
```

- [ ] **Step 3: Verify the file parses without errors**

Run:
```bash
python -c "import ast; ast.parse(open('pipeline/scripts/13_progenitor_recluster.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add pipeline/scripts/13_progenitor_recluster.py
git commit -m "chore: comment out QC overlays and lineage marker UMAPs in 13_progenitor_recluster"
```

---

## Self-Review

**Spec coverage:**
- ✅ Remove non-central plots from scripts 05, 07, 09, 13
- ✅ All plot blocks to be removed are shown with exact commented-out code
- ✅ All scripts that only produce central plots (01, 03, 04, 08, 11, 12, 14) are untouched
- ✅ No logic changes — only plot-generating lines are commented

**Placeholder scan:** No TBD, no "implement later". All edits show exact before/after.

**Risk:** Commenting (not deleting) means plots can be re-enabled by uncommenting. ✅
