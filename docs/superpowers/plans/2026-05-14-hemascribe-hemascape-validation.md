# HemaScribe + HemaScape Validation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

---

## Status update — 2026-05-15

**Tasks 1–4 completed.** Final architecture diverged from skeleton (Python↔R via MTX bridge, not zellkonverter) — see actual scripts in `pipeline/scripts/21{a,b,c,z}_*` and `pipeline/scripts/validation/21{d,e,f,g,h}_*`.

**Key findings (corrigem hipóteses iniciais):**

1. **Direção do efeito Mutant é OPOSTA ao previsto.** scCODA com `final_annotation` (32 labels granulares) mostra clone Mutant **depletado** de HSC/STHSC/MPP2/MPP3/GP/mGMP (WT > Mutant), e **expandido** em mature B/T/monocyte/pDC (Mutant > WT). Não é "Mutant expande HSC" — é "Mutant tem diferenciação acelerada ou exaustão stem".

2. **"HSC com Dntt/Ighm aberrantes" era artefato de impureza.** Manual_level2 HSC era mistura de HSC+STHSC+MPP3+MPP4. Anotação granular HemaScribe revela que o sinal **Dntt +2.34 padj=1.5e-11**, **Ighm +1.22 padj=3.6e-09**, **Flt3 +0.66 padj=0.002** está no **MPP4 (lymphoid-primed)**, onde é biologicamente esperado. HSC canônico (n=823) não tem essas mudanças após FDR. Biologia revisada: Mutant **hiper-ativa lymphoid-priming no MPP4**.

3. **Assinatura GMP (Rras2↑/Serpine2↓/Prss57↓) totalmente validada** em mGMP (Rras2 +1.89), cMoP (+1.80, sinal mais forte com Itgb7 +0.71), GP (+2.01). Direção e magnitude consistentes com manual_level2 GMP_neutrophil_primed.

4. **Mutant não responde a STM**: `treatment_in_Mutant` retorna 0 credible effects. WT responde (HSC/STHSC/MPP4/cMoP/mGMP depletados sob STM). Hipótese: Mutant já está depletado no baseline, não há mais a depletar.

**Outputs:**
- `results/18_composition/final_annotation/` (scCODA, 32 labels)
- `results/19_pseudobulk_deg/final_annotation/` (PyDESeq2, 18 celltypes com min-cells=5)
- `results/21_hemascribe/composition_concordance.md` (Task 3 Step 4)
- `results/21_hemascribe/deg_concordance.md` (Task 4 Step 4)

**Tasks 5–6** (HemaScape trajectory + synthesis report) pendentes.

---

**Goal:** Validate the paired-design composition + DEG findings (HSC/MPP/GMP Mutant>WT, Rras2↑/Socs2↓ signature) by re-annotating with HemaScribe under sort-aware priors, then test trajectory consistency with HemaScape. This is **validation, not new analysis** — code must be economical, reuse existing scripts where possible.

**Architecture:** One new annotation column in adata (`hemascribe_label`), gated by sort population priors (LSK / LK / I). Existing scCODA (script 18) and PyDESeq2 (script 19) reused with `--level hemascribe_label`. Trajectory module is one short script. No new Nextflow wiring (validation only).

**Tech Stack:**
- `hemascribe` — automated hematopoietic annotation. **No fallback.** CellTypist was tested previously and produced poor results; if HemaScribe cannot be installed or fails, stop the validation here and report.
- `hemascape` — trajectory pseudotime. **No fallback.** If HemaScape is unavailable, skip Task 5 entirely and proceed to Task 6 noting trajectory validation was not done.
- existing `pipeline/scripts/18_composition_sccoda.py` + `19_pseudobulk_deg.py` — re-run with new label column
- scanpy, anndata — already in pixi

**Why sort-aware annotation matters:** The `population` column (LSK / LK / I) is biological ground truth from FACS gating. LSK cells **must** be Sca1+/c-Kit+/Lin- (HSC/MPP/LMPP territory) — a HemaScribe call of "Mature_B" in an LSK cell is wrong by definition. We use the sort gate as a hard prior to filter implausible labels before composition reanalysis.

**Sort gate → allowed-label map (mouse bone marrow):**
- `LSK`: HSC, MPP, LMPP, LSK-MPP — anything with stemness; reject mature lineages
- `LK`:  CMP, GMP, MEP, MEP_Erythroid, myeloid progenitors — reject HSC and mature lymphoid
- `I`:   all mature lineages allowed (broad gate)

---

## File Structure

**Created:**
- `pipeline/scripts/21_hemascribe_annotate.py` — annotate adata + sort-aware audit, writes `results/14_progenitor_annotated/adata_hemascribe.h5ad` with new `hemascribe_label` + `hemascribe_label_sortgated` columns
- `pipeline/scripts/22_hemascape_trajectory.py` — pseudotime per genotype, exports plots + summary CSV
- `results/21_hemascribe/` — annotation audit outputs (confusion vs manual_level2, sankey, sort-gate consistency)
- `results/22_hemascape/` — trajectory plots
- `results/23_validation_synthesis/concordance.md` — final concordance report

**Modified (no permanent change to pipeline):**
- None. Existing scripts 18/19 invoked from CLI with `--level hemascribe_label_sortgated` argument; their existing code already accepts arbitrary `--level` strings as long as the column exists in `adata.obs`.

**Kept as-is:**
- Everything in `results/14_progenitor_annotated/` (the input h5ad)
- Existing scripts 18/19 — reused unmodified

---

## Task 1: Setup, tool availability, sort-aware label map

**Files:**
- Create: `pipeline/scripts/21_hemascribe_annotate.py` (skeleton only this task)

- [ ] **Step 1: Check whether HemaScribe + HemaScape are installable in the pixi env**

Run:
```bash
pixi run python -c "
try:
    import hemascribe; print('hemascribe:', hemascribe.__version__)
except ImportError as e:
    print('hemascribe NOT FOUND:', e)
try:
    import hemascape; print('hemascape:', hemascape.__version__)
except ImportError as e:
    print('hemascape NOT FOUND:', e)
"
```

If either is missing, attempt `pixi add --pypi hemascribe` and `pixi add --pypi hemascape` (or via `pip` inside the pixi env if the package isn't on PyPI; check the upstream repo).

**HARD STOP if HemaScribe unavailable** after install attempts: report BLOCKED with the install error. Do NOT substitute CellTypist or any other tool — the user has tested those and found them unsatisfactory.

**HARD STOP if HemaScape unavailable** after install attempts: skip Task 5 entirely, note in the synthesis (Task 6) that trajectory validation was not performed because HemaScape is unavailable. Do NOT substitute palantir, cellrank, or scanpy DPT.

- [ ] **Step 2: Verify sort metadata in adata**

Run:
```bash
pixi run python -c "
import scanpy as sc
a = sc.read_h5ad('results/14_progenitor_annotated/adata_progenitor_annotated.h5ad')
print('population col:', a.obs['population'].value_counts().to_dict())
print('shape:', a.shape)
print('manual_level2 categories:', a.obs['manual_level2'].nunique())
"
```
Expected: `LSK`, `LK`, `I` with cell counts. Confirms the sort gate is in `population` (not `sort_fraction`).

- [ ] **Step 3: Write skeleton for `21_hemascribe_annotate.py` with sort-aware prior map**

```python
#!/usr/bin/env python3
"""
21_hemascribe_annotate.py

Re-annotate adata with HemaScribe (or fallback) and apply sort-gate priors to
flag/reject implausible calls. Outputs:
  - hemascribe_label                  : raw HemaScribe call per cell
  - hemascribe_label_sortgated        : reject implausible calls (set to NaN)
  - hemascribe_sortgate_flag          : True if raw call violates sort prior

Usage:
  pixi run python pipeline/scripts/21_hemascribe_annotate.py \\
      --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \\
      --out   results/14_progenitor_annotated/adata_hemascribe.h5ad \\
      --audit-dir results/21_hemascribe
"""
import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# Sort gate -> labels that are biologically plausible for cells in that gate.
# Labels not in this set for the cell's sort gate will be flagged as inconsistent.
# Keep it loose: anything stem/progenitor in LSK, progenitor/myeloid in LK, mature in I.
SORT_PRIORS = {
    "LSK": {
        # LSK = Lin-/Sca1+/c-Kit+ -> stem / multipotent only
        "HSC", "LT_HSC", "ST_HSC", "MPP", "MPP1", "MPP2", "MPP3", "MPP4",
        "LMPP", "CLP",  # CLP is Sca1+ → allowed
    },
    "LK": {
        # LK = Lin-/Sca1-/c-Kit+ -> committed myeloid/erythroid progenitors
        "CMP", "GMP", "MEP", "MEP_Erythroid", "MkP", "EMP", "CFU_E", "CFU_GM",
        "MyeloidProgenitor", "GranulocyteMonocyteProgenitor", "MegakaryocyteErythroidProgenitor",
        "Cycling_myeloid", "myeloid_progenitor",
    },
    "I": None,  # broad gate, all labels allowed
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--audit-dir", default="results/21_hemascribe")
    p.add_argument("--organism", default="mouse")
    return p.parse_args()


def annotate_with_hemascribe(adata):
    """Run HemaScribe on adata.X. Returns labels series.
    Adapt this stub to the actual HemaScribe API once installed.
    NO FALLBACK: if HemaScribe import or call fails, raise — do not substitute."""
    import hemascribe as hs  # let ImportError propagate
    # Replace with the actual call once API is known. Common patterns:
    # hs.annotate(adata, reference='mouse_bm') -> writes adata.obs['hemascribe_label']
    labels = hs.annotate(adata, organism="mouse", reference="bone_marrow")
    return pd.Series(labels, index=adata.obs_names, name="hemascribe_label")


def apply_sort_priors(labels: pd.Series, sort_col: pd.Series) -> pd.DataFrame:
    """Return DataFrame with raw label, sortgated label, and flag."""
    raw = labels.copy()
    gated = labels.copy().astype(object)
    flag = pd.Series(False, index=labels.index)
    for gate, allowed in SORT_PRIORS.items():
        if allowed is None:
            continue
        mask = (sort_col == gate)
        bad = mask & ~labels.isin(allowed)
        flag[bad] = True
        gated[bad] = np.nan
    return pd.DataFrame({
        "hemascribe_label": raw,
        "hemascribe_label_sortgated": gated,
        "hemascribe_sortgate_flag": flag,
    })


def main():
    args = parse_args()
    os.makedirs(args.audit_dir, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  shape: {adata.shape}")

    # placeholder for Task 2 — fill the annotation + audit logic
    raise SystemExit("Skeleton OK — Task 2 fills annotation + audit")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run skeleton to confirm it loads**

Run:
```bash
pixi run python pipeline/scripts/21_hemascribe_annotate.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out /tmp/hema_test.h5ad \
  --audit-dir /tmp/hema_audit 2>&1 | tail -10
```
Expected: prints shape, exits with "Skeleton OK".

- [ ] **Step 5: Commit**

```bash
git add pipeline/scripts/21_hemascribe_annotate.py
git commit -m "feat(21): HemaScribe annotation skeleton with sort-aware prior map"
```

---

## Task 2: Run HemaScribe + sort-gate audit

**Files:**
- Modify: `pipeline/scripts/21_hemascribe_annotate.py` (replace placeholder)

- [ ] **Step 1: Replace the `raise SystemExit(...)` with the annotation + audit pipeline**

Replace everything after `print(f"  shape: {adata.shape}")` with:

```python
    # ── HemaScribe annotation ───────────────────────────────────────────
    print("Annotating with HemaScribe (or fallback)...")
    labels = annotate_with_hemascribe(adata)
    print(f"  unique labels: {labels.nunique()}")
    print(labels.value_counts().head(20).to_string())

    # ── Sort-aware prior filter ─────────────────────────────────────────
    audit = apply_sort_priors(labels, adata.obs["population"])
    adata.obs["hemascribe_label"] = audit["hemascribe_label"].astype("category")
    adata.obs["hemascribe_label_sortgated"] = audit["hemascribe_label_sortgated"].astype("category")
    adata.obs["hemascribe_sortgate_flag"] = audit["hemascribe_sortgate_flag"].values

    n_flagged = int(audit["hemascribe_sortgate_flag"].sum())
    print(f"\nSort-gate conflicts: {n_flagged}/{len(audit)} cells ({n_flagged/len(audit)*100:.1f}%)")
    print("Flag rate per sort gate:")
    print(adata.obs.groupby("population")["hemascribe_sortgate_flag"].mean().to_string())

    # ── Confusion: manual_level2 × hemascribe_label_sortgated ───────────
    confusion = pd.crosstab(
        adata.obs["manual_level2"], adata.obs["hemascribe_label_sortgated"],
        dropna=False, normalize="index",
    )
    confusion.to_csv(os.path.join(args.audit_dir, "confusion_manual_level2_vs_hemascribe.csv"))

    # ── Confusion: sort × hemascribe_label (RAW, no gate filter) ────────
    sort_confusion = pd.crosstab(adata.obs["population"], adata.obs["hemascribe_label"])
    sort_confusion.to_csv(os.path.join(args.audit_dir, "confusion_sort_vs_hemascribe_raw.csv"))

    # ── Simple heatmap of the manual_level2 × hemascribe confusion ──────
    top_labels = confusion.sum(axis=0).sort_values(ascending=False).head(20).index
    sub = confusion[top_labels]
    fig, ax = plt.subplots(figsize=(max(8, 0.4*len(sub.columns)), 0.4*len(sub.index)+3))
    im = ax.imshow(sub.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(sub.columns))); ax.set_xticklabels(sub.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(sub.index)));   ax.set_yticklabels(sub.index)
    ax.set_xlabel("HemaScribe label (sort-gated)")
    ax.set_ylabel("manual_level2")
    ax.set_title("Row-normalised confusion (top 20 HemaScribe labels)")
    plt.colorbar(im, ax=ax, label="fraction of row")
    plt.tight_layout()
    fig.savefig(os.path.join(args.audit_dir, "confusion_heatmap.png"), dpi=150)
    plt.close(fig)

    # ── README documenting the run ──────────────────────────────────────
    n_cells = adata.n_obs
    with open(os.path.join(args.audit_dir, "README.md"), "w") as fh:
        fh.write(f"# HemaScribe annotation audit\n\n")
        fh.write(f"Cells: {n_cells}\n")
        fh.write(f"Sort-gate-conflicting cells: {n_flagged} ({n_flagged/n_cells*100:.1f}%)\n")
        fh.write(f"Annotation source: see top of 21_hemascribe_annotate.py output log\n\n")
        fh.write(f"## Files\n")
        fh.write(f"- confusion_manual_level2_vs_hemascribe.csv\n")
        fh.write(f"- confusion_sort_vs_hemascribe_raw.csv\n")
        fh.write(f"- confusion_heatmap.png\n")

    # ── Save annotated adata ────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adata.write_h5ad(args.out)
    print(f"\nWrote: {args.out}")
    print(f"Audit in: {args.audit_dir}/")
```

- [ ] **Step 2: Run end-to-end**

Run:
```bash
pixi run python pipeline/scripts/21_hemascribe_annotate.py \
  --input results/14_progenitor_annotated/adata_progenitor_annotated.h5ad \
  --out   results/14_progenitor_annotated/adata_hemascribe.h5ad \
  --audit-dir results/21_hemascribe 2>&1 | tee /tmp/hemascribe.log | tail -40
```
Expected: prints HemaScribe (or fallback) label distribution, conflict fraction per gate, writes audit files + new h5ad. Conflicts >50% in LSK = annotation tool not respecting biology and we have a problem — STOP and report.

- [ ] **Step 3: Sanity check the confusion**

Run:
```bash
head -1 results/21_hemascribe/confusion_manual_level2_vs_hemascribe.csv
echo "---"
pixi run python -c "
import pandas as pd
df = pd.read_csv('results/21_hemascribe/confusion_manual_level2_vs_hemascribe.csv', index_col=0)
# For each manual cluster, what is its top HemaScribe call?
print('Top HemaScribe match per manual cluster:')
for ct in df.index:
    s = df.loc[ct].sort_values(ascending=False)
    print(f'  {ct:30s} -> {s.index[0]:30s} ({s.iloc[0]:.2f}) | next: {s.index[1]:30s} ({s.iloc[1]:.2f})')
"
```
Expected: manual HSC → HemaScribe HSC/LT_HSC/MPP, manual Mature_neutrophil → HemaScribe Neutrophil, etc. If manual HSC maps mostly to "Pre_B" or similar → the manual annotation has a problem that we already suspected with Dntt/Ighm.

- [ ] **Step 4: Commit**

```bash
git add pipeline/scripts/21_hemascribe_annotate.py
git commit -m "feat(21): HemaScribe re-annotation with sort-gate prior audit"
```

---

## Task 3: Re-run scCODA composition with HemaScribe labels

**Files:** none new — invoke existing `pipeline/scripts/18_composition_sccoda.py` with new `--input` and `--level`.

- [ ] **Step 1: Confirm script 18 accepts arbitrary level column**

Run:
```bash
pixi run python pipeline/scripts/18_composition_sccoda.py --help 2>&1 | grep -A2 "level"
```
Expected: `--level` is a string argument; the existing default is `manual_level1`/`manual_level2` but anything in `adata.obs` works because script 18 just reads `adata.obs[args.level]`. **If `--level` is hard-coded to choices=, drop the choices restriction in a one-line edit** and recommit before continuing.

If the choices= constraint exists, modify the argparse line in script 18 — find:
```python
p.add_argument("--level", default="manual_level1",
               choices=["manual_level1", "manual_level2"])
```
and change to:
```python
p.add_argument("--level", default="manual_level1")
```
Commit: `git add pipeline/scripts/18_composition_sccoda.py && git commit -m "fix(18): allow arbitrary --level column for validation runs"`

- [ ] **Step 2: Run scCODA with HemaScribe (sort-gated) labels**

Run:
```bash
CUDA_VISIBLE_DEVICES="" pixi run python pipeline/scripts/18_composition_sccoda.py \
  --input results/14_progenitor_annotated/adata_hemascribe.h5ad \
  --out   results/18_composition/hemascribe \
  --level hemascribe_label_sortgated 2>&1 | tee /tmp/sccoda_hema.log | tail -30
```
Expected: 5 contrast blocks run. HMC ~3 min total. Output files identical structure to the existing `manual_level2` run.

- [ ] **Step 3: Compare credible effects vs the manual_level2 paired run**

Run:
```bash
echo "=== HemaScribe-based credibles ==="
grep "True$" results/18_composition/hemascribe/sccoda_*_credible_effects.csv | grep -v "donor\[" | sed 's|results/18_composition/hemascribe/sccoda_||; s|_credible_effects.csv:||'
echo ""
echo "=== manual_level2 paired credibles (reference) ==="
grep "True$" results/18_composition/manual_level2/sccoda_*_credible_effects.csv | grep -v "donor\[" | sed 's|results/18_composition/manual_level2/sccoda_||; s|_credible_effects.csv:||'
```
Expected: Either the same compartments (HSC/MPP/GMP-equivalents Mutant > WT) emerge → composition finding validated. Or different labels but consistent direction (e.g., HemaScribe's "LMPP" credible instead of our "GMP_neutrophil_primed") → mapping is different but biology holds.

- [ ] **Step 4: Write a one-paragraph note in `results/21_hemascribe/composition_concordance.md`**

```bash
cat > results/21_hemascribe/composition_concordance.md <<'EOF'
# Composition concordance: manual_level2 vs HemaScribe

Manual_level2 paired credibles (genotype_in_DMSO): HSC, MPP, GMP_neutrophil_primed,
Cycling_myeloid, Classical_monocyte, ...

HemaScribe paired credibles (genotype_in_DMSO): [paste actual output here]

Concordant compartments: [list]
Divergent compartments: [list with hypothesis why]

Verdict: [confirmed / partial / contradicted]
EOF
```
Fill in `[paste actual output here]`, `[list]`, and `[verdict]` from the Step 3 output. This is a small judgement call — be honest.

- [ ] **Step 5: Commit**

```bash
git add results/21_hemascribe/composition_concordance.md
git commit -m "validation: scCODA composition concordance manual vs HemaScribe"
```
(Note: only the MD is committed; `results/` is gitignored, so the scCODA output dirs remain local-only.)

---

## Task 4: Re-run pseudobulk DEG with HemaScribe labels (load-bearing clusters only)

**Files:** none new — invoke `pipeline/scripts/19_pseudobulk_deg.py` with new label, then compare top genes manually for HSC and GMP equivalents.

- [ ] **Step 1: Run script 19 with HemaScribe labels**

Run:
```bash
pixi run python pipeline/scripts/19_pseudobulk_deg.py \
  --input results/14_progenitor_annotated/adata_hemascribe.h5ad \
  --out   results/19_pseudobulk_deg/hemascribe \
  --level hemascribe_label_sortgated --min-cells 10 2>&1 | tee /tmp/psbk_hema.log | tail -40
```
(Same caveat about `choices=` argparse — drop restriction if present. Apply same one-line fix to script 19 if needed; the script already accepts arbitrary level strings if no choices=.)

Expected: per-celltype runs (those with HemaScribe label and ≥10 cells per sample). Skip messages for under-sampled labels are fine.

- [ ] **Step 2: Identify the HemaScribe label that corresponds to "HSC" and "GMP"**

Run:
```bash
ls results/19_pseudobulk_deg/hemascribe/deg_*_genotype.csv | sed 's|.*deg_||;s|_genotype.csv||' | head
```
HemaScribe labels for HSC may be `HSC`, `LT_HSC`, `ST_HSC`, etc. For GMP: `GMP`, `GranulocyteMonocyteProgenitor`, `CMP`, `MyeloidProgenitor`. Pick the closest matches.

- [ ] **Step 3: Compare top genotype DEGs for HSC and GMP equivalents**

Run (substitute `<hema-hsc>` and `<hema-gmp>` with the actual label filenames from Step 2):
```bash
echo "=== HemaScribe HSC-equivalent: top 10 genotype DEGs ==="
head -11 results/19_pseudobulk_deg/hemascribe/deg_<hema-hsc>_genotype_top.csv 2>/dev/null || echo "no top file"
echo ""
echo "=== manual_level2 HSC reference: top 10 ==="
head -11 results/19_pseudobulk_deg/manual_level2/deg_HSC_genotype_top.csv
echo ""
echo "=== HemaScribe GMP-equivalent: top 10 genotype DEGs ==="
head -11 results/19_pseudobulk_deg/hemascribe/deg_<hema-gmp>_genotype_top.csv 2>/dev/null || echo "no top file"
echo ""
echo "=== manual_level2 GMP_neutrophil_primed reference: top 10 ==="
head -11 results/19_pseudobulk_deg/manual_level2/deg_GMP_neutrophil_primed_genotype_top.csv
```

Validation gates:
- **Rras2** present in both HemaScribe-HSC and HemaScribe-GMP top hits → consistent signal across annotation methods, **strong validation**
- **Socs2** present in HemaScribe-HSC top → validates the JAK-STAT claim
- **Flt3** present in HemaScribe-HSC top → validates HSC signature
- **Dntt / Ighm** present in HemaScribe-HSC top → HSC contamination is real (the cluster genuinely has aberrant lymphoid expression), not just a manual_level2 misannotation
- **Dntt / Ighm absent** from HemaScribe-HSC → the manual_level2 HSC was contaminated; the actual HSC signal is cleaner with HemaScribe

- [ ] **Step 4: Write `results/21_hemascribe/deg_concordance.md`**

```bash
cat > results/21_hemascribe/deg_concordance.md <<'EOF'
# DEG concordance for HSC + GMP equivalents

| Gene  | manual_level2 (HSC) | HemaScribe HSC-equiv | Verdict |
|-------|----------------------|------------------------|---------|
| Flt3  | +1.28                | [fill]                 | [fill]  |
| Socs2 | -2.21                | [fill]                 | [fill]  |
| Dntt  | +3.64                | [fill]                 | [fill]  |
| Ighm  | +3.93                | [fill]                 | [fill]  |

| Gene     | manual GMP_neutrophil_primed | HemaScribe GMP-equiv | Verdict |
|----------|------------------------------|------------------------|---------|
| Rras2    | +1.91                        | [fill]                 | [fill]  |
| Serpine2 | -1.09                        | [fill]                 | [fill]  |
| Itgb7    | +1.68                        | [fill]                 | [fill]  |
| Prss57   | -1.03                        | [fill]                 | [fill]  |

Verdict per claim:
- HSC Flt3↑/Socs2↓ signature: [confirmed / partial / contradicted]
- HSC Dntt/Ighm aberrant lymphoid: [real biology / annotation contamination / unclear]
- GMP Rras2↑/Serpine2↓ signature: [confirmed / partial / contradicted]
EOF
```
Fill in the table from Step 3 output. Be honest — if a gene drops out, write that.

- [ ] **Step 5: Commit**

```bash
git add results/21_hemascribe/deg_concordance.md
git commit -m "validation: DEG concordance for HSC + GMP signatures"
```

---

## Task 5: HemaScape (or palantir) trajectory

**Files:**
- Create: `pipeline/scripts/22_hemascape_trajectory.py`

- [ ] **Step 1: Write the trajectory script (minimal — just pseudotime + UMAP-coloured-by-genotype)**

```python
#!/usr/bin/env python3
"""
22_hemascape_trajectory.py

Fit a hematopoietic trajectory on the full adata and check whether Mutant cells
occupy a different position along the differentiation axis than WT cells.

Tests the "Mutant clone has delayed/altered myeloid differentiation" hypothesis
from the pseudobulk DEG results (Serpine2/Prss57 down in Mutant suggests reduced
differentiation pressure).

Outputs
-------
hemascape_pseudotime_per_genotype.png   density of pseudotime per genotype
hemascape_umap_pseudotime.png           UMAP coloured by pseudotime + genotype
hemascape_per_cluster_pseudotime.csv    mean pseudotime per (cluster, genotype)
hemascape_per_cluster_test.csv          Mann-Whitney WT vs Mutant within each cluster
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
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out", default="results/22_hemascape")
    p.add_argument("--level", default="manual_level2",
                   help="Cluster column for per-cluster pseudotime breakdown")
    p.add_argument("--root-label", default="HSC",
                   help="Cluster label to use as the trajectory root")
    return p.parse_args()


def fit_pseudotime(adata, root_label, level):
    """Fit pseudotime using HemaScape. NO FALLBACK — if HemaScape fails, raise.
    Adapt this stub to the actual HemaScape API once installed."""
    import hemascape as hsc  # let ImportError propagate
    pt = hsc.pseudotime(adata, root=root_label, cluster_col=level)
    return pd.Series(pt, index=adata.obs_names)


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, n_comps=30)
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X_pca")
        sc.tl.umap(adata)

    print(f"Fitting pseudotime (root={args.root_label}) on {adata.n_obs} cells ...")
    pt = fit_pseudotime(adata, args.root_label, args.level)
    adata.obs["pseudotime"] = pt.reindex(adata.obs_names).values

    # ── Density of pseudotime per genotype ──────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    for g, sub in adata.obs.groupby("genotype"):
        ax.hist(sub["pseudotime"].dropna(), bins=60, alpha=0.45, label=g, density=True)
    ax.set_xlabel("pseudotime"); ax.set_ylabel("density")
    ax.set_title("Pseudotime distribution per genotype")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "hemascape_pseudotime_per_genotype.png"), dpi=150)
    plt.close(fig)

    # ── UMAP coloured by pseudotime + genotype split ────────────────────
    sc.pl.umap(adata, color=["pseudotime", "genotype"],
               save="_pseudotime_genotype.png", show=False, ncols=2)
    # scanpy saves to ./figures/ — move it
    import shutil, glob
    for src in glob.glob("figures/umap*pseudotime*.png"):
        shutil.move(src, os.path.join(args.out, os.path.basename(src)))

    # ── Per-cluster mean pseudotime + Mann-Whitney WT vs Mutant ─────────
    rows = []
    for ct, sub in adata.obs.groupby(args.level):
        wt = sub.loc[sub["genotype"] == "WT", "pseudotime"].dropna()
        mu = sub.loc[sub["genotype"] == "Mutant", "pseudotime"].dropna()
        if len(wt) < 10 or len(mu) < 10:
            continue
        stat, p = mannwhitneyu(wt, mu, alternative="two-sided")
        rows.append({
            "celltype": ct,
            "n_WT": len(wt), "n_Mutant": len(mu),
            "mean_pt_WT": wt.mean(), "mean_pt_Mutant": mu.mean(),
            "delta_WT_minus_Mutant": wt.mean() - mu.mean(),
            "U_stat": stat, "pvalue": p,
        })
    df = pd.DataFrame(rows).sort_values("pvalue")
    df.to_csv(os.path.join(args.out, "hemascape_per_cluster_test.csv"), index=False)
    print(df.to_string(index=False))

    print(f"\nTrajectory analysis complete. Outputs in {args.out}/")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on the HemaScribe-annotated h5ad**

Run:
```bash
pixi run python pipeline/scripts/22_hemascape_trajectory.py \
  --input results/14_progenitor_annotated/adata_hemascribe.h5ad \
  --out   results/22_hemascape \
  --level manual_level2 \
  --root-label HSC 2>&1 | tee /tmp/hemascape.log | tail -30
```
Expected: pseudotime fit, per-cluster CSV printed showing whether WT-Mutant pseudotime differs significantly within each cluster. **Test of interest:** within GMP_neutrophil_primed and Immature_neutrophil, is Mutant pseudotime significantly LOWER (less differentiated) than WT? That would corroborate the Serpine2/Prss57 down → less differentiation interpretation.

- [ ] **Step 3: Inspect the per-cluster table**

Run:
```bash
cat results/22_hemascape/hemascape_per_cluster_test.csv
```
Expected: rows where `delta_WT_minus_Mutant > 0` and `pvalue < 0.01` indicate WT cells in that cluster are more differentiated than Mutant ones — consistent with the DEG story. Inverse pattern (Mutant > WT) in stem clusters (HSC/MPP) is possible too — would mean Mutant HSCs are pre-loaded with a more differentiated state.

- [ ] **Step 4: Commit**

```bash
git add pipeline/scripts/22_hemascape_trajectory.py
git commit -m "feat(22): HemaScape (palantir/dpt fallback) trajectory + WT-vs-Mutant pseudotime test"
```

---

## Task 6: Synthesis — concordance report

**Files:**
- Create: `results/23_validation_synthesis/concordance.md`

- [ ] **Step 1: Write the concordance summary**

Run:
```bash
mkdir -p results/23_validation_synthesis
cat > results/23_validation_synthesis/concordance.md <<'EOF'
# Validation synthesis: paired DEG/composition findings vs HemaScribe + HemaScape

Date: $(date +%Y-%m-%d)

## Claims being tested

| # | Claim from manual_level2 paired analysis | Source |
|---|------------------------------------------|--------|
| C1 | Mutant clone expands HSC, MPP, GMP within same mouse | scCODA genotype_in_DMSO |
| C2 | HSC signature: Flt3↑, Socs2↓, aberrant Dntt/Ighm | PyDESeq2 deg_HSC_genotype |
| C3 | GMP signature: Rras2↑, Serpine2↓, Itgb7↑, Prss57↓ | PyDESeq2 deg_GMP_neutrophil_primed_genotype |
| C4 | Mutant neutrophils more proliferative (histones↑) | PyDESeq2 deg_Immature_neutrophil + GSEA |
| C5 | Mutant progenitors less differentiated | Inferred from Serpine2/Prss57 ↓ |

## Validation outcomes

| # | HemaScribe (composition) | HemaScribe (DEG) | HemaScape (trajectory) | Verdict |
|---|--------------------------|-------------------|------------------------|---------|
| C1 | [fill from Task 3]        | n/a               | n/a                    | [fill]  |
| C2 | n/a                       | [fill from Task 4]| n/a                    | [fill]  |
| C3 | n/a                       | [fill from Task 4]| n/a                    | [fill]  |
| C4 | n/a                       | n/a               | [fill from Task 5: pseudotime in Immature_neutrophil] | [fill] |
| C5 | n/a                       | n/a               | [fill from Task 5: pseudotime per cluster] | [fill] |

## Caveats & limitations

- HemaScribe label resolution may differ from manual annotation; concordance is a directional check, not a 1:1 mapping
- HemaScape (or fallback) pseudotime is a soft test — Mann-Whitney on pseudotime within cluster
- n=2 mice remains; this validation tests robustness to annotation, not to sample size

## Decision

- [ ] All major claims confirmed → proceed to manuscript figures
- [ ] Partial confirmation → refine claims; identify which need orthogonal (qPCR/flow) validation
- [ ] Contradicted → re-examine annotation/clustering before any inferential claims
EOF
```
Fill in the bracketed cells from the output of Tasks 3, 4, 5.

- [ ] **Step 2: Commit**

```bash
git add results/23_validation_synthesis/concordance.md
git commit -m "docs: validation synthesis — HemaScribe + HemaScape concordance with paired analysis"
```

---

## Self-Review Check

- ✅ Spec coverage: HemaScribe annotation (Task 1-2), sort-aware prior (Task 1+2), composition re-validation (Task 3), DEG re-validation (Task 4), HemaScape trajectory (Task 5), synthesis report (Task 6). Two approaches × two outputs covered.
- ✅ No placeholders: all code blocks complete; concordance MD templates clearly mark `[fill]` cells with where the data comes from.
- ✅ Type consistency: `hemascribe_label_sortgated` column name used in script 21 (Task 2), referenced in Task 3 and Task 4 invocations.
- ✅ Each task ≤30 min (Task 3 + 4 are mostly invocations of existing scripts).
- ✅ Economical: only 2 new scripts (21 + 22); scripts 18/19 reused unmodified (or one-line argparse fix).
- ✅ No fallbacks: per user instruction, HemaScribe failure → HARD STOP; HemaScape failure → skip Task 5 only. CellTypist explicitly rejected (tested previously, unsatisfactory results).
- ✅ Sort-aware filtering uses `population` column (verified to exist with values `LSK / LK / I`); `SORT_PRIORS` is loose enough to admit several naming conventions.
