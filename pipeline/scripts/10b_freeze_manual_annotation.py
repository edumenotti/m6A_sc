#!/usr/bin/env python
"""Stamp the frozen canonical per-cell labels onto a fresh object by barcode.

Two per-cell labels cannot be re-derived by the pipeline and must be transferred
from the canonical object:
  - manual_level1/level2: a human, cluster-ID-based annotation made against one
    specific clustering run. Re-running scVI + Leiden on different hardware drifts
    the partition (different cluster count/IDs), so a `cluster_id -> label` map can
    never be reapplied reproducibly.
  - genotype (WT/Mutant): an in-silico per-cell assignment that is NOT encoded in
    sample_id (which is donor_treatment_population_45_replicate), so it cannot be
    recovered from sample metadata. Downstream 15/18/19 read it directly.

Cell barcodes are identical across runs (same input matrix), so we transfer these
frozen labels by barcode instead of re-deriving them. This replaces the fragile
chain APPLY_MANUAL_ANNOTATION -> SUBSET_RECLUSTER -> PROGENITOR_RECLUSTER ->
APPLY_PROGENITOR_ANNOTATION with a single deterministic join, producing the same
adata_progenitor_annotated.h5ad the HemaScribe chain consumes. Expression and
embeddings come fresh from the (reproducible) pipeline; only these labels are frozen.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc

FROZEN_COLS = ["manual_level1", "manual_level2", "manual_annotation_confidence", "genotype"]


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    frozen = pd.read_csv(args.frozen, index_col=0, dtype=str)
    missing_cols = [c for c in FROZEN_COLS if c not in frozen.columns]
    if missing_cols:
        raise ValueError(f"Frozen annotation is missing columns: {missing_cols}")

    fresh_bc = adata.obs_names
    frozen_bc = frozen.index
    inter = fresh_bc.intersection(frozen_bc)

    n_fresh, n_frozen, n_inter = len(fresh_bc), len(frozen_bc), len(inter)
    fresh_only = n_fresh - n_inter          # fresh cells with no frozen label (dropped)
    frozen_only = n_frozen - n_inter        # canonical cells absent from fresh (QC drift)
    coverage = n_inter / n_frozen if n_frozen else 0.0

    print(f"[freeze] fresh cells:            {n_fresh}")
    print(f"[freeze] frozen (canonical):     {n_frozen}")
    print(f"[freeze] barcode intersection:   {n_inter}")
    print(f"[freeze] fresh-only (dropped):   {fresh_only}")
    print(f"[freeze] frozen-only (missing):  {frozen_only}")
    print(f"[freeze] frozen coverage:        {coverage:.4%}")

    # A small amount of drift is expected if DoubletFinder calls differ slightly
    # across runs; a large gap means the upstream is not reproducing the canonical
    # cell set and the annotation transfer would be unreliable.
    if coverage < args.min_coverage:
        raise ValueError(
            f"Only {coverage:.2%} of the {n_frozen} canonical cells are present in the "
            f"fresh object (min required {args.min_coverage:.0%}). Upstream QC/DoubletFinder "
            f"is not reproducing the canonical cell set; aborting rather than annotating a "
            f"divergent subset."
        )

    # Keep only cells that carry a frozen label, preserving the fresh object's order.
    adata = adata[adata.obs_names.isin(inter)].copy()
    aligned = frozen.loc[adata.obs_names, FROZEN_COLS]
    for col in FROZEN_COLS:
        adata.obs[col] = pd.Categorical(aligned[col].to_numpy())
    adata.obs["manual_annotation_source"] = "frozen_barcode_transfer"
    adata.obs["manual_annotation_version"] = args.annotation_version

    summary = pd.DataFrame(
        {
            "metric": [
                "n_fresh", "n_frozen", "n_intersection",
                "fresh_only_dropped", "frozen_only_missing", "frozen_coverage",
            ],
            "value": [
                n_fresh, n_frozen, n_inter, fresh_only, frozen_only, round(coverage, 6),
            ],
        }
    )
    summary.to_csv(out / "freeze_manual_annotation_summary.csv", index=False)
    adata.obs[FROZEN_COLS].to_csv(out / "frozen_labels_per_cell.csv")
    adata.write_h5ad(out / "adata_progenitor_annotated.h5ad")
    print(f"[freeze] wrote adata_progenitor_annotated.h5ad with {adata.n_obs} cells")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Fresh clustered h5ad")
    parser.add_argument("--frozen", required=True, help="Frozen per-cell annotation CSV(.gz)")
    parser.add_argument("--out", required=True)
    parser.add_argument("--annotation-version", default="manual_full_r1.0_2026-05-05")
    parser.add_argument(
        "--min-coverage", type=float, default=0.95,
        help="Abort if the fraction of canonical cells found in the fresh object is below this.",
    )
    main(parser.parse_args())
