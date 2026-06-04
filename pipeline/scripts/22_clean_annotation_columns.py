#!/usr/bin/env python3
"""
22_clean_annotation_columns.py

Slim down the final progenitor AnnData (output of 21z) to a single, unambiguous
set of annotation columns. The object had accumulated 120 obs columns across
several annotation rounds (popV/scANVI global, hspc_*, mature_*, cluster
reconciliation, an old `final_cell_type` round) whose names collided with the
final HemaScribe round — e.g. `final_annotation_source/reason/confidence`
belonged to the *legacy* `final_cell_type`, NOT to `final_annotation`.

This script:
  - renames the canonical final call  final_annotation -> cell_type
                                       final_annotation_rule -> cell_type_rule
  - keeps sample metadata, QC metrics, leiden clusterings, the manual calls,
    the HemaScribe outputs, and the QC flags
  - drops every legacy/intermediate annotation column and constant bookkeeping
  - records the dropped columns + backup path in adata.uns['column_cleanup']
  - rewrites the object IN PLACE (a dated backup must already exist)

Usage:
  pixi run python pipeline/scripts/22_clean_annotation_columns.py \
    --input results/14_progenitor_annotated/adata_hemascribe.h5ad \
    --backup results/14_progenitor_annotated/adata_hemascribe.PRE_CLEANUP_20260604.h5ad
"""
import argparse
import os

import scanpy as sc

# ── Canonical rename: the ONE final call and its provenance companion ────────
RENAME = {
    "final_annotation": "cell_type",
    "final_annotation_rule": "cell_type_rule",
}

# ── Columns to KEEP (after rename is applied, so use the NEW names) ───────────
KEEP = [
    # sample metadata
    "donor", "treatment", "population", "replicate", "sample_id", "pool",
    "genotype",
    # QC metrics
    "n_genes_by_counts", "log1p_n_genes_by_counts",
    "total_counts", "log1p_total_counts",
    "total_counts_mt", "log1p_total_counts_mt", "pct_counts_mt",
    "total_counts_ribo", "log1p_total_counts_ribo", "pct_counts_ribo",
    "total_counts_hb", "log1p_total_counts_hb", "pct_counts_hb",
    "doubletfinder_score",
    # clustering
    "leiden_r0.3", "leiden_r0.5", "leiden_r0.8", "leiden_r1.0",
    "leiden_r1.5", "leiden_r2.0", "leiden_r2.5", "leiden_r3.0",
    # manual annotation (the marker-based human calls)
    "manual_level1", "manual_level2", "manual_annotation_cluster",
    "manual_annotation_basis", "manual_annotation_confidence",
    "manual_review_flag",
    # HemaScribe round
    "hemascribe_broad", "hemascribe_fine", "hemascribe_score",
    "hemascribe_broad_sortgated", "hemascribe_fine_sortgated",
    "hemascribe_sortgate_flag",
    # canonical final + QC flags
    "cell_type", "cell_type_rule", "qc_sort_leak", "qc_low_score",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--backup", default=None,
                   help="Pre-cleanup backup; required only for an in-place rewrite "
                        "(--out omitted or equal to --input). Not needed when --out "
                        "is a distinct path (e.g. under Nextflow).")
    p.add_argument("--out", default=None, help="default: overwrite --input")
    return p.parse_args()


def main():
    args = parse_args()
    out_path = args.out or args.input
    in_place = os.path.abspath(out_path) == os.path.abspath(args.input)

    if in_place:
        # Overwriting the input — a backup MUST exist (manual-workflow safety).
        if not args.backup or not os.path.exists(args.backup):
            raise SystemExit(
                "[22] Refusing in-place rewrite without an existing --backup. "
                "Pass a distinct --out (e.g. under Nextflow) or make a backup first."
            )
    elif args.backup and not os.path.exists(args.backup):
        raise SystemExit(f"[22] --backup {args.backup} does not exist.")

    print(f"[22] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    n_before = adata.obs.shape[1]
    print(f"      obs columns before: {n_before}")

    # rename canonical columns
    present_rename = {k: v for k, v in RENAME.items() if k in adata.obs.columns}
    missing_rename = set(RENAME) - set(present_rename)
    if missing_rename:
        raise SystemExit(f"[22] Expected columns to rename are missing: {missing_rename}")
    adata.obs = adata.obs.rename(columns=present_rename)

    # verify all KEEP columns exist
    missing_keep = [c for c in KEEP if c not in adata.obs.columns]
    if missing_keep:
        raise SystemExit(f"[22] KEEP columns not found (aborting): {missing_keep}")

    dropped = [c for c in adata.obs.columns if c not in KEEP]
    print(f"[22] Keeping {len(KEEP)} columns, dropping {len(dropped)}")

    # record provenance before dropping
    adata.uns["column_cleanup"] = {
        "script": "pipeline/scripts/22_clean_annotation_columns.py",
        "date": "2026-06-04",
        "canonical_annotation": "cell_type",
        "renamed": present_rename,
        "backup_with_all_columns": (
            os.path.basename(args.backup) if args.backup else "none (out != input)"
        ),
        "dropped_columns": sorted(dropped),
        "note": (
            "Legacy popV/scANVI and an old `final_cell_type` round were removed. "
            "`final_annotation_source/reason/confidence` belonged to that legacy "
            "round, NOT to the final call. The canonical annotation is `cell_type` "
            "(HemaScribe + manual validation, script 21z)."
        ),
    }

    adata.obs = adata.obs[KEEP].copy()

    print(f"[22] obs columns after: {adata.obs.shape[1]}")
    print(f"[22] Writing {out_path}")
    adata.write_h5ad(out_path)
    print("[22] Done.")
    print("\n[22] Dropped columns:")
    for c in sorted(dropped):
        print(f"      - {c}")


if __name__ == "__main__":
    main()
