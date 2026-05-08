#!/usr/bin/env python
"""Apply human-reviewed progenitor annotation map to final AnnData.

Reads pipeline/config/progenitor_annotation_map.tsv (filled by human after
reviewing results/13_progenitor_recluster/ outputs). Updates manual_level1
and manual_level2 for progenitor cells in adata_annotated_final.h5ad.
Saves the result to results/14_progenitor_annotated/adata_progenitor_annotated.h5ad.

Usage:
    pixi run python pipeline/scripts/14_apply_progenitor_annotation.py \\
        --input results/12_subset_recluster/adata_annotated_final.h5ad \\
        --assignments results/13_progenitor_recluster/leiden_subset_assignments.csv \\
        --map pipeline/config/progenitor_annotation_map.tsv \\
        --out results/14_progenitor_annotated
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1


def read_chosen_resolution(map_path: str) -> str:
    """Parse the # chosen_resolution: X.X comment from the TSV header.

    Returns the raw string (e.g. "0.4") rather than a float, to avoid
    float-string round-trip mismatches when looking up
    `leiden_sub_r{chosen_res}` columns or filtering the map TSV's
    `leiden_resolution` column (where "0.4" and "0.40" must match).
    """
    with open(map_path) as f:
        for line in f:
            m = re.match(r"#\s*chosen_resolution:\s*([0-9.]+)", line.strip())
            if m:
                return m.group(1)
    raise ValueError(
        f"No '# chosen_resolution: X.X' comment found in {map_path}.\n"
        "Add a comment line at the top of the file, e.g.:\n"
        "  # chosen_resolution: 0.4"
    )


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # ── Load ─────────────────────────────────────────────────────────────────
    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells")

    assignments = pd.read_csv(args.assignments, index_col=0)
    print(f"Leiden assignments loaded: {assignments.shape[0]} progenitor barcodes")

    chosen_res = read_chosen_resolution(args.map)
    leiden_col = f"leiden_sub_r{chosen_res}"
    if leiden_col not in assignments.columns:
        available = [c for c in assignments.columns if c.startswith("leiden_sub")]
        raise ValueError(
            f"Column '{leiden_col}' not in assignments CSV.\n"
            f"Available: {available}\n"
            "Update # chosen_resolution comment in the map TSV."
        )
    print(f"Using leiden resolution: {chosen_res} (column: {leiden_col})")

    # ── Load map ─────────────────────────────────────────────────────────────
    map_df = pd.read_csv(args.map, sep="\t", comment="#")
    # Compare resolutions as floats so "0.4" and "0.40" both match.
    map_df = map_df[
        pd.to_numeric(map_df["leiden_resolution"], errors="coerce") == float(chosen_res)
    ].copy()
    if map_df.empty:
        raise ValueError(
            f"No rows in {args.map} match leiden_resolution={chosen_res}. "
            "Check the # chosen_resolution comment and the leiden_resolution column."
        )
    map_df["subcluster_id"] = map_df["subcluster_id"].astype(str)
    required = {"leiden_resolution", "subcluster_id", "new_level1", "new_level2"}
    missing = required - set(map_df.columns)
    if missing:
        raise ValueError(f"Missing columns in annotation map: {missing}")
    empty = map_df[map_df["new_level1"].isna() | (map_df["new_level1"] == "")]
    if not empty.empty:
        raise ValueError(
            f"Empty new_level1 for {len(empty)} rows. Fill all rows before running."
        )

    subcluster_to_level1 = dict(zip(map_df["subcluster_id"], map_df["new_level1"]))
    # Drop rows where new_level2 is missing/empty so we don't silently overwrite
    # a previously valid manual_level2 with NaN/"".
    level2_valid = map_df.loc[
        map_df["new_level2"].notna() & (map_df["new_level2"].astype(str) != "")
    ]
    subcluster_to_level2 = dict(zip(level2_valid["subcluster_id"], level2_valid["new_level2"]))
    print(f"Annotation map loaded: {len(subcluster_to_level1)} cluster mappings")
    for k, v in subcluster_to_level1.items():
        print(f"  subcluster {k} → level1={v}, level2={subcluster_to_level2.get(k)}")

    # ── Apply ─────────────────────────────────────────────────────────────────
    adata.obs["manual_level1"] = adata.obs["manual_level1"].astype(str)
    adata.obs["manual_level2"] = adata.obs["manual_level2"].astype(str)

    assignments["subcluster"] = assignments[leiden_col].astype(str)
    n_updated_l1 = 0
    n_updated_l2 = 0
    for barcode, row in assignments.iterrows():
        sub = row["subcluster"]
        if barcode not in adata.obs_names:
            continue
        if sub in subcluster_to_level1:
            adata.obs.at[barcode, "manual_level1"] = subcluster_to_level1[sub]
            n_updated_l1 += 1
        if sub in subcluster_to_level2:
            adata.obs.at[barcode, "manual_level2"] = subcluster_to_level2[sub]
            n_updated_l2 += 1

    adata.obs["manual_level1"] = adata.obs["manual_level1"].astype("category")
    adata.obs["manual_level2"] = adata.obs["manual_level2"].astype("category")
    print(f"\nUpdated {n_updated_l1} barcodes for level1, {n_updated_l2} for level2")

    # ── Summary ───────────────────────────────────────────────────────────────
    counts = (
        adata.obs.groupby(["manual_level1", "manual_level2"], observed=True)
        .size().reset_index(name="n_cells")
    )
    counts["pct"] = (counts["n_cells"] / adata.n_obs * 100).round(2)
    counts.to_csv(out / "final_annotation_counts.csv", index=False)
    print("\nFinal annotation counts:")
    print(counts.to_string(index=False))

    # ── UMAP verification plot ────────────────────────────────────────────────
    sc.pl.umap(
        adata, color=["manual_level1", "manual_level2"],
        show=False, frameon=False, legend_loc="right margin",
        legend_fontsize=6, ncols=2,
    )
    plt.savefig(out / "umap_final_annotation.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"✓ UMAP saved to {out}/umap_final_annotation.png")

    # ── Save ──────────────────────────────────────────────────────────────────
    out_h5ad = out / "adata_progenitor_annotated.h5ad"
    adata.write_h5ad(out_h5ad)
    print(f"✓ Saved {out_h5ad}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",
        default="results/12_subset_recluster/adata_annotated_final.h5ad")
    parser.add_argument("--assignments",
        default="results/13_progenitor_recluster/leiden_subset_assignments.csv")
    parser.add_argument("--map",
        default="pipeline/config/progenitor_annotation_map.tsv")
    parser.add_argument("--out",
        default="results/14_progenitor_annotated")
    main(parser.parse_args())
