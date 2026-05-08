#!/usr/bin/env python
"""Reconcile HSPC-focused and mature-lineage reference annotations."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


UNKNOWN_VALUES = {"", "nan", "none", "unknown", "unassigned", "na"}
ANNOTATION_COLUMNS = {
    "_predict_cells",
    "popv_prediction",
    "popv_prediction_score",
    "popv_majority_vote_prediction",
    "popv_majority_vote_score",
    "popv_scanvi_agree",
    "scanvi_label",
    "scanvi_confidence",
}


def is_unknown(value: object) -> bool:
    return str(value).strip().lower() in UNKNOWN_VALUES


def first_existing(row: pd.Series, columns: list[str]) -> object:
    for col in columns:
        if col in row.index and not is_unknown(row[col]):
            return row[col]
    return "unknown"


def numeric_column(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index)
    return pd.to_numeric(frame[column], errors="coerce")


def consensus_label(obs: pd.DataFrame, prefix: str) -> pd.Series:
    popv_col = f"{prefix}_popv_prediction"
    majority_col = f"{prefix}_popv_majority_vote_prediction"
    scanvi_col = f"{prefix}_scanvi_label"
    labels = []
    for _, row in obs.iterrows():
        popv = row[popv_col] if popv_col in row.index else "unknown"
        scanvi = row[scanvi_col] if scanvi_col in row.index else "unknown"
        if not is_unknown(popv) and str(popv) == str(scanvi):
            labels.append(popv)
        else:
            labels.append(first_existing(row, [majority_col, popv_col, scanvi_col]))
    return pd.Series(labels, index=obs.index, dtype="object")


def confidence_score(obs: pd.DataFrame, prefix: str) -> pd.Series:
    candidates = [
        numeric_column(obs, f"{prefix}_popv_prediction_score"),
        numeric_column(obs, f"{prefix}_popv_majority_vote_score"),
        numeric_column(obs, f"{prefix}_scanvi_confidence"),
    ]
    return pd.concat(candidates, axis=1).max(axis=1, skipna=True)


def prefixed_annotation_obs(adata: ad.AnnData, prefix: str) -> pd.DataFrame:
    cols = [
        col
        for col in adata.obs.columns
        if col in ANNOTATION_COLUMNS or col.startswith("popv_") or col.startswith("scanvi_")
    ]
    return adata.obs[cols].rename(columns={col: f"{prefix}_{col}" for col in cols}).copy()


def reconcile_row(row: pd.Series, hspc_regex: re.Pattern[str], min_confidence: float) -> dict[str, object]:
    hspc_label = row["hspc_consensus_label"]
    mature_label = row["mature_consensus_label"]
    hspc_conf = row["hspc_confidence"]
    mature_conf = row["mature_confidence"]

    mature_valid = not is_unknown(mature_label)
    hspc_valid = not is_unknown(hspc_label)
    mature_high_conf = mature_valid and pd.notna(mature_conf) and mature_conf >= min_confidence
    hspc_high_conf = hspc_valid and pd.notna(hspc_conf) and hspc_conf >= min_confidence
    mature_is_hspc_like = mature_valid and bool(hspc_regex.search(str(mature_label)))

    if mature_high_conf and not mature_is_hspc_like:
        return {
            "final_cell_type": mature_label,
            "final_annotation_source": "mature_reference",
            "final_annotation_reason": "mature_non_hspc_high_confidence",
            "final_annotation_confidence": mature_conf,
        }
    if hspc_valid and (mature_is_hspc_like or hspc_high_conf):
        reason = "hspc_refines_progenitor_like_state" if mature_is_hspc_like else "hspc_high_confidence"
        return {
            "final_cell_type": hspc_label,
            "final_annotation_source": "hspc_reference",
            "final_annotation_reason": reason,
            "final_annotation_confidence": hspc_conf,
        }
    if mature_valid:
        return {
            "final_cell_type": mature_label,
            "final_annotation_source": "mature_reference",
            "final_annotation_reason": "mature_fallback",
            "final_annotation_confidence": mature_conf,
        }
    if hspc_valid:
        return {
            "final_cell_type": hspc_label,
            "final_annotation_source": "hspc_reference",
            "final_annotation_reason": "hspc_fallback",
            "final_annotation_confidence": hspc_conf,
        }
    return {
        "final_cell_type": "unknown",
        "final_annotation_source": "unresolved",
        "final_annotation_reason": "no_valid_reference_label",
        "final_annotation_confidence": np.nan,
    }


def write_cluster_summary(
    adata: ad.AnnData,
    out: Path,
    cluster_key: str,
    cluster_majority_min_fraction: float,
) -> None:
    if cluster_key not in adata.obs:
        return

    rows = []
    for cluster, obs in adata.obs.groupby(cluster_key, observed=True):
        label_counts = obs["final_cell_type"].astype(str).value_counts()
        source_counts = obs["final_annotation_source"].astype(str).value_counts()
        top_label = label_counts.index[0]
        top_fraction = float(label_counts.iloc[0] / label_counts.sum())
        rows.append(
            {
                cluster_key: cluster,
                "n_cells": int(obs.shape[0]),
                "final_cell_type_top": top_label,
                "final_cell_type_top_fraction": top_fraction,
                "final_annotation_source_top": source_counts.index[0],
                "final_annotation_source_top_fraction": float(source_counts.iloc[0] / source_counts.sum()),
                "review_flag": bool(top_fraction < cluster_majority_min_fraction),
            }
        )
    summary = pd.DataFrame(rows).sort_values(cluster_key)
    summary.to_csv(out / f"reconciled_annotation_summary_{cluster_key}.csv", index=False)
    adata.obs = adata.obs.join(
        summary.set_index(cluster_key)[["final_cell_type_top", "final_cell_type_top_fraction", "review_flag"]],
        on=cluster_key,
    )
    adata.obs = adata.obs.rename(
        columns={
            "final_cell_type_top": "cluster_reconciled_cell_type",
            "final_cell_type_top_fraction": "cluster_reconciled_fraction",
            "review_flag": "cluster_annotation_review_flag",
        }
    )


def plot_reconciled_umaps(adata: ad.AnnData, out: Path) -> None:
    if "X_umap" not in adata.obsm:
        return
    for col in ["final_cell_type", "final_annotation_source", "final_annotation_confidence"]:
        sc.pl.umap(
            adata,
            color=col,
            show=False,
            legend_loc="on data" if col != "final_annotation_confidence" else None,
        )
        plt.savefig(out / f"umap_{col}.png", dpi=150, bbox_inches="tight")
        plt.close()


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    hspc = sc.read_h5ad(args.hspc)
    mature = sc.read_h5ad(args.mature)
    if not hspc.obs_names.equals(mature.obs_names):
        raise ValueError("HSPC and mature annotation files must contain the same cells in the same order.")

    adata = hspc.copy()
    hspc_obs = prefixed_annotation_obs(hspc, "hspc")
    mature_obs = prefixed_annotation_obs(mature, "mature")
    for col in hspc_obs.columns:
        adata.obs[col] = hspc_obs[col].to_numpy()
    for col in mature_obs.columns:
        adata.obs[col] = mature_obs[col].to_numpy()

    adata.obs["hspc_consensus_label"] = consensus_label(adata.obs, "hspc")
    adata.obs["mature_consensus_label"] = consensus_label(adata.obs, "mature")
    adata.obs["hspc_confidence"] = confidence_score(adata.obs, "hspc").to_numpy()
    adata.obs["mature_confidence"] = confidence_score(adata.obs, "mature").to_numpy()

    hspc_regex = re.compile(args.hspc_label_regex)
    reconciled = adata.obs.apply(
        lambda row: reconcile_row(row, hspc_regex=hspc_regex, min_confidence=args.min_confidence),
        axis=1,
        result_type="expand",
    )
    for col in reconciled.columns:
        adata.obs[col] = reconciled[col].to_numpy()

    write_cluster_summary(adata, out, args.cluster_key, args.cluster_majority_min_fraction)
    plot_reconciled_umaps(adata, out)
    adata.obs[
        [
            "hspc_consensus_label",
            "hspc_confidence",
            "mature_consensus_label",
            "mature_confidence",
            "final_cell_type",
            "final_annotation_source",
            "final_annotation_reason",
            "final_annotation_confidence",
        ]
    ].to_csv(out / "reconciled_predictions.csv")
    metadata = {
        "hspc_annotation": args.hspc,
        "mature_annotation": args.mature,
        "hspc_label_regex": args.hspc_label_regex,
        "min_confidence": args.min_confidence,
        "cluster_majority_min_fraction": args.cluster_majority_min_fraction,
    }
    (out / "reconciliation_metadata.json").write_text(json.dumps(metadata, indent=2))
    adata.write_h5ad(out / "adata_reconciled.h5ad")
    print(f"Saved: {out / 'adata_reconciled.h5ad'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hspc", required=True)
    parser.add_argument("--mature", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cluster-key", default="leiden_r0.5")
    parser.add_argument(
        "--hspc-label-regex",
        default=r"(?i)(^|[^A-Za-z])(LT-?HSC|ST-?HSC|HSC|MPP|LMPP|CMP|GMP|MEP|HSPC|progenitor|prog|stem)([^A-Za-z]|$)",
    )
    parser.add_argument("--min-confidence", type=float, default=0.55)
    parser.add_argument("--cluster-majority-min-fraction", type=float, default=0.70)
    main(parser.parse_args())
