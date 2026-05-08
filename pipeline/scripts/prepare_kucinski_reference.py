#!/usr/bin/env python
"""Prepare the Kucinski 2024 murine bone marrow h5ad as an annotation reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


LABEL_KEY_CANDIDATES = [
    "cell_type",
    "celltype",
    "CellType",
    "annotation",
    "Annotation",
    "cell_state",
    "cell_label",
    "state",
    "population",
    "cluster_label",
    "celltype_major",
    "celltype_l2",
    "celltype.l2",
]
BATCH_KEY_CANDIDATES = [
    "reference_batch",
    "sample_id",
    "sample",
    "batch",
    "library",
    "timepoint",
    "orig.ident",
]
UNKNOWN_VALUES = {"", "nan", "none", "unknown", "unassigned", "na"}


def pick_obs_key(obs: pd.DataFrame, requested: str | None, candidates: list[str], role: str) -> str:
    if requested:
        if requested not in obs:
            raise ValueError(f"Requested {role} key {requested!r} is not present in obs.")
        return requested
    for key in candidates:
        if key in obs:
            return key
    raise ValueError(
        f"Could not auto-detect {role} key. Available obs columns: {', '.join(map(str, obs.columns))}"
    )


def normalize_labels(values: pd.Series) -> pd.Series:
    labels = values.astype(str).str.strip()
    labels = labels.mask(labels.str.lower().isin(UNKNOWN_VALUES), "unknown")
    return labels


def maybe_use_gene_symbols(adata: sc.AnnData, gene_symbol_key: str | None) -> tuple[sc.AnnData, str | None]:
    if not gene_symbol_key:
        adata.var_names_make_unique()
        return adata, None
    if gene_symbol_key not in adata.var:
        raise ValueError(f"Requested gene symbol key {gene_symbol_key!r} is not present in var.")
    symbols = adata.var[gene_symbol_key].astype(str).str.strip()
    valid = ~symbols.str.lower().isin(UNKNOWN_VALUES)
    if valid.sum() < 500:
        raise ValueError(f"Gene symbol key {gene_symbol_key!r} has too few usable symbols.")
    adata = adata[:, valid.to_numpy()].copy()
    adata.var_names = symbols[valid].to_numpy()
    adata.var_names_make_unique()
    return adata, gene_symbol_key


def recover_counts_from_log1p_cpm(adata: sc.AnnData, counts_per_cell_after: float = 10000.0) -> str:
    if "n_counts" not in adata.obs:
        raise ValueError("Cannot recover counts from log1p CPM values because obs['n_counts'] is missing.")
    ratios = pd.to_numeric(adata.obs["n_counts"], errors="coerce").to_numpy() / counts_per_cell_after
    if np.isnan(ratios).any():
        raise ValueError("Cannot recover counts because obs['n_counts'] contains non-numeric values.")

    if sparse.issparse(adata.X):
        x = adata.X.copy().tocsr()
        x.data = np.expm1(x.data)
        x = x.multiply(ratios[:, np.newaxis]).tocsr()
        x.data = np.rint(x.data).astype(np.float32)
        x.eliminate_zeros()
        adata.X = x
    else:
        x = np.expm1(np.asarray(adata.X))
        x = ratios[:, np.newaxis] * x
        adata.X = np.rint(x).astype(np.float32)
    adata.layers["counts"] = adata.X.copy()
    return "recovered_from_log1p_cpm_using_obs_n_counts"


def ensure_counts_layer(adata: sc.AnnData, recover_from_log1p_cpm: bool) -> str:
    if "counts" in adata.layers:
        return "existing_counts_layer"
    if recover_from_log1p_cpm:
        return recover_counts_from_log1p_cpm(adata)
    adata.layers["counts"] = adata.X.copy()
    return "copied_from_X"


def filter_and_downsample(
    adata: sc.AnnData,
    label_key: str,
    min_cells_per_label: int,
    max_cells_per_label: int | None,
    seed: int,
) -> sc.AnnData:
    labels = adata.obs[label_key].astype(str)
    label_counts = labels.value_counts()
    keep_labels = label_counts[label_counts >= min_cells_per_label].index
    keep = labels.isin(keep_labels) & (labels != "unknown")
    filtered = adata[keep.to_numpy()].copy()

    if not max_cells_per_label:
        return filtered

    rng = np.random.default_rng(seed)
    selected = []
    for _, obs in filtered.obs.groupby(label_key, observed=True):
        idx = obs.index.to_numpy()
        if len(idx) > max_cells_per_label:
            idx = rng.choice(idx, size=max_cells_per_label, replace=False)
        selected.extend(idx.tolist())
    return filtered[selected].copy()


def main(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    source = sc.read_h5ad(args.input)
    if args.use_raw and source.raw is not None:
        adata = source.raw.to_adata()
    else:
        adata = source

    label_key = pick_obs_key(adata.obs, args.label_key, LABEL_KEY_CANDIDATES, "label")
    batch_key = pick_obs_key(adata.obs, args.batch_key, BATCH_KEY_CANDIDATES, "batch")
    adata, gene_symbol_source = maybe_use_gene_symbols(adata, args.gene_symbol_key)

    adata.obs["cell_type"] = normalize_labels(adata.obs[label_key])
    adata.obs["reference_batch"] = adata.obs[batch_key].astype(str).to_numpy()
    prepared = filter_and_downsample(
        adata,
        label_key="cell_type",
        min_cells_per_label=args.min_cells_per_label,
        max_cells_per_label=args.max_cells_per_label,
        seed=args.seed,
    )
    counts_status = ensure_counts_layer(prepared, args.recover_counts_from_log1p_cpm)

    metadata = {
        "source_h5ad": args.input,
        "output_h5ad": str(out),
        "selected_label_key": label_key,
        "selected_batch_key": batch_key,
        "used_raw": bool(args.use_raw and source.raw is not None),
        "gene_symbol_source": gene_symbol_source,
        "counts_status": counts_status,
        "n_cells_input": int(adata.n_obs),
        "n_genes_input": int(adata.n_vars),
        "n_cells_output": int(prepared.n_obs),
        "n_genes_output": int(prepared.n_vars),
        "min_cells_per_label": args.min_cells_per_label,
        "max_cells_per_label": args.max_cells_per_label,
        "label_counts": prepared.obs["cell_type"].astype(str).value_counts().to_dict(),
    }
    prepared.write_h5ad(out)
    out.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--label-key")
    parser.add_argument("--batch-key")
    parser.add_argument("--gene-symbol-key")
    parser.add_argument("--use-raw", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--recover-counts-from-log1p-cpm", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--min-cells-per-label", type=int, default=20)
    parser.add_argument("--max-cells-per-label", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=777)
    main(parser.parse_args())
