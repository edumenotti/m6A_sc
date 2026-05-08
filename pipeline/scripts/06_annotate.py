#!/usr/bin/env python
"""Reference-based annotation with popV primary and scANVI secondary.

Primary annotation uses popV with a labelled mouse HSPC reference, excluding
CellTypist by default. Secondary annotation trains a scVI/scANVI model on the
same reference/query gene intersection and predicts query labels independently.
"""

from __future__ import annotations

import argparse
import inspect
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from scipy import sparse


DEFAULT_POPV_METHODS = [
    "KNN_BBKNN",
    "KNN_SCVI",
    "SCANVI_POPV",
    "Support_Vector",
    "XGboost",
]
UNKNOWN_LABEL = "unknown"


def configure_torch_runtime() -> None:
    """Compatibility/performance setup for scvi-tools + PyTorch >= 2.6."""
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("high")

    original_load = torch.load

    def torch_load_scvi_compat(*args, **kwargs):
        kwargs.setdefault("weights_only", False)
        return original_load(*args, **kwargs)

    torch.load = torch_load_scvi_compat


def install_scanpy_pca_layer_compat() -> None:
    """Allow popV versions that call scanpy.pp.pca(layer=...) on older Scanpy."""
    if "layer" in inspect.signature(sc.pp.pca).parameters:
        return
    original_pca = sc.pp.pca

    def pca_with_layer_compat(adata, *args, layer=None, **kwargs):
        if layer is None:
            return original_pca(adata, *args, **kwargs)
        if layer not in adata.layers:
            raise KeyError(f"Layer {layer!r} not found for scanpy.pp.pca compatibility shim.")
        original_x = adata.X
        try:
            adata.X = adata.layers[layer]
            return original_pca(adata, *args, **kwargs)
        finally:
            adata.X = original_x

    sc.pp.pca = pca_with_layer_compat


def parse_methods(value: str) -> list[str]:
    if value.lower() in {"none", "skip", ""}:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def ensure_counts_layer(adata: ad.AnnData, layer: str = "counts") -> None:
    if layer not in adata.layers:
        adata.layers[layer] = adata.X.copy()


def load_reference(path: str, labels_key: str, batch_key: str) -> ad.AnnData:
    ref = sc.read_h5ad(path)
    if labels_key not in ref.obs:
        raise ValueError(f"Reference is missing obs['{labels_key}'].")
    if batch_key not in ref.obs:
        ref.obs[batch_key] = "reference"
    ensure_counts_layer(ref, "counts")
    ref.obs_names = "ref_" + ref.obs_names.astype(str)
    return ref


def copy_query_for_reference_mapping(adata: ad.AnnData) -> ad.AnnData:
    query = adata.copy()
    ensure_counts_layer(query, "counts")
    query.X = query.layers["counts"].copy()
    return query


def run_popv(
    query: ad.AnnData,
    ref: ad.AnnData,
    out: Path,
    labels_key: str,
    ref_batch_key: str,
    query_batch_key: str,
    methods: list[str],
    hvg: int,
    n_samples_per_label: int,
    prediction_mode: str,
    cl_obo_folder: str | bool,
) -> tuple[pd.DataFrame, dict]:
    if not methods:
        return pd.DataFrame(index=query.obs_names), {"status": "skipped", "methods": []}

    install_scanpy_pca_layer_compat()

    from popv.annotation import annotate_data
    from popv.preprocessing import Process_Query

    model_dir = out / "popv_models"
    model_dir.mkdir(parents=True, exist_ok=True)

    processor = Process_Query(
        query_adata=copy_query_for_reference_mapping(query),
        ref_adata=ref.copy(),
        ref_labels_key=labels_key,
        ref_batch_key=ref_batch_key,
        cl_obo_folder=cl_obo_folder,
        query_batch_key=query_batch_key if query_batch_key in query.obs else None,
        query_layer_key="counts",
        ref_layer_key="counts",
        prediction_mode=prediction_mode,
        unknown_celltype_label=UNKNOWN_LABEL,
        n_samples_per_label=n_samples_per_label,
        save_path_trained_models=str(model_dir),
        hvg=hvg,
    )
    popv_adata = processor.adata
    annotate_data(popv_adata, methods=methods, save_path=str(out / "popv"))

    query_obs = popv_adata.obs.loc[popv_adata.obs["_dataset"] == "query"].copy()
    popv_cols = [
        col
        for col in query_obs.columns
        if col.startswith("popv_") or col in {"_predict_cells"}
    ]
    predictions = query_obs[popv_cols].copy()
    predictions = predictions.reindex(query.obs_names)
    metadata = {
        "status": "completed",
        "methods": methods,
        "prediction_mode": prediction_mode,
        "hvg": hvg,
        "n_samples_per_label": n_samples_per_label,
        "n_query_cells_returned": int(predictions.shape[0]),
    }
    return predictions, metadata


def train_scanvi_secondary(
    query: ad.AnnData,
    ref: ad.AnnData,
    out: Path,
    labels_key: str,
    ref_batch_key: str,
    query_batch_key: str,
    n_latent: int,
    scvi_epochs: int,
    scanvi_epochs: int,
    batch_size: int,
) -> tuple[pd.DataFrame, dict]:
    common_genes = ref.var_names.intersection(query.var_names)
    if len(common_genes) < 500:
        raise ValueError(f"Too few common genes for scANVI label transfer: {len(common_genes)}")

    ref_sub = ref[:, common_genes].copy()
    query_sub = copy_query_for_reference_mapping(query[:, common_genes].copy())
    ref_sub.obs["scanvi_label_input"] = ref_sub.obs[labels_key].astype(str).to_numpy()
    query_sub.obs["scanvi_label_input"] = UNKNOWN_LABEL
    ref_sub.obs["scanvi_batch"] = ref_sub.obs[ref_batch_key].astype(str).to_numpy()
    if query_batch_key in query_sub.obs:
        query_sub.obs["scanvi_batch"] = "query_" + query_sub.obs[query_batch_key].astype(str)
    else:
        query_sub.obs["scanvi_batch"] = "query"
    ref_sub.obs["scanvi_dataset"] = "reference"
    query_sub.obs["scanvi_dataset"] = "query"

    combined = ad.concat(
        {"reference": ref_sub, "query": query_sub},
        join="inner",
        label="scanvi_source",
        index_unique=None,
    )
    combined.layers["counts"] = combined.layers["counts"].astype(np.float32)
    scvi.settings.seed = 777
    scvi.model.SCVI.setup_anndata(
        combined,
        layer="counts",
        batch_key="scanvi_batch",
        labels_key="scanvi_label_input",
    )
    scvi_model = scvi.model.SCVI(
        combined,
        n_latent=n_latent,
        n_layers=2,
        gene_likelihood="nb",
    )
    scvi_model.train(
        max_epochs=scvi_epochs,
        early_stopping=True,
        batch_size=batch_size,
        accelerator="auto",
        devices="auto",
    )
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key="scanvi_label_input",
        unlabeled_category=UNKNOWN_LABEL,
    )
    scanvi_model.train(
        max_epochs=scanvi_epochs,
        n_samples_per_label=100,
        early_stopping=True,
        batch_size=batch_size,
        accelerator="auto",
        devices="auto",
    )

    scanvi_dir = out / "scanvi_model"
    scanvi_model.save(str(scanvi_dir), overwrite=True)
    labels = scanvi_model.predict(combined)
    probabilities = scanvi_model.predict(combined, soft=True)
    if isinstance(probabilities, pd.DataFrame):
        confidence = probabilities.max(axis=1).to_numpy()
    else:
        confidence = np.asarray(probabilities).max(axis=1)

    combined.obs["scanvi_label"] = labels
    combined.obs["scanvi_confidence"] = confidence
    combined.obsm["X_scanvi"] = scanvi_model.get_latent_representation()

    query_mask = combined.obs["scanvi_dataset"] == "query"
    predictions = combined.obs.loc[query_mask, ["scanvi_label", "scanvi_confidence"]].copy()
    predictions = predictions.reindex(query.obs_names)
    pd.DataFrame(
        {
            "cell": predictions.index,
            "scanvi_label": predictions["scanvi_label"].astype(str).to_numpy(),
            "scanvi_confidence": predictions["scanvi_confidence"].to_numpy(),
        }
    ).to_csv(out / "scanvi_predictions.csv", index=False)

    metadata = {
        "status": "completed",
        "common_genes": int(len(common_genes)),
        "n_latent": n_latent,
        "scvi_epochs": scvi_epochs,
        "scanvi_epochs": scanvi_epochs,
    }
    return predictions, metadata


def write_cluster_summary(adata: ad.AnnData, out: Path, cluster_key: str) -> None:
    if cluster_key not in adata.obs:
        return
    aggregations = {"n_cells": (cluster_key, "count")}
    for col in ["popv_prediction", "popv_majority_vote_prediction", "scanvi_label"]:
        if col in adata.obs:
            aggregations[f"{col}_top"] = (col, lambda x: x.astype(str).value_counts().index[0])
            aggregations[f"{col}_top_fraction"] = (
                col,
                lambda x: float(x.astype(str).value_counts(normalize=True).iloc[0]),
            )
    summary = adata.obs.groupby(cluster_key, observed=True).agg(**aggregations)
    summary.to_csv(out / f"annotation_summary_{cluster_key}.csv")
    print(summary)


def plot_annotation_umaps(adata: ad.AnnData, out: Path) -> None:
    colors = [
        col
        for col in [
            "popv_prediction",
            "popv_majority_vote_prediction",
            "popv_prediction_score",
            "scanvi_label",
            "scanvi_confidence",
        ]
        if col in adata.obs
    ]
    for col in colors:
        sc.pl.umap(adata, color=col, show=False, legend_loc="on data" if adata.obs[col].dtype.name == "category" else None)
        plt.savefig(out / f"umap_{col}.png", dpi=150, bbox_inches="tight")
        plt.close()


def main(args: argparse.Namespace) -> None:
    configure_torch_runtime()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.input)
    ref = load_reference(args.ref, args.ref_labels_key, args.ref_batch_key)
    ensure_counts_layer(adata, "counts")

    metadata = {
        "reference": args.ref,
        "ref_labels_key": args.ref_labels_key,
        "ref_batch_key": args.ref_batch_key,
        "query_batch_key": args.query_batch_key,
    }

    popv_predictions, popv_metadata = run_popv(
        query=adata,
        ref=ref,
        out=out,
        labels_key=args.ref_labels_key,
        ref_batch_key=args.ref_batch_key,
        query_batch_key=args.query_batch_key,
        methods=parse_methods(args.popv_methods),
        hvg=args.popv_hvg,
        n_samples_per_label=args.popv_samples_per_label,
        prediction_mode=args.popv_mode,
        cl_obo_folder=False if args.cl_obo_folder.lower() == "false" else args.cl_obo_folder,
    )
    metadata["popv"] = popv_metadata
    for col in popv_predictions.columns:
        adata.obs[col] = popv_predictions[col].to_numpy()

    if args.skip_scanvi:
        metadata["scanvi"] = {"status": "skipped"}
    else:
        scanvi_predictions, scanvi_metadata = train_scanvi_secondary(
            query=adata,
            ref=ref,
            out=out,
            labels_key=args.ref_labels_key,
            ref_batch_key=args.ref_batch_key,
            query_batch_key=args.query_batch_key,
            n_latent=args.scanvi_n_latent,
            scvi_epochs=args.scvi_epochs,
            scanvi_epochs=args.scanvi_epochs,
            batch_size=args.batch_size,
        )
        metadata["scanvi"] = scanvi_metadata
        for col in scanvi_predictions.columns:
            adata.obs[col] = scanvi_predictions[col].to_numpy()

    if "popv_prediction" in adata.obs and "scanvi_label" in adata.obs:
        adata.obs["popv_scanvi_agree"] = (
            adata.obs["popv_prediction"].astype(str) == adata.obs["scanvi_label"].astype(str)
        )

    write_cluster_summary(adata, out, args.cluster_key)
    plot_annotation_umaps(adata, out)
    (out / "annotation_metadata.json").write_text(json.dumps(metadata, indent=2))
    adata.write_h5ad(out / "adata_annotated.h5ad")
    print(f"Saved: {out / 'adata_annotated.h5ad'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--ref-labels-key", default="cell_type")
    parser.add_argument("--ref-batch-key", default="reference_batch")
    parser.add_argument("--query-batch-key", default="sample_id")
    parser.add_argument("--cluster-key", default="leiden_r0.5")
    parser.add_argument("--popv-methods", default=",".join(DEFAULT_POPV_METHODS))
    parser.add_argument("--popv-mode", default="retrain", choices=["retrain", "inference", "fast"])
    parser.add_argument("--popv-hvg", type=int, default=4000)
    parser.add_argument("--popv-samples-per-label", type=int, default=300)
    parser.add_argument("--cl-obo-folder", default="false")
    parser.add_argument("--skip-scanvi", action="store_true")
    parser.add_argument("--scanvi-n-latent", type=int, default=30)
    parser.add_argument("--scvi-epochs", type=int, default=200)
    parser.add_argument("--scanvi-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=256)
    main(parser.parse_args())
