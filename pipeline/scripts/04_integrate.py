#!/usr/bin/env python
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import scvi


def main(in_path: str, out_dir: str, n_latent: int, n_layers: int, max_epochs: int) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    scvi.settings.seed = 777
    adata = sc.read_h5ad(in_path)

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
        gene_likelihood="nb",
    )
    model.train(max_epochs=max_epochs, early_stopping=True, early_stopping_patience=20, batch_size=256)

    fig, ax = plt.subplots()
    for key, label in [("elbo_train", "train ELBO"), ("elbo_validation", "validation ELBO")]:
        if key in model.history:
            ax.plot(model.history[key], label=label)
    ax.set_xlabel("Epoch")
    ax.legend()
    plt.savefig(out / "training_loss.png", dpi=150, bbox_inches="tight")
    plt.close()

    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30, metric="euclidean")
    sc.tl.umap(adata, min_dist=0.3)

    for color in ["sample_id", "pool", "donor", "treatment", "population"]:
        sc.pl.umap(adata, color=color, show=False)
        plt.savefig(out / f"umap_{color}.png", dpi=150, bbox_inches="tight")
        plt.close()

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
