#!/usr/bin/env python
"""Convert the Nestorowa Matrix Market reference export to H5AD."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import io


def main(export_dir: str, out_path: str) -> None:
    export = Path(export_dir)
    out = Path(out_path)
    counts_path = export / "counts.genes_by_cells.mtx"
    genes_path = export / "genes.tsv"
    obs_path = export / "obs.csv"
    for path in [counts_path, genes_path, obs_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    counts = io.mmread(counts_path).tocsr()
    genes = pd.read_csv(genes_path, sep="\t")
    obs = pd.read_csv(obs_path, index_col=0)

    var_names = genes["gene_symbol"].astype(str).to_numpy()
    obs.index = obs.index.astype(str)
    adata = ad.AnnData(X=counts.T.tocsr(), obs=obs, var=pd.DataFrame(index=var_names))
    adata.var["gene_id"] = genes["gene_id"].astype(str).to_numpy()
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    if "cell_type" not in adata.obs:
        raise ValueError("Reference export is missing obs['cell_type'].")
    if "reference_batch" not in adata.obs:
        adata.obs["reference_batch"] = "Nestorowa2016"

    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"Saved {out}: {adata.n_obs} cells x {adata.n_vars} genes")
    print(adata.obs["cell_type"].value_counts())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    main(args.export_dir, args.out)
