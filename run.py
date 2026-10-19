import os
import json
import argparse
import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser(description="Random clustering baseline")
parser.add_argument("--output_dir", "-o", type=str, required=True)
parser.add_argument("--name", "-n", type=str, required=True)
parser.add_argument("--data.h5ad", dest="data_h5ad", type=str, required=True)
parser.add_argument(
    "--data.clusters_truth_num",
    dest="clusters_truth_num_path",
    type=str,
    required=True,
)
parser.add_argument("--method_name", type=str)
parser.add_argument("--seed", type=int, default=42)
args, _ = parser.parse_known_args()

adata = sc.read_h5ad(args.data_h5ad)
with open(args.clusters_truth_num_path) as f:
    k = int(f.read())

rng = np.random.default_rng(args.seed)

# clusters: metrics reads the columns "leiden" and "louvain" by name
clusters = pd.DataFrame(
    {
        "cell_id": adata.obs_names,
        "louvain": rng.integers(0, k, size=adata.n_obs),
        "leiden": rng.integers(0, k, size=adata.n_obs),
    }
)
clusters.to_csv(
    os.path.join(args.output_dir, f"{args.name}.clusters.tsv"), sep="\t", index=False
)

# pca: cell_id + PC* columns (random, since the baseline has no real embedding)
pcs = rng.normal(size=(adata.n_obs, 30))
pca_df = pd.DataFrame(
    pcs, index=adata.obs_names, columns=[f"PC{i + 1}" for i in range(30)]
)
pca_df.to_csv(
    os.path.join(args.output_dir, f"{args.name}.pca.tsv"),
    sep="\t",
    index=True,
    index_label="cell_id",
)

# timings and clustering_info: any JSON dict is accepted by metrics
with open(os.path.join(args.output_dir, f"{args.name}.timings.json"), "w") as f:
    json.dump({"random_labels": 0.0}, f, indent=2)

with open(os.path.join(args.output_dir, f"{args.name}.clustering_info.json"), "w") as f:
    json.dump({"seed": args.seed, "k": k}, f, indent=2)
