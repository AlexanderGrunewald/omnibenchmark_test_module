import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--output_dir", "-o", required=True)
parser.add_argument("--name", "-n", required=True)
parser.add_argument("--data.h5ad", dest="data_h5ad", required=True)
parser.add_argument("--data.clusters_truth_num", dest="truth_num_path", required=True)
parser.add_argument("--seed", type=int, default=42)
args, _ = parser.parse_known_args()

adata = sc.read_h5ad(args.data_h5ad)            
with open(args.truth_num_path) as f:
    k = int(f.read())                          

rng = np.random.default_rng(args.seed)
labels = rng.integers(0, k, size=adata.n_obs)

pd.DataFrame({"cell_id": adata.obs_names, "random": labels}).to_csv(
    os.path.join(args.output_dir, f"{args.name}.clusters.tsv"),
    sep="\t", index=False,
)
