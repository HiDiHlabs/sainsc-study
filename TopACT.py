#!/usr/bin/env python


def main():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description="")

    parser.add_argument("training", help="Path to h5ad file for training.", type=Path)
    parser.add_argument("spatial", help="Path of spatial data (tsv).", type=Path)
    parser.add_argument("out", help="Path of output directory", type=Path)
    parser.add_argument(
        "--n_processes",
        help="Number of processes for classification",
        type=int,
        default=1,
        required=False,
    )

    args = parser.parse_args()

    import shutil
    from tempfile import TemporaryDirectory

    import pandas as pd
    from skimage.morphology import convex_hull_image
    from topact.spatial import CountGrid

    df_spatial = pd.read_table(args.spatial, comment="#", dtype={"geneID": "category"})

    clf, genes = train(args.training, df_spatial["geneID"].cat.categories)

    sd = CountGrid.from_coord_table(
        df_spatial,
        genes=genes,
        count_col="MIDCounts",
        gene_col="geneID",
    )

    mask = convex_hull_image(sd.count_matrix())

    with TemporaryDirectory(suffix="_topact") as dir:
        tmp_file = Path(dir) / "confidence.npy"
        sd.classify_parallel(
            clf,
            min_scale=3,
            max_scale=9,
            num_proc=args.n_processes,
            outfile=str(tmp_file),
            mask=mask,
            verbose=False,
        )

        shutil.move(tmp_file, args.out / "confidence.npy")


def train(h5ad_path, genes, label_col="subclass_label"):
    import anndata as ad
    from sklearn.preprocessing import normalize
    from topact.classifier import SVCClassifier, train_from_countmatrix
    from topact.countdata import CountMatrix

    adata = ad.read_h5ad(h5ad_path)
    adata = adata[:, adata.var_names.isin(genes)]

    mtx = normalize(adata.X, norm="l1")
    genes = adata.var_names.tolist()
    celltypes = adata.obs[label_col]

    sc = CountMatrix(mtx, genes=genes)
    sc.add_metadata("celltype", celltypes)

    clf = SVCClassifier()
    train_from_countmatrix(clf, sc, "celltype")

    return clf, adata.var_names.tolist()


if __name__ == "__main__":
    main()