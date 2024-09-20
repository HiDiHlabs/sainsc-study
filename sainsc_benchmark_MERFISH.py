#!/usr/bin/env python


def main():
    from argparse import ArgumentParser
    from pathlib import Path

    parser = ArgumentParser(description="")

    parser.add_argument("input", help="Path to input file.", type=Path)
    parser.add_argument("signatures", help="Path to signature file (tsv).", type=Path)
    parser.add_argument(
        "--n_threads", help="Number of threads for parallelization.", type=int
    )

    args = parser.parse_args()
    print(args)

    import os

    os.environ["POLARS_MAX_THREADS"] = str(args.n_threads)

    import pandas as pd
    import polars as pl

    from sainsc import LazyKDE

    merfish = LazyKDE.from_dataframe(
        pl.scan_csv(
            args.input,
            schema_overrides={"target_gene": pl.Categorical},
        )
        .select(["target_gene", "global_x", "global_y"])
        .rename({"target_gene": "gene", "global_x": "x", "global_y": "y"})
        .filter(~pl.col("gene").cast(pl.Utf8).str.contains("blank"))
        .collect(),
        binsize=0.5,
        resolution=1_000,
        n_threads=args.n_threads,
    )

    merfish.gaussian_kernel(8)
    print(merfish)

    signatures = pd.read_csv(args.signatures, index_col=0, sep="\t").loc[
        lambda df: df.index.isin(merfish.genes)
    ]
    print(signatures.shape)

    merfish.assign_celltype(signatures)

    print("Done")


if __name__ == "__main__":
    main()
