#!/usr/bin/env python


def main():
    from argparse import ArgumentParser
    from pathlib import Path

    parser = ArgumentParser(description="")

    parser.add_argument("input", help="Path to input file.", type=Path)
    parser.add_argument("signatures", help="Path to signature file (tsv).", type=Path)
    parser.add_argument("out", help="Output file (pickle).", type=Path)
    parser.add_argument(
        "--n_threads", help="Number of threads for parallelization.", type=int
    )
    parser.add_argument(
        "--truncate", help="Truncation of kernel", type=float, default=2
    )

    args = parser.parse_args()
    print(args)

    import pickle

    import pandas as pd
    import polars as pl

    from sainsc import LazyKDE

    xenium = LazyKDE.from_dataframe(
        pl.read_csv(
            args.input,
            columns=["feature_name", "x_location", "y_location"],
            schema_overrides={"feature_name": pl.Categorical},
            n_threads=args.n_threads,
        )
        .rename({"feature_name": "gene", "x_location": "x", "y_location": "y"})
        .filter(~pl.col("gene").cast(pl.Utf8).str.contains("(BLANK|NegControl)")),
        binsize=1,
        resolution=1_000,
        n_threads=args.n_threads,
    )

    xenium.gaussian_kernel(2.5, unit="um", truncate=args.truncate)
    xenium.calculate_total_mRNA_KDE()

    signatures = pd.read_table(args.signatures, index_col=0).loc[
        :, lambda df: ~df.columns.str.startswith(("Missegmented", "Noise"))
    ]

    xenium.assign_celltype(signatures)
    print(xenium)

    with open(args.out, "wb") as file:
        pickle.dump(xenium, file)

    print("Done")


if __name__ == "__main__":
    main()
