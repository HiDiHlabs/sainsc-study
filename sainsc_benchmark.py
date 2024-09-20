#!/usr/bin/env python


def main():
    from argparse import ArgumentParser
    from pathlib import Path

    parser = ArgumentParser(description="")

    parser.add_argument("input", help="Path to input file.", type=Path)
    parser.add_argument("signatures", help="Path to signature file (tsv).", type=Path)
    parser.add_argument(
        "--genes", help="Path to gene file.", type=Path, default=None, required=False
    )
    parser.add_argument(
        "--n_threads",
        help="Number of threads for parallelization.",
        type=int,
        default=1,
        required=False,
    )

    args = parser.parse_args()

    print(args)

    import pandas as pd

    from sainsc import read_StereoSeq

    stereo = read_StereoSeq(args.input, resolution=500, n_threads=args.n_threads)

    stereo.gaussian_kernel(8)

    print(stereo)

    signatures = pd.read_csv(args.signatures, index_col=0, sep="\t").loc[
        lambda df: df.index.isin(stereo.genes)
    ]

    if args.genes is not None:
        signatures = signatures.loc[
            lambda df: df.index.isin(pd.read_csv(args.genes, header=None).iloc[:, 0])
        ]
    print(signatures.shape)

    stereo.assign_celltype(signatures, log=True)

    print("Done")


if __name__ == "__main__":
    main()
