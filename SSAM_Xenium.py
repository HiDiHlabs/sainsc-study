#!/usr/bin/env python


def main():
    from argparse import ArgumentParser
    from pathlib import Path

    parser = ArgumentParser(description="")

    parser.add_argument("input", help="Path to input file.", type=Path)
    parser.add_argument("signatures", help="Path to signature file (tsv).", type=Path)
    parser.add_argument("out", help="Path to output file (npz).", type=Path)
    parser.add_argument(
        "--n_threads", help="Number of threads for parallelization.", type=int
    )
    parser.add_argument(
        "--fast_kde", help="Use fast_kde implementation", action="store_true"
    )
    parser.add_argument(
        "--temp_dir",
        help="Path to use as temp_dir, can be None to use system default",
        type=Path,
    )

    args = parser.parse_args()
    print(args)

    import math
    from tempfile import TemporaryDirectory

    import numpy as np
    import pandas as pd
    import ssam

    transcripts = (
        pd.read_csv(
            args.input,
            usecols=["feature_name", "x_location", "y_location"],
            dtype={"feature_name": "category"},
        )
        .rename(columns={"feature_name": "gene", "x_location": "x", "y_location": "y"})
        .loc[lambda df: ~df["gene"].str.startswith(("BLANK", "NegControl"))]
        .assign(feature_name=lambda df: df["gene"].cat.remove_unused_categories())
    )
    # To get the same "bins" (and therefore shape) as sainsc does so we can match pixels
    transcripts[["x", "y"]] -= transcripts[["x", "y"]].min().astype(int)

    width, height = (math.ceil(i) for i in transcripts[["x", "y"]].max())

    signatures = pd.read_table(args.signatures, index_col=0).loc[
        :, lambda df: ~df.columns.str.startswith(("Missegmented", "Noise"))
    ]

    genes, mrna_loci = zip(
        *(
            (gene, df[["x", "y"]].to_numpy())
            for gene, df in transcripts.groupby("gene", observed=True)
            if gene in signatures.index
        )
    )

    del transcripts

    signatures = signatures.loc[genes, :].to_numpy().transpose()
    print(signatures.shape)

    ds = ssam.SSAMDataset(genes, mrna_loci, width, height)

    with TemporaryDirectory(dir=args.temp_dir) as tmp_dir:
        analysis = ssam.SSAMAnalysis(
            ds, ncores=args.n_threads, save_dir=tmp_dir, verbose=True
        )

        if args.fast_kde:
            analysis.run_fast_kde(bandwidth=2.5)
        else:
            analysis.run_kde(bandwidth=2.5)
        print("Finished KDE")

        analysis.map_celltypes(signatures)
        print("Finished cell-type map")

        np.savez(
            args.out, ct_map=ds.celltype_maps[:, :, 0], vf_norm=ds.vf_norm[:, :, 0]
        )

    print("Done")


if __name__ == "__main__":
    main()
