import os
import re
from collections.abc import Callable
from datetime import datetime, timedelta
from io import StringIO

import anndata as ad
import pandas as pd


def celltype_signatures(
    adata: ad.AnnData,
    *,
    celltype_col: str = "leiden",
    agg_method: str | Callable = "mean",
) -> pd.DataFrame:
    """
    Calculate gene expression signatures per cluster.

    Parameters
    ----------
    adata : anndata.AnnData
        :py:class:`anndata.AnnData` object of kernel density estimates and cluster /
        celltype assignments.
    cluster : str, optional
        Name of column in :py:attr:`anndata.AnnData.obs` containing clustering
        information.
    agg_method : str or collections.abc.Callable, optional
        Function to aggregate gene expression per cluster used by
        :py:meth:`pandas.DataFrame.agg`.

    Returns
    -------
    pandas.DataFrame
        :py:class:`pandas.DataFrame` of gene expression aggregated per `cluster`.
    """
    signatures = (
        adata.to_df()
        .merge(adata.obs[celltype_col], left_index=True, right_index=True)
        .groupby(celltype_col, observed=True, sort=False)
        .agg(agg_method)
        .transpose()
        .rename_axis(adata.var_names.name)
    )

    signatures /= signatures.sum(axis=0)

    return signatures


def parse_sacct_time(time_str):
    pattern = r"(?:(?P<day>\d+)-)?(?:(?P<hour>\d{1,2}):)?(?P<min>\d{2}):(?P<sec>\d{2})(?:.(?P<usec>\d+))?"
    match = re.match(pattern, time_str)

    if not match:
        raise ValueError("Invalid SLURM time format")

    days = int(match.group("day")) if match.group("day") else 0
    hours = int(match.group("hour")) if match.group("hour") else 0
    minutes = int(match.group("min"))
    seconds = int(match.group("sec"))
    microseconds = int(match.group("usec")) if match.group("usec") else 0

    time = timedelta(
        days=days,
        hours=hours,
        minutes=minutes,
        seconds=seconds,
        microseconds=microseconds,
    )

    return time.total_seconds()


def get_runtime_stats(ids: dict):
    job_ids = pd.Series(ids, name="tool").to_frame()

    job_stats = os.popen(
        (
            "sacct "
            f"-j {','.join(job_ids.index.astype(str))} "
            "--format='JobID,Jobname%50,TotalCPU,ElapsedRaw,MaxRSS' "
            "-P --delimiter=$'\t' "
            "--units=M "
        )
    ).read()

    job_stats = pd.read_table(StringIO(job_stats))

    cpu_stats = (
        job_stats.loc[
            lambda df: ~df["JobID"].str.contains(".", regex=False),
            ["JobID", "JobName", "TotalCPU", "ElapsedRaw"],
        ]
        .assign(
            JobID=lambda df: df["JobID"].astype(int),
            TotalCPU=lambda df: df["TotalCPU"].map(parse_sacct_time),
        )
        .set_index("JobID")
        .rename(columns={"TotalCPU": "CPU time [s]", "ElapsedRaw": "wall time [s]"})
    )

    memory_stats = (
        job_stats.loc[
            lambda df: df["JobID"].str.contains(".batch", regex=False),
            ["JobID", "MaxRSS"],
        ]
        .assign(
            JobID=lambda df: df["JobID"].str.extract("(\d+)").astype(int),
            MaxRSS=lambda df: df["MaxRSS"].str.extract("([\d\\.]+)").astype(float),
        )
        .set_index("JobID")
        .rename(columns={"MaxRSS": "max memory [MB]"})
    )

    return (
        cpu_stats.join(memory_stats)
        .join(job_ids, how="inner")
        .drop(columns=["JobName"])
    )
