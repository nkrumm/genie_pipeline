#!/usr/bin/env python3

import pandas as pd

# cli tool dependencies
import typer
from pathlib import Path

THRESHOLDS = ["1X", "10X", "20X", "100X"]

app = typer.Typer()

@app.command()
def main(filename: Path, outpath: Path = '/dev/stdout'):
    df = pd.read_csv(str(filename), sep="\t")
    df["length"] = df["end"] - df["start"]
    # for threshold in THRESHOLDS:
    #     s = df[threshold].sum()/df["length"].sum()*100
    #     print(f"{threshold}:\t{s}")

    for threshold in THRESHOLDS:
        df[f"{threshold}_pct"] = df[threshold]/df["length"].astype(float) * 100

    g = df.groupby(["#chrom", "region"]).agg(
        covered_bp_at_20x = ("20X", sum), 
        n_failing_targets = ("20X_pct", lambda x: sum(x<100)),
        target_length = ("length", sum)
    )
    g["pct_passing_bases"] = g["covered_bp_at_20x"]/g["target_length"]*100
    g["n_failing_bases"] = g["target_length"]-g["covered_bp_at_20x"]
    g = g.reset_index()
    mask = g["pct_passing_bases"] < 100
    cols = ["#chrom", "region", "pct_passing_bases", "n_failing_targets"]
    g[mask][cols].to_csv(str(outpath), index=False)


if __name__ == "__main__":
    app()