#!/usr/bin/env python3
import re
import os
import glob
import numpy as np
import pandas as pd
from cbs import segment, validate
from itertools import cycle
import matplotlib.pyplot as plt
import matplotlib

# cli tool dependencies
import typer
from pathlib import Path

app = typer.Typer()

VALID_CHROMS = ["chr%s" % c for c in range(1,23)] + ["chrX", "chrY"]


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def read_df(path, sample="sample"):
    df = pd.read_csv(path,
       compression="infer", 
       sep="\t", 
       names=["chrom", "start", "end", "name", sample])
    df = df.set_index(["chrom", "start", "end", "name"])
    return df

def write_df(path, df, index=False):
    if path.endswith(".csv"):
        df.to_csv(path, sep=",", index=index)
    elif path.endswith(".tsv"):
        df.to_csv(path, sep="\t", index=index)
    elif path.endswith(".feather"):
        df.reset_index().to_feather(path)
    else:
        raise RuntimeError("Unknown output format, use .csv, .tsv or .feather.")

def load_baseline(inpath):
    if inpath.endswith(".csv"):
        df = pd.read_csv(inpath, sep=",", index=False)
    elif inpath.endswith(".tsv"):
        df = pd.read_csv(inpath, sep="\t", index=False)
    elif inpath.endswith(".feather"):
        df = pd.read_feather(inpath)
    else:
        raise RuntimeError("Unknown baseline format, use .csv, .tsv or .feather.")
    df = df.set_index(["chrom", "start", "end", "name"])
    return df


def get_gene_positions(df):
    # get starts
    d = (df.name != df.name.shift()).astype(int).values
    gene_start_ixs = d.nonzero()[0]
    d = (df.name != df.name.shift(-1)).astype(int).values
    gene_stop_ixs = d.nonzero()[0]
    return list(zip(gene_start_ixs, gene_stop_ixs))

def get_segments(vals):
    L = segment(vals, p=0.0001)
    df = pd.DataFrame()
    df["start"] = validate(vals, L, p=0.0001, shuffles=1000)
    df["end"] = df.start.shift(-1) - 1
    df = df.dropna() # drop last row, which just has last position
    df["mean"] = df.apply(lambda row: np.mean(vals[int(row.start):int(row.end)]), axis=1)
    return df

@app.command()
def baseline(input: Path, output: Path, min_rpkm: float = 5.0):
    if input.is_dir():
        paths = [str(s) for s in input.iterdir()]
        print("directory of %d inputs: %s" % (len(paths), str(input)))
    else:
        paths = glob.glob(str(input))
        print("glob of %d inputs: %s" % (len(paths), str(input)))  

    
    for ix, filepath in enumerate(paths):
        sample_id = os.path.basename(filepath).split(".")[0]
        print(sample_id, filepath)
        if ix == 0:  # first file
            df = read_df(filepath, sample_id)        
        else:  # append columns after
            df[sample_id] = read_df(filepath, sample_id)[sample_id]

    # save file
    write_df(str(output), df)
    
@app.command()
def transform(baseline: Path, sample: Path, output: Path = '/dev/stdout', 
                components_removed: int = 5, min_rpkm: float = 5.0):
    df = load_baseline(str(baseline))
    sample_id = os.path.basename(sample).split(".")[0]
    sample_df = read_df(str(sample), sample_id)

    print("baseline df: ", df.shape)
    print("sample df: ", sample_df.shape)
    #df[sample_id] = sample_df[sample_id]
    #df.loc[sample_df.index, sample_id] = sample_df[sample_id]
    df = pd.merge(df, sample_df, left_index=True, right_index=True)


    rpkm = df / df.sum(axis=0) * 10e9 / 100
    bad_probe_mask = rpkm.median(axis=1) < min_rpkm
    rpkm = rpkm[~bad_probe_mask]
    medians = rpkm.median(axis=1)
    std_devs = rpkm.std(axis=1)
    
    zrpkm = pd.DataFrame(
        data=np.apply_along_axis(lambda val, median, sd: (val - median) / sd, 0, rpkm, medians, std_devs), 
        index=rpkm.index, 
        columns=rpkm.columns
    )

    U, S, Vt = np.linalg.svd(zrpkm.values,full_matrices=False)
    new_S = np.diag(np.hstack([np.zeros([components_removed]),S[components_removed:]]))

    # reconstruct matrix
    new_values = np.dot(U, np.dot(new_S, Vt))
    out = pd.DataFrame(new_values, index=zrpkm.index, columns=zrpkm.columns)
    
    # write single transformed sample to new file
    write_df(str(output), out[sample_id], index=True)

@app.command()
def call(sample: Path, output: Path = '/dev/stdout'):
    df = pd.read_csv(sample, compression="infer", sep=",")
    df = df.set_index(["chrom", "start", "end", "name"])
    sample_id = df.columns[0]
    chromosomes = df.index.get_level_values("chrom").unique()
    chromosomes = natural_sort((set(chromosomes) & set(VALID_CHROMS)))  
        #out_f.write("\t".join(["chrom", "start_ix", "stop_ix", "logr_mean", "num_probes", "genomic_start", "genomic_end"]) + "\n")
    out = []        
    for CHROM in chromosomes:
        chr_df = df.loc[CHROM]
        vals = chr_df[sample_id].values
        for ix, call in get_segments(vals).iterrows():
            start_ix = int(call.start)
            stop_ix = int(call.end)
            region_df = chr_df.iloc[start_ix:stop_ix].reset_index()
            genomic_start, genomic_end = region_df["start"].values[0], region_df["end"].values[-1]
            out.append({
                "chrom": CHROM,
                "start_ix": call["start"],
                "stop_ix": call["end"],
                "logr_mean": call["mean"],
                "num_probes": len(region_df),
                "genomic_start": genomic_start,
                "genomic_end": genomic_end
                })
    out_df = pd.DataFrame(out)
    write_df(str(output), out_df, index=True)



@app.command()
def plot(sample: Path, calls: Path, prefix: str = "", window: int = 20, logr_threshold: float = 2):
    df = pd.read_csv(sample, compression="infer", sep="\t")
    calls_df = pd.read_csv(calls, sep=",")
    sample_id = df.columns[-1]
    y_positions = cycle(np.linspace(-4,-5.75,4))

    for ix, call in calls_df.iterrows():
        if abs(call["logr_mean"]) < logr_threshold:
            continue
        fig, ax = plt.subplots(figsize=(16,6))
        
        start_ix = int(call["start_ix"])
        stop_ix = int(call["stop_ix"])
        window_start_ix = start_ix - window
        window_stop_ix = stop_ix + window
        region_df = df[df["chrom"] == call["chrom"]].iloc[window_start_ix:window_stop_ix]
        vals = region_df.sort_values("start").reset_index()[sample_id].values
        
        zero_stack = np.zeros(len(vals))
        positions = np.repeat(np.arange(0, len(vals)), 3)
        logr = np.vstack([zero_stack, vals.flatten(), zero_stack]).transpose().ravel()

        fig, ax = plt.subplots(figsize=(28,7))        
        ax.plot(positions, logr, color='0.5', marker=None, linewidth=1)
        ax.set_title("%s | %s:%d-%d" % (sample_id, call["chrom"], call["genomic_start"], call["genomic_end"]))
        ax.set_ylim(-6,4)
        # locations
        xtick_locs = np.linspace(0,region_df.shape[0], 10).astype(int)[0:-1]
        xtick_labels = region_df.iloc[xtick_locs]["start"].values
        _ = ax.set_xticks(xtick_locs)
        _ = ax.set_xticklabels(xtick_labels)

        for gene_start_ix, gene_end_ix in get_gene_positions(region_df):
            gene_name = region_df.iloc[gene_start_ix]["name"]
            ypos = next(y_positions)
            ax.add_line(matplotlib.lines.Line2D([gene_start_ix-0.5,gene_end_ix+0.5],[ypos, ypos],color = (102/255.,33/255.,168/255.,0.6), lw=6))
            ax.text(gene_start_ix-0.5,ypos+0.25, gene_name, ha='left',va='center',fontsize=12)
        plt.savefig(f"{prefix}.{str(ix)}.png")
        fig.clf()
        plt.close()

if __name__ == "__main__":
    app()