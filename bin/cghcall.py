import numpy as np
import pandas as pd
import rpy2
import rpy2.rlike.container as rlc
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

cghcall = importr("CGHcall")
cghbase = importr("CGHbase")

rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda f: None

CHR2INT = {"chr%d" % d: d for d in range(1,23)}
CHR2INT["chrX"] = 23
CHR2INT["chrY"] = 24

callint2colname = {-2: "p_dloss", -1: "p_loss", 1: "p_gain", 2:"p_amp"}

VALID_CHROMS = ["chr%s" % c for c in range(1,23)] + ["chrX", "chrY"]

def make_calls(df, sample_id):
    in_df = df.reset_index()
    in_df = in_df[in_df.chrom.isin(VALID_CHROMS)]
    in_df["chrom"] = in_df.chrom.apply(lambda c: CHR2INT[c])
    od = rlc.OrdDict([
        ("NAME", robjects.vectors.IntVector(in_df.index.values)),
        ("CHROMOSOME",  robjects.vectors.IntVector(in_df.chrom.values)),
        ("START_POS",   robjects.vectors.IntVector(in_df.start.values)),
        ("STOP_POS",    robjects.vectors.IntVector(in_df.end.values)),
        (str(sample_id), robjects.vectors.FloatVector(in_df[sample_id].values))
    ])

    rdf = robjects.DataFrame(od)

    t  = cghbase.make_cghRaw(rdf)
    t2 = cghcall.preprocess(t)
    t3 = cghcall.segmentData(t2,**{'alpha':0.02, 'undo.splits': "sdundo", 'undo.SD':2})
    t4 = cghcall.postsegnormalize(t3)
    t5 = cghcall.CGHcall(t4)
    try:
        t6 = cghcall.ExpandCGHcall(t5,t4)
    except:
        print("Exception")
        return (t4, t5)
    out=np.vstack([np.array(cghbase.chromosomes(t6)),
                   np.array(cghbase.bpstart(t6)),
                   np.array(cghbase.bpend(t6)),
                   np.array(cghbase.copynumber(t6)).flatten(),
                   np.array(cghbase.calls(t6)).flatten(),
                   np.array(cghbase.probdloss(t6), dtype=np.float).flatten(),
                   np.array(cghbase.probloss(t6), dtype=np.float).flatten(),
                   np.array(cghbase.probnorm(t6), dtype=np.float).flatten(),
                   np.array(cghbase.probgain(t6), dtype=np.float).flatten(),
                   np.array(cghbase.probamp(t6), dtype=np.float).flatten()]).transpose()

    out_df = pd.DataFrame(out,columns=["chromosome","start","stop","svdzrpkm","call","p_dloss","p_loss","p_norm","p_gain","p_amp"])
    # post process 
    out_df["call"] = np.sign(out_df["call"])
    out_df["p_gain"] = np.maximum(out_df["p_gain"], out_df["p_amp"])
    out_df["p_loss"] = np.maximum(out_df["p_loss"], out_df["p_dloss"])

    return out_df


def findBreakPoints(x):
    #bkpts = np.where(np.logical_xor(x[1:],x[:-1],))[0]
    bkpts = np.where(np.logical_xor(np.hstack([x[0],x]),np.hstack([x,x[-1]])))[0]
    if x[0] == True:
        bkpts = np.hstack([0,bkpts])
    if x[-1] == True:
        bkpts = np.hstack([bkpts,len(x)])
    if len(bkpts) <= 2:
        if bkpts[1]-bkpts[0] == 1:
            return [tuple([bkpts[0],bkpts[0]])]
        return [tuple(bkpts)]
    else:
        return [(b[0],b[0]) if (b[1]-b[0])==1 else tuple(b) for b in bkpts.reshape(int(bkpts.shape[0]/2),2)]


def segment(cghcall_output, sample_id):
    c = []
    for chrom, out in cghcall_output.groupby("chromosome"):
        call_states_present = set(out["call"])-set([0])
        for i in call_states_present:
            calls = findBreakPoints(np.array(1 * (out["call"] == i)))
            for call in calls:
                prob_vals = out.iloc[call[0]:call[1]-1][callint2colname[i]]
                svdzrpkm_vals = out.iloc[call[0]:call[1]-1]["svdzrpkm"]
                start_bp = out.iloc[call[0]]["start"]
                stop_bp = out.iloc[call[1]-1]["stop"]
                chromosome_vals = np.unique(out.iloc[call[0]:call[1]-1]["chromosome"])
                if len(chromosome_vals) == 1:
                    chromosome = chromosome_vals[0]
                else:
                    print("[ERROR] multiple chromosomes in this call!, Ignoring")
                    print(chromosome_vals, start_bp, stop_bp, call)
                    continue
                start_exon= np.where(out[out["chromosome"]==chromosome]["start"]>=start_bp)[0][0]
                stop_exon = np.where(out[out["chromosome"]==chromosome]["start"]<=stop_bp)[0][-1]
                #np.median(prob_vals), np.median(svdzrpkm_vals), np.std(svdzrpkm_vals))
                t = {"sample_id": sample_id,\
                     "chromosome":int(chromosome),\
                     "start":int(start_bp),\
                     "stop":int(stop_bp),\
                     "state":int(i),\
                     "start_exon":int(start_exon),\
                     "stop_exon":int(stop_exon),\
                     "num_probes":int(stop_exon-start_exon+1),\
                     "size_bp":int(stop_bp - start_bp),\
                     "probability": np.median(prob_vals),
                     "median_svdzrpkm": np.median(svdzrpkm_vals),
                     "stdev_svdzrpkm":np.std(svdzrpkm_vals)}
                c.append(t)
    return pd.DataFrame(c)