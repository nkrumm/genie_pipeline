#!/usr/bin/env python3

# vcf/xls dependencies
import yaml
import allel
import pandas as pd

# cli tool dependencies
import typer
from pathlib import Path

MAX_ANN_TO_KEEP = 100
ANN_FIELDS = ["ANN_Allele", "ANN_Annotation", "ANN_Annotation_Impact", "ANN_Gene_Name", "ANN_Gene_ID", 
              "ANN_Feature_Type", "ANN_Feature_ID", "ANN_Transcript_BioType", "ANN_Rank", 
              "ANN_HGVS_c", "ANN_HGVS_p", "ANN_cDNA_pos", "ANN_cDNA_length", 
              "ANN_CDS_pos", "ANN_CDS_length", "ANN_AA_pos", "ANN_AA_length", "ANN_Distance"]
ANN_FIELDS_TO_KEEP = ["ANN_Annotation", "ANN_Annotation_Impact", "ANN_Gene_Name", "ANN_Feature_ID", 
                      "ANN_Gene_ID", "ANN_HGVS_c", "ANN_HGVS_p"]


def process_vcf(path):
    # read VCF file; note that the `numbers={'ANN': 100}` means we read in max 100 ANNotations for each variant
    path_str = str(path)
    df = allel.vcf_to_dataframe(path_str, fields=['*'], alt_number=1, numbers={'ANN': MAX_ANN_TO_KEEP},
        fills={'cpdx_AC': 0})
    # read in sample GT/AD data
    callset = allel.read_vcf(path_str, fields=['calldata/GT', 'calldata/AD'], alt_number=1)
    df["GT"] = allel.GenotypeArray(callset['calldata/GT'])[:,0]
    df["GT"] = df.GT.apply(lambda l: "/".join([str(x) for x in l]))

    # calculate allele balance for 0/1 sites and 1/1 sites
    df["AD"] = allel.GenotypeArray(callset['calldata/AD'])[:,0]
    df.loc[(df.GT=="0/1"), "AB"] = df[df.GT=="0/1"].AD.apply(lambda depths: round(depths[0] / (depths[0] + depths[1]), 2))
    df.loc[(df.GT=="1/1"), "AB"] = df[df.GT=="1/1"].AD.apply(lambda depths: round(depths[1] / (depths[0] + depths[1]), 2))
    df["AD"] = df.AD.apply(lambda l: ",".join([str(x) for x in l]))

    # then, coalesce all the ANN_1, ANN_2, ANN_{MAX_ANNMAX_ANN_TO_KEEP} into a single "ANN" column (as a list)
    ANN_COLS = ["ANN_%d" % i for i in range (1,MAX_ANN_TO_KEEP+1)]
    df["ANN"] = df[ANN_COLS].values.tolist()
    # drop broken out ANN_1..ANN_n col
    df = df.drop(ANN_COLS, axis=1)

    # then, coalesce all the FILTER_* fields in the same way
    FILTER_COLS = df.columns[df.columns.str.startswith("FILTER")]
    df["FILTER"] = df[FILTER_COLS].values.tolist()
    # replace "False" with . and report only the unique filter set (as in original VCF)
    df["FILTER"] = df["FILTER"].apply(lambda l: ", ".join(set([str(x) if x is not False else "." for x in l])))
    df = df.drop(FILTER_COLS, axis=1)

    # "explode" the dataframe so each ANN gets its own row
    df = df.explode("ANN")
    # most variants have fewer than MAX_ANN_TO_KEEP annotations, so those have nan and we drop them
    df = df.dropna(subset=["ANN"])

    # reset the index, then mark the first ANN record for each variant as "primary"
    df = df.reset_index(drop=False)
    df.loc[df.groupby(by='index')['POS'].head(1).index, 'primary'] = True
    # secondary transcripts/effects get the False flag here
    df['primary'] = df['primary'].fillna(False)

    # next, for each field, extract the ANN subfields into a 'ANN_parts', which contains a dictionary
    df["ANN_parts"] = df.ANN.apply(lambda ann: dict(zip(ANN_FIELDS, ann.split("|"))))
    for f in ANN_FIELDS_TO_KEEP:
        df[f] = df.ANN_parts.apply(lambda p: p[f])

    # add alignment link
    #https://ncgl.uwcpdx.org/redir/bam/?name=20-10241&path=/samples/20-10241/exome/analyses/02e79501-fb9c-4af9-9eb6-12c11d3fd740/output/alignments/bwa.gatk/20-10241.bwa.gatk.bam&locus=3:98304417-98304517
    df['alignment_link'] = '=HYPERLINK("https://ncgl.uwcpdx.org/redir/bam", "BAM")'

    return df


app = typer.Typer()

@app.command()
def main(in_vcf: Path, out: Path = '/dev/stdout'):
    process_vcf(in_vcf).to_csv(out, sep=",", index=False)


if __name__ == "__main__":
    app()