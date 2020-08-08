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
    print(path_str)
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

def write_xls(df, coverage_df, config, out):
    
    out_cols = config["fields"].keys()
    nrows, ncols = df[out_cols].shape

    # set up xls file
    writer = pd.ExcelWriter(str(out), engine='xlsxwriter')
    
    # TODO: parametrize these
    df_styled = df[out_cols].style\
    .apply(lambda x: ['background-color: #B0C4DE' if x.primary == False else '' for i in x], axis=1)\
    .apply(lambda x: ['number-format: #,##0']*len(x), subset=['POS'], axis=0)\
    .apply(lambda x: ['number-format: 0; text-align: center']*len(x), subset=['QUAL'], axis=0)\
    .apply(lambda x: ['text-align: center']*len(x), subset=['AD', 'GT', 'AB'], axis=0)\
    .apply(lambda x: x.map({True: 'background-color: red', False: ''}), subset=['clinvar_rev_flag'], axis=0)

    # write main table
    df_styled.to_excel(writer, engine='xlsxwriter', sheet_name='Sheet1', index=False)

    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']
    header_format = workbook.add_format(config["header"]["format"])

    # set col width + name
    for ix, col in enumerate(out_cols):
        specs = config["fields"].get(col, None)
        if specs is None:
            specs = {}

        col_name = specs.get('name', col)
        col_width = specs.get('width', 10)

        # format the column
        worksheet.set_column(ix, ix, col_width)
        # format/rename the header row
        worksheet.write(0, ix, col_name, header_format)


    # set row height of header
    worksheet.set_row(0, 30)

    # add a filter for all columns
    worksheet.autofilter(0, 0, nrows, ncols-1)
    # freeze top row + gene columns
    worksheet.freeze_panes(1, 1)

    coverage_df.to_excel(writer, engine='xlsxwriter', sheet_name='Sheet2', index=False)

    writer.save()



app = typer.Typer()

@app.command()
def main(filename: Path, configpath: Path, coveragepath: Path, outpath: Path):
    # read data
    vcf_df = process_vcf(filename)
    coverage_df = pd.read_csv(str(coveragepath))

    config = yaml.load(open(str(configpath)))
    print(config)
    write_xls(vcf_df, coverage_df, config, outpath)




if __name__ == "__main__":
    app()