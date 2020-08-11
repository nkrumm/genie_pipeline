#!/usr/bin/env python3

import yaml
import json
import datetime
import pandas as pd

# cli tool dependencies
import typer
from pathlib import Path

def info_sheet(data, sheet_name, writer):
    workbook  = writer.book
    worksheet =  workbook.add_worksheet(sheet_name)
    header = workbook.add_format({'bold': True, 'font_size': 20})
    row_ix = 0
    worksheet.merge_range(row_ix, 0, row_ix, 1, "UW Clinical Exome Variant Report: %s" % data["samples"], header)
    row_ix += 1
    worksheet.write(row_ix, 0, "Started at")
    worksheet.write(row_ix, 1, data["started_at"])
    row_ix += 1
    worksheet.write(row_ix, 0, "Generated at")
    worksheet.write(row_ix, 1, datetime.datetime.now().strftime("%b %d, %Y at %H:%M"))
    row_ix += 1
    worksheet.write(row_ix, 0, "Pipeline version")
    worksheet.write(row_ix, 1, data["pipeline_version"])
    row_ix += 1
    worksheet.write(row_ix, 0, "Nextflow hash")
    worksheet.write(row_ix, 1, f"{data['nextflow_script_id']}")
    row_ix += 1
    worksheet.write(row_ix, 0, "Nextflow version")
    worksheet.write(row_ix, 1, f"{data['nextflow_version']['major']}.{data['nextflow_version']['minor']}.{data['nextflow_version']['patch']}")
    row_ix += 2
    worksheet.write(row_ix, 0, "Input VCF")
    worksheet.write(row_ix, 1, data["input_vcf"])
    row_ix += 2
    worksheet.write(row_ix, 0, "Samples Included")
    worksheet.write(row_ix, 1, data["samples"])
    
    worksheet.set_column(0, 0, 20)
    worksheet.set_column(1, 1, 100)

    return writer


def snv_sheet(df, sheet_name, writer, config):
    out_cols = config["fields"].keys()
    nrows, ncols = df[out_cols].shape

    
    # TODO: parametrize these
    df_styled = df[out_cols].style\
    .apply(lambda x: ['background-color: #B0C4DE' if x.primary == False else '' for i in x], axis=1)\
    .apply(lambda x: ['number-format: #,##0']*len(x), subset=['POS'], axis=0)\
    .apply(lambda x: ['number-format: 0; text-align: center']*len(x), subset=['QUAL'], axis=0)\
    .apply(lambda x: ['text-align: center']*len(x), subset=['AD', 'GT', 'AB'], axis=0)\
    .apply(lambda x: x.map({True: 'background-color: red', False: ''}), subset=['clinvar_rev_flag'], axis=0)

    # write main table
    df_styled.to_excel(writer, engine='xlsxwriter', sheet_name=sheet_name, index=False)

    workbook  = writer.book
    worksheet = writer.sheets[sheet_name]
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
    return writer


def coverage_sheet(df, sheet_name, writer):
    df.to_excel(writer, engine='xlsxwriter', sheet_name=sheet_name, index=False)
    return writer

app = typer.Typer()

@app.command()
def main(variants: Path = typer.Option(...), 
         config: Path = typer.Option(...), 
         info: Path = None,
         coverage: Path = None, 
         out: Path = '/dev/stdout'):

    # set up xls file
    writer = pd.ExcelWriter(str(out), engine='xlsxwriter')

    if info:
        writer = info_sheet(
            data = json.load(open(info)),
            sheet_name = "Info",
            writer = writer
        )
    
    writer = snv_sheet(
        df = pd.read_csv(str(variants)),
        sheet_name = "SNVs+Indels",
        writer = writer,
        config = yaml.load(open(str(config)))
    )
    if coverage:
        writer = coverage_sheet(
            df = pd.read_csv(str(coverage)),
            sheet_name="Coverage",
            writer = writer
        )

    writer.save()



if __name__ == "__main__":
    app()