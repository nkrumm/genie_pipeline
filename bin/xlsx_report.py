#!/usr/bin/env python3

import yaml
import pandas as pd

# cli tool dependencies
import typer
from pathlib import Path

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
         coverage: Path = None, 
         out: Path = '/dev/stdout'):

    # set up xls file
    writer = pd.ExcelWriter(str(out), engine='xlsxwriter')
    
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