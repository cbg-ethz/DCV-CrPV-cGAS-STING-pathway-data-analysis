#!/usr/bin/env python3

import pandas as pd


def main(fname_crpv, fname_dcv, fname_out):
    # load dataframes
    df = pd.concat([pd.read_csv(fname_crpv), pd.read_csv(fname_dcv)]).reset_index()

    # great pivot table with
    df_pivot = df.pivot_table(values='coverage', index=['genome', 'sample'], aggfunc='mean')

    df_pivot.to_csv(fname_out)

if __name__ == "__main__":
    main(
        snakemake.input.fname_crpv,
        snakemake.input.fname_dcv,
        snakemake.output.fname_out,
    )
