#!/usr/bin/env python3

import pandas as pd
import math

def f_windows_covering_positions(row):
    return len([0 for x in [row['Freq1'], row['Freq2'], row['Freq3']] if not math.isnan(x)])

def f_windows_pass_posterior(row):
    posterior_threshold = 0.8
    return len([0 for x in [row['Post1'], row['Post2'], row['Post3']] if x>=posterior_threshold])

def main(fnames_csv, fname_out):

    # load dataframes
    df = pd.concat([pd.read_csv(fname) for fname in fnames_csv]).reset_index()

    df['coverage'] = df['Ftot'] + df['Rtot']
    df['n_var'] = df['Rvar'] + df['Fvar']
    df['freq'] = df['n_var'] / df['coverage']
    df['sample']=df['file'].str.split('/').str[-2]

    df['windows_pass_posterior'] = df.apply(f_windows_pass_posterior, axis=1).astype(float)
    df['windows_covering_positions'] = df.apply(f_windows_covering_positions, axis=1)

    # filtering
    df = df[(df['windows_covering_positions']>1)  & (df['windows_pass_posterior']>1)]

    # write to csv
    df.to_csv(fname_out)



if __name__ == "__main__":
    main(
        snakemake.input.fnames_csv,
        snakemake.output.fname_out,
    )
