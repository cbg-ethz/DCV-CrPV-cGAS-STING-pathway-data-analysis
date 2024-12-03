#!/usr/bin/env python3

import pandas as pd
import math

def f_windows_covering_positions(row):
    return len([0 for x in [row['Freq1'], row['Freq2'], row['Freq3']] if not math.isnan(x)])

def f_windows_pass_posterior(row):
    posterior_threshold = 0.8
    return len([0 for x in [row['Post1'], row['Post2'], row['Post3']] if x>=posterior_threshold])

def main(fname_crpv, fname_dcv, fname_out):

    # load dataframes
    df = pd.concat([pd.read_csv(fname_crpv), pd.read_csv(fname_dcv)])

    #INFO field to dataframe
    info_strings = '{"' + df.INFO.str.split(';').str.join('","').str.replace('=','":"').str.replace("\"\",", "") + '"}'
    info_df = pd.json_normalize(info_strings.apply(eval))
    df = pd.concat([df, info_df], axis=1)

    # fiilter columns
    cols_of_interest = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'virus', 'file', 'Fvar', 'Rvar', 'Ftot', 'Rtot',
                        'RefCodon', 'AltCodon', 'RefAminoAcid', 'AltAminoAcid', 'CodonPosition', 'SNPCodonPosition',
                        'AminoAcidChange', 'IsSynonymous',  'Product', 'ProteinID', 'VariantType', 'FeatureType',
                        'Freq1', 'Post1', 'Freq2', 'Post2', 'Freq3', 'Post3']
    df = df[cols_of_interest]

    # convert columns to float
    columns_to_float = ['Freq1', 'Post1','Freq2', 'Post2', 'Freq3', 'Post3', 'Fvar', 'Rvar', 'Ftot', 'Rtot']
    df[columns_to_float] = df[columns_to_float].astype(float)

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
        snakemake.input.fname_crpv,
        snakemake.input.fname_dcv,
        snakemake.output.fname_out,
    )
