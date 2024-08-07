#!/usr/bin/env python3
"""
Script aggregating diversity measures from all samples into one file.
"""
import pandas as pd

def main(fnames_snv_csv, fout_all_mutations_csv):

    tmp = []

    for f_snv_vcf in fnames_snv_csv:

        f_snv_vcf = str(f_snv_vcf)

        df_tmp = pd.read_csv(f_snv_vcf)
        df_tmp["window_size"] = f_snv_vcf.split("/variants_w")[1].split("/")[0]
        df_tmp["treatment"] = f_snv_vcf.split("/variants")[0].split("/")[-3]
        df_tmp["rep"] = f_snv_vcf.split("/variants")[0].split("/")[-2]
        df_tmp["passage"] = f_snv_vcf.split("/variants")[0].split("/")[-1]

        tmp.append(df_tmp)

    merged_div_csv = pd.concat(
        tmp
    )
    merged_div_csv.to_csv(fout_all_mutations_csv)

if __name__ == "__main__":
    main(
        snakemake.input.fnames_snv_csv,
        snakemake.output.fname_all_mutations,
    )
