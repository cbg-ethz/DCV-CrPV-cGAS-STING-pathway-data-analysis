#!/usr/bin/env python3
"""
Collect snpgenie results from all samples.
"""
import pandas as pd


def concat_files(in_files_list, out_fname):

    tmp = []

    for file in in_files_list:
        df = pd.read_csv(file, sep='\t')
        df["genotype"] = file.split("/variants")[0].split("/")[-4]
        df["replicate"] = file.split("/variants")[0].split("/")[-2]
        df["passage"] = file.split("/variants")[0].split("/")[-1]

        tmp.append(df)

    pd.concat(tmp).to_csv(out_fname)

def main(
    in_fname_codon,
    in_fname_population_summary,
    in_fname_product_results,
    in_fname_site_results,
    out_fname_codon,
    out_fname_population_summary,
    out_fname_product_results,
    out_fname_site_results,
):

    concat_files(in_fname_codon, out_fname_codon)
    concat_files(in_fname_population_summary, out_fname_population_summary)
    concat_files(in_fname_product_results, out_fname_product_results)
    concat_files(in_fname_site_results, out_fname_site_results)


if __name__ == "__main__":
    main(
        snakemake.input.fname_codon,
        snakemake.input.fname_population_summary,
        snakemake.input.fname_product_results,
        snakemake.input.fname_site_results,
        snakemake.output.fname_codon,
        snakemake.output.fname_population_summary,
        snakemake.output.fname_product_results,
        snakemake.output.fname_site_results,
    )
