#!/usr/bin/env python3
from fuc import pyvcf


def add_vcf_DP_and_AF(vcf_in, vcf_out):
    """
    # alter vcf file such that it matches the necessary input of SNPGenie
    # https://github.com/chasewnelson/SNPGenie#snpgenie-input
    """
    vf = pyvcf.VcfFrame.from_file(vcf_in)

    Rtot = vf.df["INFO"].str.split("Rtot=").str[1].str.split(";").str[0]
    Ftot = vf.df["INFO"].str.split("Ftot=").str[1].str.split(";").str[0]
    DP = Rtot.astype("int") + Ftot.astype("int")

    Rvar = vf.df["INFO"].str.split("Rvar=").str[1].str.split(";").str[0]
    Fvar = vf.df["INFO"].str.split("Fvar=").str[1].str.split(";").str[0]
    AF = (Rvar.astype("int") + Fvar.astype("int")) / DP

    vf.df["INFO"] = vf.df["INFO"] + ";DP=" + DP.astype("str")
    vf.df["INFO"] = vf.df["INFO"] + ";AF=" + AF.astype("str")

    vf.to_file(vcf_out)


def main(fname_snvs_vcf_in, fname_snvs_vcf_out):

    # add DP4 to vcf file
    add_vcf_DP_and_AF(fname_snvs_vcf_in, fname_snvs_vcf_out)


if __name__ == "__main__":
    main(
        snakemake.input.fname_snvs_vcf,
        snakemake.output.fname_snvs_vcf
    )
