#!/usr/bin/env python3
import vcf

def f_add_DP_and_AF(record_INFO):
    Rtot = record_INFO['Rtot']
    Ftot = record_INFO['Ftot']
    DP = int(Rtot) + int(Ftot)

    Rvar = record_INFO['Rvar']
    Fvar = record_INFO['Fvar']
    AF = (int(Rvar) + int(Fvar)) / DP

    record_INFO["DP"]= DP
    record_INFO["AF"]= AF

    return record_INFO

def f_windows_covering_positions(record_INFO):
    return len([key for key in record.INFO.keys() if key.startswith("Freq")])

def f_windows_pass_posterior(record_INFO):
    posterior_threshold = 0.8
    post_keys = [key for key in record.INFO.keys() if key.startswith("Post")]
    post_values = [record_INFO[key][0] for key in post_keys]

    return len([0 for x in post_values if float(x)>=posterior_threshold])


def main(fname_snvs_vcf_in, fname_snvs_vcf_out):
    vcf_reader = vcf.Reader(open(fname_snvs_vcf_in, 'r'))
    vcf_writer = vcf.Writer(open(fname_snvs_vcf_out, 'w'), vcf_reader)

    for record in vcf_reader:
        n_windows_covering_positions = f_windows_covering_positions(record.INFO)
        n_windows_pass_posterior = f_windows_pass_posterior(record.INFO)
        if (n_windows_covering_positions>1) & (n_windows_pass_posterior>1):
            record.INFO = f_add_DP_and_AF(record.INFO)
            vcf_writer.write_record(record)


if __name__ == "__main__":
    main(
        snakemake.input.fname_snvs_vcf,
        snakemake.output.fname_snvs_vcf
    )
