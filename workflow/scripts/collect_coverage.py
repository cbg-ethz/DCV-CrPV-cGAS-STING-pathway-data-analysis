import pandas as pd
import yaml

def collect_coverage(fnames, fname_out):
    tmp =[]

    for fname in fnames:
        df_tmp =  pd.read_csv(fname, sep="\t")
        df_tmp['genome'] = fname.split("/")[-4]
        df_tmp['sample'] = fname.split("/")[-3]
        sample_id = fname.split("/")[-4]+"/"+fname.split("/")[-3]
        df_tmp['coverage'] = df_tmp[sample_id]
        tmp.append(df_tmp[['genome', 'sample', 'coverage', 'ref', 'pos']])

    pd.concat(tmp).to_csv(fname_out)

def collect_read_len(fnames, fname_out):

    genome = []
    sample = []
    mean_read_len = []

    for fname in fnames:
        with open(fname, 'r') as file:
            info = yaml.safe_load(file)
            id = fname.split("/")[-4] + "/" + fname.split("/")[-3]
            genome.append(fname.split("/")[-4])
            sample.append(fname.split("/")[-3])
            mean_read_len.append(info[id]['readlen_mean'])

    pd.DataFrame({
        "genome": genome,
        "sample": sample,
        "mean_read_len": mean_read_len
    }).to_csv(fname_out)


def main(fnames, fname_out, fname_out_read_len):

    collect_coverage(fnames, fname_out)

    fnames_yaml = [fname.split("coverage")[0]+"REF_aln_stats.yaml" for fname in fnames]
    collect_read_len(fnames_yaml, fname_out_read_len)

if __name__ == "__main__":
    main(
        snakemake.input.fnames_coverage,
        snakemake.output.fname_coverage,
        snakemake.output.fname_read_len,
    )
