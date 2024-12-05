#!/usr/bin/env python3
import subprocess


def run_snpgenie(fname_reference, in_vcf, gtffile_CDS_annotations, dname_work):
    """
    WITHIN-POOL ANALYSIS. Use snpgenie.pl, the original SNPGenie.
    Analyzes within-sample πN/πS from pooled NGS SNP data.

    fname_reference: fasta of reference sequence wrt which the SNVs where called
    """
    subprocess.run(
        [
            "snpgenie.pl",
            "--vcfformat=2",
            "--snpreport=" + str(in_vcf),
            "--fastafile=" + str(fname_reference),
            "--gtffile=" + str(gtffile_CDS_annotations),
            "--minfreq=0.0001",
            "--outdir="+str(dname_work)
        ],
        check=True,
        #cwd=dname_work,
    )


def main(fname_reference, fname_snv_in, gtffile_CDS_annotations, dname_work):

    from pathlib import Path
    import shutil

    test = Path(dname_work)
    if test.exists() and test.is_dir():
        shutil.rmtree(test)

    # run snpgenie
    run_snpgenie(fname_reference, fname_snv_in, gtffile_CDS_annotations, dname_work)


if __name__ == "__main__":
    main(
        snakemake.input.fname_reference,
        snakemake.input.fname_snvs_vcf,
        snakemake.input.fname_gtffile,
        snakemake.output.dname_work
    )
