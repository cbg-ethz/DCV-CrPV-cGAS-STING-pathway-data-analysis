import subprocess
from pathlib import Path

def main(
    fname_bam,
    fname_reference_dcv,
    fname_reference_crpv,
    fname_reference_p0_dcv,
    fname_reference__p0_crpv,
    window_size,
    fname_snv_vcf,
    fname_snv_csv,
    fname_csv,
    dname_work):

    # get the correct reference
    if "parental" in str(fname_bam):
        if "dcv" in str(fname_bam):
            fname_reference = fname_reference_dcv.resolve()
        elif "crpv" in str(fname_bam):
            fname_reference = fname_reference_crpv.resolve()
    else:
        if "dcv" in str(fname_bam):
            fname_reference = fname_reference_p0_dcv.resolve()
        elif "crpv" in str(fname_bam):
            fname_reference = fname_reference__p0_crpv.resolve()

    alpha = 0.0001
    n_max_haplotypes = 100
    window_size = 50

    dname_work.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "viloca",
            "run",
            "-b",
            fname_bam.resolve(),
            "-f",
            fname_reference,
            "--mode",
            "use_quality_scores",
            "--alpha",
            str(alpha),
            "--n_max_haplotypes",
            str(n_max_haplotypes),
            "-w",
            str(window_size),
            "--min_windows_coverage",
            "2",
            "--NO-strand_bias_filter",
        ],
        cwd=dname_work,
    )

    (dname_work / "snv" / "SNVs_0.010000_final.vcf").rename(fname_snv_vcf)
    (dname_work / "snv" / "SNVs_0.010000_final.csv").rename(fname_snv_csv)
    (dname_work / "snv" / "cooccurring_mutations.csv").rename(fname_csv)


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference_dcv),
        Path(snakemake.input.fname_reference_crpv),
        Path(snakemake.input.fname_reference_p0_dcv),
        Path(snakemake.input.fname_reference__p0_crpv),
        Path(snakemake.output.fname_snv_vcf),
        Path(snakemake.output.fname_snv_csv),
        Path(snakemake.output.fname_csv),
        Path(snakemake.output.dname_work),
    )
