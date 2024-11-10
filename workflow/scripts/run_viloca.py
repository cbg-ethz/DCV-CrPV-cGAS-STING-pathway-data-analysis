import subprocess
from pathlib import Path

def main(
    fname_bam,
    fname_reference_dcv,
    fname_reference_crpv,
    fname_reference_p0_dcv,
    fname_reference__p0_crpv,
    fname_snv_vcf,
    fname_snv_csv,
    fname_csv,
    dname_work,
    window_size):

    # get the correct reference
    print(fname_bam)
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
    # n_starts = 1 as default
    # unique modus is on per default

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
            "--threshold",
            "0.0",
            "--min_windows_coverage",
            "1",
            "--NO-strand_bias_filter",
        ],
        cwd=dname_work,
    )


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
        snakemake.params.window_size,
    )
