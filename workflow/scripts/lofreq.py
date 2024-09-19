import subprocess
from pathlib import Path
import pandas as pd

def main(
    fname_bam,
    fname_reference_dcv,
    fname_reference_crpv,
    fname_reference_p0_dcv,
    fname_reference__p0_crpv,
    fname_results_snv,
    dname_work):

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


    dname_work.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "lofreq",
            "call",
            "-f",
            fname_reference.resolve(),
            "-o",
            fname_results_snv.resolve(),
            fname_bam.resolve(),
        ],
        cwd=dname_work,
        check=True,
    )


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference_dcv),
        Path(snakemake.input.fname_reference_crpv),
        Path(snakemake.input.fname_reference_p0_dcv),
        Path(snakemake.input.fname_reference__p0_crpv),
        Path(snakemake.output.fname_snv_vcf),
        Path(snakemake.output.dname_work),
    )
