"""
 input an alignment file at outputs a dataframe where each row is one read with
 columns of the full lenght of the read and the size of the aligned parts of
 the sequence
 """

import pysam
import pandas as pd
from pathlib import Path

# Function to parse the BAM file and create a dataframe
def bam_to_dataframe(bam_file):

    alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    data = []
    bam_file = str(bam_file)
    for read in alignment_file:
        data.append({'read_id': read.query_name,
                     'passage': bam_file.split("/variants")[0].split("/")[-3],
                     'treatment': bam_file.split("/variants")[0].split("/")[-5],
                     'rep': bam_file.split("/variants")[0].split("/")[-4],
                     'full_len': read.query_length,
                     'aligned_len': read.reference_length})

    df = pd.DataFrame(data)
    return df

def main(
    fname_bam,
    fname_csv):

    bam_to_dataframe(fname_bam).to_csv(fname_csv)


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.output.fname_csv),
    )
