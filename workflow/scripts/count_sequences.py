import os

def main(data_dir, output_file):
    # Open output file for writing
    with open(output_file, 'w') as outfile:
        outfile.write("File\tNumber of Sequences\n")

        # Walk through the directory and subdirectories
        for root, dirs, files in os.walk(data_dir):
            for file in files:
                if file.endswith(".fastq"):
                    # Run seqtk to count sequences and capture output
                    count = os.popen("seqtk seq -A " + os.path.join(root, file) + " | wc -l").read().strip()
                    # Write filename and sequence count to output file
                    outfile.write(f"{file}\t{count}\n")


if __name__ == "__main__":
    main(
        snakemake.input.datadir,
        snakemake.output.fname,
    )
