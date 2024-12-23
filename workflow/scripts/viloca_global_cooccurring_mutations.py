from Bio import SeqIO
from dataclasses import dataclass
from collections import namedtuple
from re import search
import pandas as pd
import os
"""
Returns csv file.
-- each row is one mutation from one haplotype
-- if a mutation occurrs in several haplotypes there are several rows for this mutation for each haplotype one.
"""


def _deletion_length(seq, char):
    """Determines the length of the deletion. Note that a sequence migth have
    more than one deletion
    seq: substring of the reconstructed haplotype
    char: character that is used to mark a deletion
    """
    count = 0
    for c in seq:
        if c == char:
            count += 1
        else:
            break
    return count

def _count_double_X(ref, seq, x):
    for i in reversed(range(x + 1)):
        if not ref[i] == seq[i] == "X":
            return x - i
    return 0

def __mutations_in_haplotype(ref: str, seq: str, start, av, post, chrom, haplotype_id):
    """
    This function is adapted from _compare_ref_to_read.
    """

    assert len(ref) == len(seq)

    pos = start
    tot_snv = 0
    aux_del = -1

    tmp = [] # collect dicts of mutations on haplotype

    change_in_reference_space = 0

    for idx, v in enumerate(ref):  # iterate on the reference

        if v != seq[idx]:  # SNV detected, save it
            assert not (v != "X" and seq[idx] == "X")

            if seq[idx] == "-" or v == "X": # TODO what is if window starts like that?
                char = "-"
                relevant_seq = seq
                secondary_seq = ref
                if v == "X":
                    char = "X"
                    relevant_seq = ref
                    secondary_seq = seq
                # Avoid counting multiple times a long deletion in the same haplotype
                if idx > aux_del:
                    tot_snv += 1
                    # Check for gap characters and get the deletion
                    # length
                    del_len = _deletion_length(relevant_seq[idx:], char)
                    aux_del = idx + del_len

                    pos_prev = pos - 1
                    num_double_X = _count_double_X(ref, seq, pos_prev - start)
                    secondary_seq = secondary_seq[pos_prev - start - num_double_X] + secondary_seq[
                        (pos - start) : (pos_prev + del_len - start + 1)
                    ] # TODO pos_prev - 1 - beg might be out of range

                    # add deletion to list
                    snv_dict = {
                            "haplotype_id": haplotype_id,
                            "chrom": chrom,
                            "start": start,
                            "position": pos_prev - change_in_reference_space,
                            "ref": ref[pos_prev - start - num_double_X] if v =="X" else secondary_seq,
                            "var": secondary_seq if v =="X" else relevant_seq[pos_prev - start - num_double_X],
                            "reads": av,
                            "support": post,
                        }

                    tmp.append(snv_dict)
            else:
                tot_snv += 1

                snv_dict = {
                            "haplotype_id": haplotype_id,
                            "chrom": chrom,
                            "start": start,
                            "position": pos - change_in_reference_space,
                            "ref": v,
                            "var": seq[idx],
                            "reads": av,
                            "support": post,
                            #"support_normalize": post * av, # TODO Lara: not only post ??
                        }

                tmp.append(snv_dict)

        if v == "X":
            change_in_reference_space += 1

        pos += 1

    return tmp

def get_cooccuring_muts_df(haplo_filename, ref_filename, beg, end, chrom):

    reads = 0.0
    start = int(beg)
    max_snv = -1

    tmp_df = []

    with open(haplo_filename, "rt") as window, open(ref_filename, "rt") as ref:
        d = dict([[s.id, str(s.seq).upper()] for s in SeqIO.parse(ref, "fasta")])
        refSlice = d[chrom]

        for s in SeqIO.parse(window, "fasta"):
            seq = str(s.seq).upper()
            haplotype_id = str(s.id.split("|")[0]) + "-" + beg + "-" + end
            match_obj = search(r"posterior=(.*)\s*ave_reads=(.*)", s.description)
            post, av = float(match_obj.group(1)), float(match_obj.group(2))
            reads += av

            tot_snv = __mutations_in_haplotype(
                refSlice, seq, start, av, post, chrom, haplotype_id
            )
            if tot_snv == []:
                # haplotype is reference
                tot_snv = [{
                            "haplotype_id": "reference"+"-" + beg + "-" + end,
                            "chrom": chrom,
                            "start": start,
                            "reads": av,
                            "support": post,
                        }]

            tmp_df.append(pd.DataFrame(tot_snv))

    if len(tmp_df)>0:
        df = pd.concat(tmp_df)
        df['coverage']= reads
    else:
        df = pd.DataFrame(columns=['coverage'])

    return df

def main(fname, fname_coocc_csv_out):
    base_dir = fname.split("snv/cooccurring_mutations.csv")[0]

    tmp_df =[]
    for root, dirs, files in os.walk(base_dir + "support/"):
        for file in files:
            if file.endswith(".fas"):
                fname_haplo=base_dir+"support/"+file.split("/")[-1]
                fname_ref=base_dir+"raw_reads/"+file.split("/")[-1].split("reads-support")[0]+"ref.fas"
                beg=file.split("/")[-1].split("-")[2]
                end=file.split("/")[-1].split("-")[3].split(".reads")[0]
                chrom="pCrPV" #file.split("/")[-1].split("-")[1]
                tmp_df.append(get_cooccuring_muts_df(fname_haplo, fname_ref, beg,end,chrom))

    pd.concat(tmp_df).to_csv(fname_coocc_csv_out)


if __name__ == "__main__":
    main(
        snakemake.input.fname_csv,
        snakemake.output.fname_csv,
    )
