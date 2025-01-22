#! /usr/bin/env python
import pysam
import pandas as pd
import numpy as np
import pyranges as pr
from cigar import Cigar

# position of SNV 5764
start_region = 5764 -1
end_region = 5774
print('start ', start_region, ' end ', end_region)
deletion_postion = list(range(5764-1,5774))

def f_filter_region(row):
    if (float(row['Start'])<5764-1) & (float(row['End'])>5774):
        return True
    else:
        return False

def postiion_of_deletion(cigar, start_position):
    c = Cigar(cigar)
    del_list = []  # list of tuples (start_position, del_length)
    # iterate through cigar
    for c_item in list(c.items()):
        if c_item[1]=="M":
            start_position+=c_item[0]
        elif c_item[1]=="D":
            del_list = del_list + [(start_position, c_item[0])]
            start_position+=c_item[0]
        elif c_item[1]=="I":
            start_position+=0
        elif c_item[1]=="S":
            start_position+=0
        else:
            print(c_item[1])
        # ignore I = "insertion", ignore "S"
    return del_list

def f_del_in_our_region(row):
    del_list = postiion_of_deletion(row["Cigar"], row["Start"])
    if len(del_list)>0:
        for deletion in del_list:
            if deletion[0] in list(range(start_region-2, end_region+2)):
                return True
    return False

def f_qual_del_in_region(row):
    len_region = end_region+2 - (start_region-2)
    index_start_region = start_region-2 - row['Start']
    if index_start_region<0:
            return np.mean(row['Quality'][row['QueryStart']+0:row['QueryStart']+0+len_region+1])
    else:
        return np.mean(row['Quality'][row['QueryStart']+index_start_region:row['QueryStart']+index_start_region+len_region+1])


def f_del_region_close_to_read_end(row):
    if row['Start'] in range(start_region-5, start_region+5):
        return True
    if row['End'] in range(end_region-5, end_region+5):
        return True
    else:
        return False


def get_df_bam(fname_bam):

    fname_bam = str(fname_bam)
    df_bam = pr.read_bam(fname_bam, sparse=False, as_df=True)
    df_bam['covers_region'] = df_bam.apply(f_filter_region, axis=1)
    df_bam = df_bam[df_bam['covers_region']==True]

    df_bam['del_in_region'] = df_bam.apply(f_del_in_our_region, axis=1)
    df_bam['qual_del_in_region'] = df_bam.apply(f_qual_del_in_region, axis=1)
    df_bam['del_region_close_to_read_end'] = df_bam.apply(f_del_region_close_to_read_end, axis=1)

    df_bam = df_bam[['Start', 'End', 'Strand', 'Flag', 'QueryStart', 'QueryEnd', 'Name', 'Cigar', 'covers_region', 'del_in_region', 'qual_del_in_region', 'del_region_close_to_read_end']]

    # add sample information
    df_bam["file"] = fname_bam
    df_bam["snv_position"] = "5764"

    return df_bam



def main(fname_bams, fname_out):

    #df = pd.DataFrame()

    #for fname_bam in fname_bams:
    #    df = pd.concat([df, get_df_bam(fname_bam)], ignore_index=True)

    get_df_bam(fname_bams).to_csv(fname_out)

if __name__ == "__main__":
    main(
        snakemake.input.fname_bams,
        snakemake.output.fname_out,
    )
