#!/usr/bin/env python

# Written by Maxime Borry and released under the MIT license. 
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
import os
import pandas as pd
import numpy as np

def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='merge_kraken_res',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Merging csv count files in one table')
    parser.add_argument(
        '-or',
        dest="readout",
        default="kraken_read_count_table.csv",
        help="Read count output file. Default = kraken_read_count_table.csv")
    parser.add_argument(
        '-ok',
        dest="kmerout",
        default="kraken_kmer_unicity_table.csv",
        help="Kmer unicity output file. Default = kraken_kmer_unicity_table.csv")
    parser.add_argument(
        '-okc',
        dest="kmercountsout",
        default="kraken_kmer_counts_table.csv",
        help="unique kmer counts output file. Default = kraken_kmer_counts_table.csv")    
    parser.add_argument(
        '--bracken', 
        dest='bracken', 
        default=False, 
        action='store_true')

    args = parser.parse_args()

    readout = args.readout
    kmerout = args.kmerout
    kmercountsout = args.kmercountsout
    bracken = args.bracken

    return(readout, kmerout, kmercountsout, bracken)

def get_csv(bracken):
    if bracken:
        tmp = [i for i in os.listdir() if ".csv" in i]
        read = [i for i in tmp if '.read_' in i]
        return(read)
    else:
        tmp = [i for i in os.listdir() if ".csv" in i]
        kmer = [i for i in tmp if '.kmer_kraken' in i]
        kmercounts = [i for i in tmp if '.kmer_counts_kraken' in i]
        read = [i for i in tmp if '.read_kraken' in i]
        return(read, kmer, kmercounts)

def _get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)

def merge_csv(all_csv):
    df = pd.read_csv(all_csv[0], index_col=0)
    for i in range(1, len(all_csv)):
        df_tmp = pd.read_csv(all_csv[i], index_col=0)
        df = pd.merge(left=df, right=df_tmp, on='TAXID', how='outer')
    df.fillna(0, inplace=True)
    return(df)

def write_csv(pd_dataframe, outfile):
    pd_dataframe.to_csv(outfile)

if __name__ == "__main__":
    READOUT, KMEROUT, KMERCOUNTSOUT, BRACKEN = _get_args()
    if BRACKEN:
        reads = get_csv(bracken=True)
        read_df = merge_csv(reads)
        write_csv(read_df, READOUT)
    else:
        reads, kmers, kmercounts = get_csv(bracken=False)
        read_df = merge_csv(reads)
        kmer_df = merge_csv(kmers)
        kmercounts_df = merge_csv(kmercounts)
        write_csv(read_df, READOUT)
        write_csv(kmer_df, KMEROUT)
        write_csv(kmercounts_df, KMERCOUNTSOUT)