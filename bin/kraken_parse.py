#!/usr/bin/env python

# Written by Maxime Borry and released under the MIT license. 
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
import csv

def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='kraken_parse',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Parsing kraken')
    parser.add_argument('krakenReport', help="path to kraken report file")
    parser.add_argument(
        '-c',
        dest="count",
        default=50,
        help="Minimum number of hits on clade to report it. Default = 50")
    parser.add_argument(
        '-or',
        dest="readout",
        default=None,
        help="Read count output file. Default = <basename>.read_kraken_parsed.csv")
    parser.add_argument(
        '-ok',
        dest="kmerout",
        default=None,
        help="Kmer Output file. Default = <basename>.kmer_kraken_parsed.csv")
    parser.add_argument(
        '-okc',
        dest="kmercountsout",
        default=None,
        help="Unique kmer counts output file. Default = <basename>.kmer_counts_kraken_parsed.csv")
    parser.add_argument('--bracken', dest='bracken', default=False, action='store_true')

    args = parser.parse_args()

    infile = args.krakenReport
    countlim = int(args.count)
    readout = args.readout
    kmerout = args.kmerout
    kmercountsout = args.kmercountsout
    bracken = args.bracken

    return(infile, countlim, readout, kmerout, kmercountsout, bracken)


def _get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)


def parse_kraken(infile, countlim, bracken):
    '''
    INPUT:
        infile (str): path to kraken report file
        countlim (int): lowest count threshold to report hit
    OUTPUT:
        resdict (dict): key=taxid, value=readCount
    BRACKEN:
        bracken (logical): is infile the output from bracken?

    '''
    with open(infile, 'r') as f:
        read_dict = {}
        kmer_dict = {}
        kmer_counts = {}
        csvreader = csv.reader(f, delimiter='\t')
        for line in csvreader:
            reads = int(line[1])
            if reads >= countlim:
                if bracken:
                    taxid = line[4].strip().replace(",","")
                    read_dict[taxid] = reads
                else:
                    taxid = line[6].strip().replace(",","")
                    kmer = line[3]
                    unique_kmer = line[4]
                    try:
                        kmer_duplicity = float(kmer)/float(unique_kmer)
                    except ZeroDivisionError:
                        kmer_duplicity = 0
                    kmer_dict[taxid] = kmer_duplicity
                    read_dict[taxid] = reads
                    kmer_counts[taxid] = unique_kmer
        if bracken:
            return(read_dict)
        else:
            return(read_dict, kmer_dict, kmer_counts)


def write_output(resdict, infile, outfile):
    with open(outfile, 'w') as f:
        basename = _get_basename(infile)
        f.write(f"TAXID,{basename}\n")
        for akey in resdict.keys():
            f.write(f"{akey},{resdict[akey]}\n")


if __name__ == '__main__':
    INFILE, COUNTLIM, readout, kmerout, kmercountsout, bracken = _get_args()

    if not readout:
        read_outfile = _get_basename(INFILE)+".read_kraken_parsed.csv"
    else:
        read_outfile = readout
    if not kmerout:    
        kmer_outfile = _get_basename(INFILE)+".kmer_kraken_parsed.csv"
    else:
        kmer_outfile = kmerout
    if not kmercountsout:
        kmer_countsfile = _get_basename(INFILE)+".kmer_counts_kraken_parsed.csv"
    else:
        kmer_countsfile = kmercountsout

    if bracken:
        read_dict = parse_kraken(infile=INFILE, countlim=COUNTLIM, bracken=True)
        write_output(resdict=read_dict, infile=INFILE, outfile=read_outfile)
    else:
        read_dict, kmer_dict, kmer_counts = parse_kraken(infile=INFILE, countlim=COUNTLIM, bracken=False)
        write_output(resdict=read_dict, infile=INFILE, outfile=read_outfile)
        write_output(resdict=kmer_dict, infile=INFILE, outfile=kmer_outfile)
        write_output(resdict=kmer_counts, infile=INFILE, outfile=kmer_countsfile)

