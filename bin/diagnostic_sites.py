###This script will identify positions in an alignment that are fixed for differential alleles between two populations
###Takes a fasta alignment as input, and a file assigning sequences to one of two populations
###The popfile should have the sample/header name in the first column and the population assignment in the second.
###Samples that doesn't have a population assignment will be ignored in this analysis


#Import the modules to be used
import os
import sys
import pysam #this is a module for parsing alignment files, will need to be installed
from Bio import SeqIO
import argparse

############ Defining the functions we will use


##### these would be reference mt genome alignment

#This is the function that open and transforms the fasta file to a dictionary in its simplest use (it also subsets by position and sample if wanted, that's why its long and messy.).
def FastaFetch(infasta, sequences=None, start=None, stop=None):
    f = pysam.infasta(infasta)
    seqDict = []
    if sequences is not None:
        if os.path.isfile(sequences):
            seqList = []
            with open(sequences) as s:
                for line in s.readlines():
                    if not line == "":
                        seqList.append(line.rstrip().lstrip(">")) 
        else:
            seqList = []
            for s in sequences.split(","):
                if not s == "":
                    seqList.append(s.rstrip().lstrip(">"))
        if start and stop:
            for seq in seqList:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start, stop)})
        elif start:
            for seq in seqList:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start)})
        else:
            for seq in seqList:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq)})
    else:
        if start and stop:
            seqDict = []
            for seq in f.references:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start, stop)})
        elif start:
            for seq in f.references:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq, start)})
        else:
            for seq in f.references:
                seqDict.append({'header':seq, 'sequence':f.fetch(seq)})
    return seqDict

#This function checks pop1's genotypes and when that's fixed, compares it to pop2's genotype and outputs the position and bases if different
def FindDiagnosticPositions(seqDict, popfile, reference):
    with open(popfile) as pf:
        popmap = pf.read().split("\n")
        popdict = []
        for sample in popmap:
            if len(sample.split("\t")) == 2:
                popdict.append({'header':sample.split("\t")[0], 'pop':sample.split("\t")[1]})
    
    #a dictionary with pop as key and samples assigned as values
    #first get the names of the populations
    pop1 = popdict[0]['pop']
    for pop in popdict:
        if not pop['pop'] == pop1:
            pop2 = pop['pop']
            break
    #then put in the samples:
    popmap = {pop1: [], pop2: []}
    for i in popdict:
        if i['pop'] == pop1:
            popmap[pop1].append(i['header'].lstrip(">"))
        elif i['pop'] == pop2:
            popmap[pop2].append(i['header'].lstrip(">"))
    print("assigning the following samples as references: ")
    print(popmap)
    #go through all the bases in the alignment and check them
    diagnosticSites = []
    for base in range(0, len(seqDict[0]['sequence'])):
        pop1_alleles = []
        pop2_alleles = []
        for seq in seqDict:
            if seq['header'] in popmap[pop1]:
                if not seq['sequence'][base] in ["N","n","-","?"]:
                    pop1_alleles.append(seq['sequence'][base])
            elif seq['header'] in popmap[pop2]:
                if not seq['sequence'][base] in ["N","n","-","?"]:
                    pop2_alleles.append(seq['sequence'][base])
        if len(set(pop1_alleles)) == 1:
            if len(set(pop2_alleles)) == 1:
                if not pop1_alleles[0] == pop2_alleles[0]:
                    diagnosticSites.append((base, {pop1:pop1_alleles[0], pop2:pop2_alleles[0]}))
    if len(diagnosticSites) == 0:
        print("No fixed differences were found")
    #then we write those to a tsv reference:
    with open(reference, "w") as of:
        of.write("Position\t" + pop1 + "\t" + pop2 + "\n")
        for site in diagnosticSites:
            of.write(str(site[0] + 1) + "\t" + site[1][pop1] + "\t" + site[1][pop2] + "\n")
        of.close()

#This function checks pop1's genotypes and when that's fixed, compares it to pop2's genotype and outputs the position and bases if different
def CountDiagnosticPositions(seqDict, reference, outputcounts):
    #First line should hold the popnames
    with open(reference) as f:
        header = f.readline().split("\t")
        pop1 = header[1]
        pop2 = header[2].rstrip()
        #and line 1- should have the sites
        f.seek(0)
        sites = f.read().split("\n")[1:]
        distanceDict = []
        for seq in seqDict:
            #count number of similar sites per sample, start at 0 for each sample
            pop1_similarity = 0
            pop2_similarity = 0
            other = 0
            missing_data = 0
            for site in sites:
                if site.split("\t")[0] == "": #if position is a blank line, break and start with next sample
                    break
                pos = int(site.split("\t")[0]) - 1 #zero based position in python, so if input is not one-based this needs changing (it is one-based when using FindDiagnosticSites.py output)
                pop1_allele = site.split("\t")[1]
                pop2_allele = site.split("\t")[2]
                #count the similarities: for now it's still case sensitive! 
                if seq['sequence'][pos] == pop1_allele:
                    pop1_similarity += 1
                elif seq['sequence'][pos] == pop2_allele:
                    pop2_similarity += 1
                elif seq['sequence'][pos] in ['N', 'n', '-', '?']:
                    missing_data += 1
                else:
                    other += 1
            #add the samples stats to dictionary
            distanceDict.append({'sample':seq['header'], pop1:pop1_similarity,
            pop2:pop2_similarity, 'other':other, 'missing_data':missing_data})
    #then we write those to a tsv outfile:
    with open(outputcounts, "w") as of:
        of.write("Sample\t" + pop1 + "\t" + pop2 + "\t" + "other" + "\t" + "missing" + "\n")
        for sample in distanceDict:
            of.write(sample['sample'] + "\t" + str(sample[pop1]) + "\t" + str(sample[pop2]) + "\t" + str(sample['other']) + "\t" + str(sample['missing_data']) + "\n")
        of.close()

# TODO: check that this writes properly
def ExtractDiagSites(diagsites,allelecounts,extractedsites):
    with open(diagsites, "r") as diag_sites:
            lines=[line.rstrip() for line in diag_sites]

    allele_count=open(allelecounts,"r")
    for line in allele_count:
            values=line.split(" ")
            if values[0] in lines:
                    with open(extractedsites, "w") as of:
                        of.append(print(line))
                        of.close()


##### these would be angsd fasta output
## might not need this
# def parse_fasta(infasta):
#     with open(infasta) as f:
#         seqDict = []
#         seqNr = 0
#         for line in f.readlines():
#             if line.startswith(">"):
#                 seqDict.append({'header':line.rstrip().lstrip(">"), 'sequence':""})
#                 seqNr += 1
#             else:
#                 seqDict[seqNr - 1]['sequence'] = seqDict[seqNr - 1]['sequence'] + line.rstrip()
#     return seqDict

# def change_headers(headers, seqDict):
#     with open(headers) as h:
#         if len(h.readline().split("\t")) == 2:
#             h.seek(0)
#             for line in h.readlines():
#                 for i in range(0, len(seqDict)):
#                     if seqDict[i]['header'].lstrip(">") == line.split()[0].lstrip(">"):
#                         seqDict[i]['header'] = line.split()[1].rstrip().lstrip(">")
#         else:
#             h.seek(0)
#             new_headers = []
#             for line in h.readlines():
#                 new_headers.append(line.rstrip().lstrip(">"))
#             for i in range(0, len(seqDict)):
#                     seqDict[i]['header'] = new_headers[i]
#     return seqDict

def WriteSeqDictToFasta(seqDict, outfasta, append=False):
    with open(outfasta ,"wt") as f:
        for seq in seqDict:
            f.write('>' + seq['header'] + "\n" + seq['sequence'] + "\n")
        f.close()

##### the rest are for odd samples(?)

# this script takes a translation file (one column with the aligned positions, one column with the 
# translated positions) and gives me a translated allele count file. Resulting file can then be used
# with the 2extract_sites_from_allele_counts.py to extract diagnostic sites.

def translate(allele_count,translated_sites):
    allele_count=open(sys.argv[1],"r")
    for line in allele_count:
            values=line.split(" ")
            with open(translated_sites,"r") as translated_sites:
                    for sites in translated_sites:
                            sites=sites.strip()
                            translation=sites.split("\t")
                            if values[0] == translation[1]:
                                    translated_sites.write(translation[0]+" "+values[1]+" "+values[2] + "\n")
                                    translated_sites.close()

def diag_sites_from_fasta(concatenated_seqs_refs2,concatenated_seqs_refs_diag_sites_only):
    # each record in a fastaiterator as in SeqIO.parse is actually a sample, so we want to
    # loop through all of these, one by one
    #you can create a list of dictionaries (how Axel does it) and append each base to the 
    # sample key
    seqDictDiag=[] # creates an empty list for us to append to later on
    for record in SeqIO.parse(concatenated_seqs_refs2,"fasta"):
        #and for each sample/record loop through all positions
        seqDictDiag.append({'sample':record.name, 'sequence':""}) #this creates a list of
    # dictionaries, with the 'sample' key being set to the record id and the sequence is
    #an empty string so far, which we will appen to in the following loop
        with open(concatenated_seqs_refs_diag_sites_only, "r") as positions: # generally good
    # to use the "with" satement when looping through files like this, makes sure they
    # open and close when they should (it will close once the with-loop has finished
            for x in positions:
                print(record.id,int(x),record.seq[int(x)-1])
                seqDictDiag[-1]['sequence']=seqDictDiag[-1]['sequence'] + record.seq[int(x)-1]

#Input arguments
parser = argparse.ArgumentParser(description="Identify and count sites per sample that are fixed for differential alleles between two predefined populations")
# parser.add_argument('--infasta', type=str, help="Fasta alignment input", required=True)
# parser.add_argument('--headers', type=str, help="This file should have either one column with the new headers only, must then be in exactly the same order as in the infasta. can also be two columns with old name in  first column and new name in second, columns should be separated by a tab.", required=True)
# parser.add_argument('--outfasta', type=str, help="The re-headered fasta file to be written", required=True)
parser.add_argument('-i', '--input', type=str, help='Input fasta alignment.', required=True)
parser.add_argument('--popfile', type=str, help="File with two columns separated by a tab. First column should have sample name and second the population it belongs to. Second column can be empty and then sample will be ignored.")
parser.add_argument('--reference', type=str, help="File with reference's diagnostic sites", required=False)
parser.add_argument('-oc', '--outputcounts', type=str, help="Path to output file containing the diagnostic counts per sample", required=True)
parser.add_argument('--allele_count', type=str, help="",required=True)
parser.add_argument('--translated_sites', type=str, help="",required=True)
parser.add_argument('--concatenated_seqs_refs2', type=str, help="",required=True)
parser.add_argument('--concatenated_seqs_refs_diag_sites_only', type=str, help="",required=True)
parser.add_argument('--seqdicttofasta',type=str,help="",required=True)

args = parser.parse_args()

# first convert fasta to dictionary:
# seqDict = parse_fasta(args.infasta)

#next we rename the headers according to the headerfile
# seqDictRename = change_headers(args.headers, seqDict)

#and write this to a new file:
# WriteSeqDictToFasta(seqDictRename, args.outfasta, append=False)

## handle the input fasta to a dictionary
seqDictAln = FastaFetch(args.input)
print("Found " + str(len(seqDictAln)) + " sequences in the alignment.")
#then run this seqDictAln through the next function along with the popfile and output from input arguments
FindDiagnosticPositions(seqDictAln, args.popfile, args.outputsites)

print("First sequence length: " + str(len(seqDictAln[0]['sequence'])))

#then count all diagnostic sites per sample and output to tsv
CountDiagnosticPositions(seqDictAln, args.reference, args.reference)

ExtractDiagSites(diagsites,allelecounts,extractedsites)

# get the right positions
translate(args.allele_count,args.translated_sites)

# make seq dictionary from diagnostic sites
seqDictDiag = diag_sites_from_fasta(args.concatenated_seqs_refs2,args.concatenated_seqs_refs_diag_sites_only)

# write dictionary to file
WriteSeqDictToFasta(seqDictDiag, seqdicttofasta)
