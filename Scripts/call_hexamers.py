# Import libraries
import argparse
import pandas as pd
from Bio import SeqIO

parser =  argparse.ArgumentParser()

parser.add_argument('-g', '--genome', help="Reference genome")
parser.add_argument('-f', '--filename', help="Filename")
parser.add_argument('-o', '--output', help="Output folder")

args = parser.parse_args()

genome_parsed = args.genome
filename_path_parsed = args.filename
filename_parsed = filename_path_parsed.split("/")[1]  
folder_output_parsed = args.output

# Import reference genome
my_seqlist = []
for seq_record in SeqIO.parse(genome_parsed, 'fasta'):
    my_seqlist.append(seq_record)

# Test reference genome
print("Number of sequences loaded:")
print(len(my_seqlist))
print("Sequences loaded succesfully\n")

# Obtain CG kmers list
from itertools import product

kmers = list(product('ATCG', repeat=6))
kmers = ["".join(x) for x in kmers]

kmers_CG = []
for i in kmers:
    if i[2:4] == "CG":
        kmers_CG.append(i)

# Test CG-kmer list
print("Number of CG K-mers loaded:")
print(len(kmers_CG))
print("CG K-mers loaded succesfully\n")

# Dictionary key location
dict_chr_loc = {}
for i in range(len(my_seqlist)):
    dict_chr_loc[my_seqlist[i].id] = i

# K-mer analysis
def kmer_analysis():    
    file_input = filename_path_parsed
    print(file_input)
    df = pd.read_csv(file_input,sep='\t',skiprows=1,header=None)

    pd_kmers_CG = pd.DataFrame(0, index=kmers_CG, columns=['Met','Unmet'])

    for i in range(len(df)):
        position = my_seqlist[dict_chr_loc.get(str(df[2][i]))].seq[(df[3][i])-1]
        if position == "C":
            substring = my_seqlist[dict_chr_loc.get(str(df[2][i]))].seq[(df[3][i])-3:df[3][i]+3]
        elif position =="G":
            substring = (my_seqlist[dict_chr_loc.get(df[2][i])].seq[(df[3][i])-4:df[3][i]+2]).reverse_complement()
        if (substring[2:4] == "CG") & (str(substring) in kmers_CG):
            if df[4][i] == "Z":
                pd_kmers_CG.loc[str(substring)][0] = (pd_kmers_CG.loc[str(substring)][0]) + 1
            else:
                pd_kmers_CG.loc[str(substring)][1] = (pd_kmers_CG.loc[str(substring)][1]) + 1

    pd_kmers_CG['Total'] = pd_kmers_CG.sum(axis=1)
    pd_kmers_CG['Per_met'] = pd_kmers_CG['Met']*100/pd_kmers_CG['Total']
    file_output = str(folder_output_parsed) + "/" + str(filename_parsed.split(".")[0]) + ".csv"
    print(file_output)
    pd_kmers_CG.to_csv(file_output)

kmer_analysis()
