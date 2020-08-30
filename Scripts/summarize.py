# Import libraries
import argparse
import pandas as pd
from os import listdir

parser =  argparse.ArgumentParser()

parser.add_argument('-f', '--folder', help="Folder input")
parser.add_argument('-p', '--prefix', help="Prefix to remove")
parser.add_argument('-s', '--suffix', help="Suffix to remove")
parser.add_argument('-e', '--extension', help="Extension")
parser.add_argument('-o', '--output', help="Ouput name")

args = parser.parse_args()

folder_parsed = args.folder
prefix_parsed = args.prefix
suffix_parsed = args.suffix  
extension_parsed = args.extension
output_parsed = args.output

files = listdir(args.folder)
inputs = []
for i in range(len(files)):
    if files[i].endswith(args.extension) == True:
        inputs.append(files[i])

names = []
for i in range(len(inputs)):
    name = inputs[i].replace(args.prefix,'')
    name = name.replace(args.suffix,'')
    names.append(name)

def summary_file(file_name,name):
    df = pd.read_csv(file_name,sep=',',skiprows=1,header=None)
    df.columns = ["kmer","Met","Unmet","Total","Perc_Met"]
    kmer_temp = df.loc[:,["kmer","Total","Perc_Met"]]
    met_sum = df.loc[:,"Met"].sum()
    tot_sum = df.loc[:,"Total"].sum()
    per_met = met_sum*100/tot_sum
    kmer_temp.loc["Sum"] = ["Total", tot_sum, per_met]
    kmer_temp["Perc_Nor"] = kmer_temp["Perc_Met"] - per_met
    kmer_temp.columns = ["kmer",str(name)+"_Total",str(name)+"_Perc_met",str(name)+"Perc_Norm"]
    return(kmer_temp)

folder_input=args.folder

output_file = summary_file(str(folder_input)+"/"+str(inputs[0]),names[0])
for i in range(1,len(inputs)):
    temp_output_file = summary_file(str(folder_input)+"/"+str(inputs[i]),names[i])
    output_file = output_file.merge(temp_output_file, left_on='kmer', right_on='kmer')

output_file_t = output_file.T
output_file_t.columns = output_file_t.iloc[0]
output_file_t = output_file_t.drop(output_file_t.index[0])

output_file_total = output_file_t.iloc[list(range(0,len(output_file_t),3)),:]
output_file_total.T.to_csv(str(args.output) + "_total.csv")
output_file_meth = output_file_t.iloc[list(range(1,len(output_file_t),3)),:]
output_file_meth.T.to_csv(str(args.output) + "_meth.csv") 
output_file_norm = output_file_t.iloc[list(range(2,len(output_file_t),3)),:]
output_file_norm.T.to_csv(str(args.output) + "_norm.csv")
