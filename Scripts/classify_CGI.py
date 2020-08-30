# Import libraries
import argparse
import pybedtools
import pandas as pd
import subprocess
from os import listdir

parser =  argparse.ArgumentParser()

parser.add_argument('-b', '--bed', help="BED file")
parser.add_argument('-i', '--input', help="Input filename")
parser.add_argument('-o', '--output', help="Output filename")

args = parser.parse_args()

bed_parsed = args.bed
filename_parsed = args.input
output_parsed = args.output

# Import BED file
BED_CGI = pybedtools.BedTool(bed_parsed)

# Obtain CGI and non CGI calls
df = pd.read_csv(filename_parsed,sep='\t',skiprows=1,header=None)
list_positions = pd.DataFrame(df[[2,3,3,4]]).values.tolist()
temp_bed = pybedtools.BedTool(list_positions)
temp_CGI = subprocess.call(temp_bed.intersect(BED_CGI))
temp_non_CGI = (temp_bed - temp_CGI)

# Save results
temp_CGI.saveas("CGI_" + output_parsed + ".csv")
temp_non_CGI.saveas("non_CGI_" + output_parsed + ".csv")
