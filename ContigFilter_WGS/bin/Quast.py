__author__ = 'florianmauffrey'

import os
import re

'''Iterate Quast analysis with filtered fasta files as input. Extract information from Quast reports in a single output file. Please specify the path to quast.py before running'''

SpecList = os.listdir("Output/Sequence_filtered")

for spec in SpecList:
    '''Path to Quast script'''
    os.system("/home/florian/Documents/App/WGS_tools/QUAST/quast.py Output/Sequence_filtered/{0}  -o Output/QUAST/{0}_output".format(spec))

'''Extraction of information from report.txt'''
with open("Output/Stats_quast.txt", "w") as output:
    output.write("Specimen\tAssembly's_length\tLargest_contig\tN50\n")
    for spec in SpecList:
        report = open("Output/QUAST/{}_output/report.txt".format(spec), "r").read()
        assembly_length = re.search(r'Total length\s+[0-9]{1,10}', report).group()
        assembly_length_value = re.search(r'[0-9]{1,10}$', assembly_length).group()
        largest_contig = re.search(r'Largest contig\s+[0-9]{1,10}', report).group()
        largest_contig_value = re.search(r'[0-9]{1,10}$', largest_contig).group()
        N50 = re.search(r'N50\s+[0-9]{1,10}', report).group()
        N50_value = re.search(r'[0-9]{1,10}$', N50).group()
        output.write("{}\t{}\t{}\t{}\n".format(spec, assembly_length_value, largest_contig_value, N50_value))
