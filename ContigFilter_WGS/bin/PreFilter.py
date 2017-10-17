import os
import re
from Bio import SeqIO

'''
Filter contigs according to length and coverage. 
Make an Input and an Output folder. Put fasta sequence in the Input folder.
'''

#List of specimen
SpecList=['MA080904']

#Coverage and length parameters
cov_min = 5
cov_max = 1000
len_min = 500
len_max = 10000000

#Filtering according to previous parameters
for spec in SpecList:
    fasta_sequences = SeqIO.parse("Input/{}.fasta".format(spec),'fasta')
    with open("Output/{}_filtered_cov[{}-{}]_length[{}-{}].fasta".format(spec, str(cov_min), str(cov_max), str(len_min), str(len_max)), "w") as output:
        for fasta in fasta_sequences:
            ID, sequence = str(fasta.id), str(fasta.seq) 
            Cov = float(re.search(r'(?<=cov_)[0-9]{1,4}', ID).group())
            Len = float(re.search(r'(?<=length_)[0-9]{1,10}', ID).group())
  
            if Cov <= cov_max and Cov >= cov_min and Len <= len_max and Len >= len_min:
                output.write(">{}\n{}\n".format(ID, sequence))
