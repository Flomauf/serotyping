__autor__ = "florianMauffrey"


import os
import sys
import re
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

'''
This scripts allows blasting whole genome assemblies against a local cps database and analyzing output files to determine the serotype.

INSTRUCTIONS
All draft genome fasta files must be placed in the "Query" folder. The cps list was made from the list in Camargo et al. 2015.
This script will create a database from each fasta files and blast cps sequences on these databases. Parsing of the XML files will extract best candidates for serotype identifition based on HSP length. Require Biopython and Blast+.
'''

class Blast:

    def __init__ (self):

        '''Creation of queries list'''
        self.filesList = os.listdir("Query")

        '''Databases creation from queries'''
        print ("Creating databases ...")
        for file in self.filesList:
            os.system("makeblastdb -in Query/{0} -out Data/temp/{0} -dbtype nucl > /dev/null".format(file))

    def local_blast(self):

        '''
        Blast of each cps sequence on databases
        '''
        nb_file = 0

        for file in self.filesList:
            '''test if blast already done. If not, perform blast.'''
            try:
                test_file = open("Data/output_blast/Local_output_{}".format(file), "r")
                test_file.close()
                nb_file += 1
                pass
            except:
                os.system("blastn -db Data/temp/{0} -query Data/Liste_cps.fasta -out Data/output_blast/Local_output_{0} -outfmt 5 -max_hsps 3 -num_alignments 5".format(file))
                nb_file += 1
            print ("Local blast : {}/{}".format(nb_file, len(self.filesList)))


    def parsing_local(self):

        '''
        Analysis of XML output files
        '''

        print ("Local parsing ...")
        with open("Local_BLAST_results", "w") as results:
            for file in self.filesList:
                results.write("{}\nQuery\tHit\tNucleotides_identity\tIdentity(%)\tSerotype\tScore\n".format(file))
                results_handle = open("Data/output_blast/Local_output_{}".format(file))
                
                '''Parsing of the handle'''
                blast_records = NCBIXML.parse(results_handle)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:

                            '''Filtering results with an minimum HSP of 5000 pb'''
                            if hsp.align_length > 5000 and hsp.align_length < 30000: 

                                '''Extraction of identity length and alignment length'''
                                identity = float(hsp.identities)
                                length = float(hsp.align_length)

                                '''Calculation of the % of identity'''
                                perc_ident = round((identity/length)*100, 2)
                                score = int(hsp.score)
                                try:
                                    serotype = re.search(r'(?<=serotype_)[0-9]{1,2}[A-Za-z]{0,1}', record.query).group()
                                except:
                                    serotype="N/A"

                                '''Saving results'''
                                results.write("{}\t{}\t{}/{}\t{}\t{}\t{}\n".format(record.query, alignment.title, str(hsp.identities), str(hsp.align_length), str(perc_ident), serotype, str(hsp.score)))
            results.write("\n")
        print ("Local parsing done")


#Script execution
os.system("rm -rf Data/temp/*")
programme = Blast()
programme.local_blast()
programme.parsing_local()
