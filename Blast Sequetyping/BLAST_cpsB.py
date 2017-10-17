__autor__ = "florianMauffrey"

import os
import re
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

'''
This script allows blasting cpsB sequences against databases and analyzing the output files to determine serotype.

INSTRUCTIONS
All fasta files must be placed in a "Query" folder. Take care to not put | in file name or fasta header. For both local and online parsing, only the 10 best results are saved. Biopython and Blast+ are required.
'''

class Blast:
    
    def __init__ (self):

        '''List of fasta sequences in Query folder'''
        self.list_fasta = os.listdir("Query")

    def makeLocalDB(self):

        '''Create DB from cpsB sequences from Jin et al. Has to be done only once.'''
        os.system("makeblastdb -in Data/cps_list_JIN.txt -out Data/database/cps_JIN -dbtype nucl")
                  
    def localBlast(self):
        
        '''Blast of local Jin database'''
        nb_file = 0
        for fasta in self.list_fasta:
            '''If the blast is already performed, it is skipped.'''
            try:
                test_file = open("Data/output_blast/Local_output_{}".format(fasta), "r")
                test_file.close()
                nb_file += 1
                pass
            except:
                os.system("blastn -db Data/database/cps_JIN -query Query/{0} -out Data/output_blast/Local_output_{0} -outfmt 5".format(fasta))
                nb_file += 1
            print ("Local blast : {}/{}".format(nb_file, len(self.list_fasta)))

                        
                
    def onlineBlast(self):

        '''Blast on the NCBI database'''
        nb_file = 0
        for fasta in self.list_fasta:
            '''If the blast is already performed, it is skipped.'''
            try:
                test_file = open("Data/output_blast/Online_output_{}".format(fasta), "r")
                test_file.close()
                nb_file += 1
                pass
            except:
                fasta_string = open("Query/{}".format(fasta)).read()
                result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=fasta_string, hitlist_size=10, alignments = 100, descriptions = 100)
                with open("Data/output_blast/Online_output_{}".format(fasta), "w") as save:
                    save.write(result_handle.read())
                result_handle.close()
                nb_file += 1
            print ("Online blast : {}/{}".format(nb_file, len(self.list_fasta)))

    def localParsing(self):

        '''
        Analyze local blast XML output files
        '''

        '''Parsing on the XML files from the local blast.'''
        print ("Local parsing ...")
        with open("Local_BLAST_results", "w") as results:
            for fasta in self.list_fasta:
                results_handle = open("Data/output_blast/Local_output_{}".format(fasta))
                blast_record = NCBIXML.read(results_handle)
                results_handle.close()
                results.write("{}\n".format(fasta))

                '''Only the 10 best alignments are saved'''
                i = 0
                for alignment in blast_record.alignments:
                    hsp = alignment.hsps[0]
                    if i < 10:
                        identity = str(hsp.identities) + "/" + str(hsp.align_length)
                        results.write("{}\t{}\n".format(alignment.title, identity))
                        i += 1
                    else:
                        break
                results.write("\n")
        results.close()
        print ("Local parsing done")



    def onlineParsing(self):

        '''
        Analyze online blast XML output files
        '''

        '''Parsing on the XML files from the online blast.'''
        print ("Online parsing ...")
        with open("Online_BLAST_results", "w") as results:
            for fasta in self.list_fasta:
                results_handle = open("Data/output_blast/Online_output_{}".format(fasta))
                blast_record = NCBIXML.read(results_handle)
                results_handle.close()
                results.write("{}\n".format(fasta))

                '''Only the 10 best alignments are saved'''
                i = 0
                for alignment in blast_record.alignments:
                    hsp = alignment.hsps[0]
                    if i < 10:
                        identity = str(hsp.identities) + "/" + str(hsp.align_length)
                        results.write("{}\t{}\n".format(alignment.title, identity))
                        i += 1
                    else:
                        break
                results.write("\n")
        results.close()
        print ("Online parsing done")



'''
Script execution
'''      
option = 0
while option != (1 or 2 or 3 or 4 or 5):
    option = input("\nOptions\n1 - Local Blast\n2 - Online Blast\n3 - Local and Online Blast\n4 - Clean Output_blast\n5 - Create local database\n6 - Exit\n\n")
    if option == 1:
        os.system("rm -rf Data/Temp/*")
        script = Blast()
        script.localBlast()
        script.localParsing()
    elif option == 2:
        os.system("rm -rf Data/Temp/*")
        script = Blast()
        script.onlineBlast()
        script.onlineParsing()
    elif option == 3:
        os.system("rm -rf Data/Temp/*")
        script = Blast()
        script.localBlast()
        script.parsing_local()
        script.onlineBlast()
        script.onlineParsing()
    elif option == 4:
        os.system("rm -rf Data/Output_blast/*")
    elif option == 5:
        script = Blast()
        script.makeLocalDB()
    elif option == 6:
        break
    else:
        print ("Invalid")    
