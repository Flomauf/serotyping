__author__ = 'ericfournier'

import re
import sys
import yaml
import os
from Bio import SeqIO
from Bio.SeqUtils import GC


class CovStat():
    """
    This class prepare recursively (accross all specimen subdir) a contig coverage file for histogram visualisation in R
    """

    def __init__(self,yamlFile,project,run,assembler,lengthT,specList):

        self.myYamFile=""

        '''NGS project name'''
        self.projectDir=project+"/"

        '''run name'''
        self.runDir="run_"+run+"/"

        '''assembler directory'''
        self.assemDir=assembler+"/"

        '''length threshold'''
        self.lengthT=lengthT

        self.SpecList=specList

        '''For each specimen, generate length vs coverage file'''
        map(self.ParseContigFile,self.SpecList)


    def ParseYamlFile(self):
        """
        parse the yaml file
        """

        '''number of specimen'''
        self.NbSpec=len(SpecList)

    def ParseContigFile(self,spec):
        """
        Generate coverage vs contigs length file from assembler fasta files
        """
        '''Remove spaces in specimen names'''
        spec=re.sub(" ","",spec)
        spec=str(spec)


        '''Output file'''
        with open("Sequences/{0}/Filtered/{0}_KeepOver_500pb_and_5x_stat.txt".format(spec),'w') as writef:

            '''header'''
            writef.write("NODE\tCoverage\tLength\n")

            try:
                '''Read assembling fasta file'''
                with open("Sequences/{0}/Filtered/{0}_KeepOver_500pb_and_5x.fasta".format(spec)) as readf:

                    for line in readf:
                            '''contig header'''
                            if re.search(r'>',line):
                                '''extract name, length and coverage'''
                                param=self.ParseHeader(line)
                                node=param[0]
                                length=param[1]
                                coverage=param[2]

                                '''save in output file if contig length is greater than the threshold'''
                                writef.write(node+"\t"+coverage+"\t"+length+"\n")

                readf.close()
            except IOError:
                print "Can't find file ", spec,".fasta\n"
        writef.close()


    def ParseHeader(self,line):
        """
        Extract name, length and coverage
        """
        HeaderPattern=re.search(r'NODE_(\d+)_length_(\d+)_cov_(\d+)',line)

        param=[]

        param.append(HeaderPattern.group(1))
        param.append(HeaderPattern.group(2))
        param.append(HeaderPattern.group(3))

        return  param

    def ComputeGlobalStat(self,spec):

        stat_list=[]

        self.SpecList = [spec]

        self.stat_file = open("Sequences/{0}/Filtered/{0}_GlobalStat.txt".format(spec), 'w')
        self.stat_file.write("Isolate\tTotal length\tPercent GC\tMean coverage\n")

        for spec in self.SpecList:
            stat_list.append(spec)
            stat_list.append(str(self.ComputeTotalLength(spec)))
            stat_list.append(str(self.ComputeGC(spec)))
            stat_list.append(str(self.ComputeMeanCov(spec)))
            stat_list.append('\n')


        self.stat_file.write('\t'.join(stat_list))

        self.stat_file.close()

    def ComputeTotalLength(self,spec):
        tot_length=0

        stat_file =open("Sequences/{0}/Filtered/{0}_KeepOver_500pb_and_5x_stat.txt".format(spec))

        for line in stat_file:
            if not re.search(r'^NODE',line):
                try:
                    tot_length+=int(re.search(r'^\d+\t\d+\t(\d+)$',line).group(1))
                except:
                    print "Prob with ",spec

        stat_file.close()

        return tot_length


    def ComputeGC(self,spec):
        tot_gc=0
        nb_contig=0
        mean_gc=0

        '''path to filtered contig file'''
        for myrec in SeqIO.parse("Sequences/{0}/Filtered/{0}_KeepOver_500pb_and_5x.fasta".format(spec),'fasta'):
            nb_contig+=1
            tot_gc+=GC(myrec.seq)

        mean_gc=round(float(tot_gc)/nb_contig,0)

        return mean_gc

    def ComputeMeanCov(self,spec):
        total_cov=0
        nb_contig=0
        mean_cov=0

        '''path to stat file from ParseContigFile'''
        stat_file = open("Sequences/{0}/Filtered/{0}_KeepOver_500pb_and_5x_stat.txt".format(spec))

        for line in stat_file:
            if not re.search(r'^NODE',line):
                nb_contig+=1
                try:
                    total_cov+=int(re.search(r'^\d+\t(\d+)\t\d+$',line).group(1))
                except:
                    print "Prob with ",spec

        stat_file.close()

        mean_cov=round(float(total_cov)/nb_contig,0)

        return mean_cov


if __name__ == "__main__":


    #the input yaml file
    ListSpecYaml=""

    project=""

    run=""

    assembler=""

    #contig length threshold
    LengthThreshold=1000

    #specimen list
    SpecList = open("SpecList", "r").read().split("\n")
    SpecList.remove("")

    myCovStat=CovStat(ListSpecYaml,project,run,assembler,LengthThreshold,SpecList)

    for spec in  SpecList:
        myCovStat.ComputeGlobalStat(spec)

    #information are gathered in a single file and all filtered fasta file are grouped in a single directory
    with open("Output/GlobalStats", "w") as output:
        for spec in SpecList:
            stat = open("Sequences/{0}/Filtered/{0}_GlobalStat.txt".format(spec), "r")
            output.write(stat.read() + "\n")
    os.mkdir("Output/Sequence_filtered")
    for spec in SpecList:
        os.system("cp Sequences/{0}/Filtered/{0}_KeepOver_500pb_and_5x.fasta Output/Sequence_filtered/{0}_filtered.fasta".format(spec))


