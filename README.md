# Serotyping
Scripts for serotyping S. pneumoniae.

These scripts are used to determine S. pneumoniae serotype from DNA data.

Blast Sequetyping

This script performs blast of cpsB sequences on local cpsB database or online NCBI database. All fasta files must be placed in a "Query" directory located in the "Blast Sequetyping" directory. Then just run the BLAST_cpsB.py script.

Requirements:
Python v3.5
Biopython v1.70
Blast+ v2.2.18+

Blast WGS

This script performs blast of WGS assemblies on local cps database. All fasta files must be placed in a "Query" directory located in the "Blast WGS" directory. Then just run the BLAST_WGS.py script.

Requirements:
Python v3.5
Biopython v1.70
Blast+ v2.2.18+

ContigFilter_WGS

This script allows filtrating contigs in fasta files resulting from a SPAdes assembly. All fasta files must be placed in a "Sequences" directory located in the "ContigFilter_WGS" directory.

Requirements:
Python v3.5
Biopython v1.70
Yaml Python package
Quast v4.5
