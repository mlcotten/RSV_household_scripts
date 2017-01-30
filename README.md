Scripts used in publication "Transmission patterns and evolution of respiratory syncytial virus in a community outbreak identified by genomic analysis" Virus Evolution 2017 (in press). Contact Matthew Cotten  (m.cotten@erasmusmc.nl) for details.

1. HiLiter_RSVHH.py generates plots of nucleotide changes in RSV genomes occuring across households. It takes as input a list of nucleotide alignments in fasta format and outputs a pdf of the figure (see Figure 2 in the manuscript or example commands and outputs in the alignment_lists_examples.zip directory). The set of 8 hoiusehold alignments is also provided (eight_RSV_HH_alignments.zip). Typical command is python HiLiter_RSVHH.py hh_alignment_list_normal.txt normal

2. cartman_einfach.py is a simple nucleotide motif counting script using ack, a faster version of grep (http://beyondgrep.com/why-ack/). It takes as inputs a table of motifs to count (usually 21 nt with the variant nucleotide at nucleotide 11) and a list of paired fastq read files to examine. The paired reads are combined, examined for motif content and the numbers of occurences of the forward and reverse complement of the motif in the combined fastq file are reported in a csv table. Typical command is python cartman_einfach.py test_reads_10K.txt test_snp_table.csv 
