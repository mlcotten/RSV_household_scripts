Scripts used in publication "Transmission patterns and evolution of respiratory syncytial virus in a community outbreak identified by genomic analysis" by Charles N. Agoti et al. Virus Evolution 2017 (in press).
email Matthew Cotten  (m.cotten@erasmusmc.nl) for details.

1. hiliter.py shows nucleotide changes in RSV genomes occuring across households. It takes as input a list of nucleotide alignments in fasta format and outputs a pdf of the figure. 

2. cartman_einfach.py is a simple grep script using the much faster ack (http://beyondgrep.com/why-ack/). It takes as input a table of 21 nt motifs to count and a list of fastq read files to examine. Number of reads in the fastq file bearing the forward or reverse complement of the motif are reported.
