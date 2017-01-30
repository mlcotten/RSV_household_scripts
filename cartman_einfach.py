#!/usr/local/bin/python
from __future__ import division
import sys
import csv
import os
import subprocess
import glob
path_to_ack= "/Users/mc13/Scripts_to_be_called/ack.pl" #modify with local path to ack.pl

# Takes table of paired fastq files, table of polymorphic sites, counts reads containing each polymorphic site.
# written by M.Cotten WTSI and EMC

print 'Running cartman_einfach.py'	
if len(sys.argv)!=3:
	print "Usage: python cartman.py input_read_files.txt snps.txt"
	sys.exit()
input_read_files = sys.argv[1]
outprefix = os.path.splitext(input_read_files)[0]
snp_table = sys.argv[2]

#Set of functions to generate rev_com
def reverse(s): 
	letters = list(s) 
	letters.reverse() 
	return ''.join(letters) 
def complement(s): 
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y', 'K':'M','M':'K', 'S':'S', 'W':'W', 'N':'N'} 
	letters = list(s) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters) 
def rev_com(s): 
	s = reverse(s) 
	s = complement(s) 
	return s

#create list of reads from table
combined_read_list=[]
read_reader = csv.reader(open(input_read_files, "rU"), delimiter=',')
for row in read_reader:
	reads_1 = str(row[0])
	outprefix_reads = os.path.splitext(reads_1)[0]
	reads_2 = str(row[1])
	combined_file_name = outprefix_reads+"_zusammen.fq"
	cat_call = "cat "+reads_1+" "+reads_2+" > "+outprefix_reads+"_zusammen.fq"
	subprocess.call(cat_call, shell=True)	
	combined_read_list.append(combined_file_name)
	
snp_counts_file = open(outprefix+'_snp_counts_list.csv','a')
print_string_header= "reads_file, position,snpA,snpA_nt,snpA_count,snpB,snpB_nt,snpB_count"
print >> snp_counts_file,print_string_header

for combined_read_file in combined_read_list:
	snp_reader = csv.reader(open(snp_table, "rU"), delimiter=',')
	snp_reader.next() # Skip header line.
	for row in snp_reader:
		position = row[0]
		snpA_seq = row[1]
		snpA_seq_rc = rev_com(snpA_seq)
		snpA_nt =row[2]
		snpB_seq = row[3]
		snpB_nt =row[4]
		snpB_seq_rc = rev_com(snpB_seq)
		snpA_count = subprocess.Popen("perl /Users/mc13/Scripts_to_be_called/ack.pl -c "+snpA_seq+" "+combined_read_file, stdout=subprocess.PIPE, shell=True).communicate()[0].strip()
		snpA_count_rc = subprocess.Popen("perl /Users/mc13/Scripts_to_be_called/ack.pl -c "+snpA_seq_rc+" "+combined_read_file, stdout=subprocess.PIPE, shell=True).communicate()[0].strip()
		snpB_count = subprocess.Popen("perl /Users/mc13/Scripts_to_be_called/ack.pl -c "+snpB_seq+" "+combined_read_file, stdout=subprocess.PIPE, shell=True).communicate()[0].strip()
		snpB_count_rc = subprocess.Popen("perl /Users/mc13/Scripts_to_be_called/ack.pl -c "+snpB_seq_rc+" "+combined_read_file, stdout=subprocess.PIPE, shell=True).communicate()[0].strip()
		total_snpA_counts= int(snpA_count)+int(snpA_count_rc)
		total_snpB_counts= int(snpB_count)+int(snpB_count_rc)
		snp_counts_file = open(outprefix+'_snp_counts_list.csv','a')
		print_string_snp=str(combined_read_file)+","+str(position)+","+str(snpA_seq)+","+str(snpA_nt)+","+str(total_snpA_counts)+","+str(snpB_seq)+","+str(snpB_nt)+","+str(total_snpB_counts)
		print >> snp_counts_file,print_string_snp
		snp_counts_file.close()

combined_files_to_delete = glob.glob('*zusammen*')
for file in combined_files_to_delete: 
		os.remove(file)
		
print 'That\'s All Folks!'