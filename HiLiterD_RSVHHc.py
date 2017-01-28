#!/usr/local/bin/pythonimport sys
sys.path.insert(0, "/Users/matthewcotten/biopython-1.60")# new version of biopython
import os.pathfrom Bio import SeqIO
from Bio.Seq import Seqfrom Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
import sys
import os.path
from Bio.Blast import NCBIXML
from Bio import Entrez
import time
import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#Compares all sequences in aligned fasta file to first sequence in set

if len(sys.argv)!= 5:
	print("Usage: python HiLiter.py alignment_table list_yes_no, features.csv, panel_size")
	sys.exit()

print "Running HiLiter.py really fast"

alignment_table= sys.argv[1]
outprefix = os.path.splitext(alignment_table)[0]
list_yes_no = sys.argv[2]
feature_table = sys.argv[3]
panel_size =sys.argv[4]

#prepare list of alignments to pyplot
alignment_list=[]
starting_csv = open(alignment_table, "rU")
reader1 = csv.reader(starting_csv)
# reader1.next() # Skip header line.
for row in reader1:
	alignment_list.append(row[0])

fig = plt.figure()
plot_index=1
for alignment in alignment_list:
	seq_names =[]
	sequences =[]
	for record in SeqIO.parse(open(alignment, "rU"), "fasta"):
		seq_names.append(record.id)
		sequence_string = str(record.seq)
		sequences.append(sequence_string)
	total_number_sequences = len(sequences)
	genome_length = len(sequences[0])

	#width should be proportional to length of sequence
	#if genome_length <= 500:
	#	width_unit = 0.5
	#elif genome_length >500 and genome_length <= 5000:
	#	width_unit = 5
	#elif genome_length >5000 and genome_length <= 10000:
	#	width_unit =7
	#elif genome_length >10000 and genome_length <= 20000:
	#	width_unit =15
	#if genome_length >20000:
	#	width_unit =30


	if total_number_sequences >18 and total_number_sequences <=21:
		ytick_label_size = 4
	else:
		ytick_label_size = 6

	#genome_length = len(sequences[0])
	#if genome_length >= 0 and genome_length < 150:
	#	marker_size = 1
	#elif genome_length >= 150 and genome_length < 300:
	#	marker_size = 2
	#elif genome_length >= 300 and genome_length < 600:
	#	marker_size = 3
	#elif genome_length >= 600 and genome_length < 1500:
	#	marker_size = 6
	#elif genome_length >= 1500 and genome_length < 3000:
	#	marker_size = 15
	#else:
	#	marker_size = 45
	#graph_pad = 5 * marker_size
	graph_names=[]
	graph_names.append(seq_names[0])

	width_unit = 15
	gap_unit=1
	# ytick_label_size = 6
	marker_size = 45
	graph_pad = 5 * marker_size

	#Generate lists of differences between ref and test sequences
	#function for comparing sequence a to b
	def vergleichen(seqA, seqB, index):
	#	width_unit = 25 #may need to change the default width to 25 for larger genomes
		gap_end_marker = 500000 #need value that would never be encountered.
		gap_start_position = 0
		current_gap_width = 0
		current_gap_marker = 0
		gap_dictionary = {}
		for i in range (len(seqA)):
			position = i
			if seqB[i]!=seqA[i] and seqB[i] == "A":
				color = "orangered"
				line = index
				diff_list.append((position,width_unit,line,color))
				if list_yes_no == "yes":
					difference_table = open(outprefix +"_difference_table.csv", "a")
					print >> difference_table, [i], seqA[i], seqB[i], seqA[i-10:i+11], seqB[i-10:i+11]
					difference_table.close()
			elif seqB[i]!=seqA[i] and seqB[i] == "T":
				color = "crimson"
				line = index
				diff_list.append((position,width_unit,line,color))
				if list_yes_no == "yes":
					difference_table = open(outprefix +"_difference_table.csv", "a")
					print >> difference_table, [i], seqA[i], seqB[i], seqA[i-10:i+11], seqB[i-10:i+11]
					difference_table.close()
			elif seqB[i]!=seqA[i] and seqB[i] == "G":
				color = "indigo"
				line = index
				diff_list.append((position,width_unit,line,color))
				if list_yes_no == "yes":
					difference_table = open(outprefix +"_difference_table.csv", "a")
					print >> difference_table, [i], seqA[i], seqB[i], seqA[i-10:i+11], seqB[i-10:i+11]
					difference_table.close()
			elif seqB[i]!=seqA[i] and seqB[i] == "C":
				color = "slateblue"
				line = index
				diff_list.append((position,width_unit,line,color))
				if list_yes_no == "yes":
					difference_table = open(outprefix +"_difference_table.csv", "a")
					print >> difference_table, [i], seqA[i], seqB[i], seqA[i-10:i+11], seqB[i-10:i+11]
					difference_table.close()
	#		elif seqB[i]!=seqA[i] and seqB[i] == "-": #original version currently this maps each nt gap.
	#			color = "DarkGray"
	#			line = index
	#			diff_list.append((i,line,color))
	#			if list_yes_no == "yes":
	#				difference_table = open(outprefix +"_difference_table.csv", "a")
	#				print >> difference_table, [i], seqA[i], seqB[i], seqA[i-10:i+11], seqB[i-10:i+11]
	#				difference_table.close()

			elif seqB[i] == "N": #currently this maps each nt gap.
				if position-1 != gap_end_marker: #if the previous position is not gap then this is start of new gap
					gap_start_position = position
					gap_end_marker = position
					current_gap_width = gap_unit
					print "position"
					print position
					print "start of a gap"
					color = "#e6e6e6"
					line = index
					print "line"
					print line
					gap_dictionary[gap_start_position]=(current_gap_width,line,color) #move to
				else: #then it must be == to gap end marker and is an extension of a gaps
					current_gap_width = current_gap_width+gap_unit #extend current gap	by 1 unit
					gap_end_marker = position
					print "extend a gap"
					print "gap_start_position"
					print gap_start_position
					print "current_gap_width"
					print current_gap_width

					color = "#e6e6e6"
					line = index
					print "line"
					print line
	#				diff_list.append((gap_start_position,current_gap_width,line,color))
					gap_dictionary[gap_start_position]=(current_gap_width,line,color) #move to gap_dictionary
	#				if list_yes_no == "yes":
	#					difference_table = open(outprefix +"_difference_table.csv", "a")
	#					print >> difference_table, [i], seqA[i], seqB[i], seqA[i-10:i+11], seqB[i-10:i+11]
	#					difference_table.close()
	#			else: #if the previous position not a gap then this is start of new gap
	#				gap_start_position = position
	#				print "start of a gap"

		for k,v in gap_dictionary.iteritems():
			gap_start_position = k
			current_gap_width = v[0]
			line = v[1]
			color = v[2]
			print "dictionary dump"
			print gap_start_position,current_gap_width,line,color
			diff_list.append((gap_start_position,current_gap_width,line,color))
			if list_yes_no == "yes":
				difference_table = open(outprefix +"_difference_table.csv", "a")
				print >> difference_table, gap_start_position,current_gap_width
				difference_table.close()

	diff_list=[]
	#call vergleichen on each pair
	for i in range(1, len(sequences)):
		vergleichen(sequences[0], sequences[i], i)

	#now generate graph of differences.
	ticks=[]
	for x in range(len(sequences)):
		ticks.append(x)
	# fig = plt.figure()
	#fig.subplots_adjust(left=0.25, bottom=None, right=None, wspace=None, top = None)
	# ax2 = fig.add_subplot(4,2, int(plot_index))
	#http://stackoverflow.com/questions/3330137/adjusting-heights-of-individual-subplots-in-matplotlib-in-python
	# ax2 = plt.subplot2grid((8,1), (int(plot_index),0), rowspan=3)
	if panel_size == "normal":
		ax2 = fig.add_subplot(4,2, int(plot_index))
	elif panel_size == "large":
		ax2 = fig.add_subplot(2,2, int(plot_index))

	ax2.set_ylim(0.9, len(sequences))
	ax2.set_xlim(-200, (genome_length+200))
	ax2.set_yticks(ticks)
	ax2.set_yticklabels(seq_names, va='bottom', size = ytick_label_size)
	#ax2.set_xticklabels(genome_positions,size = ytick_label_size)
	ax2.tick_params(axis='x', which='major', labelsize=6)
	ax2.set_xlabel(seq_names[0]+' Genome Position', fontsize=8)

	#for i in range(len(diff_list)):
	for i in range(0, len(diff_list)):

		color = diff_list[i][3]
		ax2.broken_barh([(int(diff_list[i][0]), (diff_list[i][1]))], (int(diff_list[i][2]), 0.85), facecolors=(color),edgecolors=(color))

	plot_index = plot_index+1


#add feature graph at top of each column
# feature_table = csv.reader(open(feature_table, "rU"))
# feature_names =[]
# feature_map_positions = []
# feature_lengths = []
# feature_colors=[]
# feature_table.next()#skip header
# for row in feature_table:
# 	feature_name = str(row[0])
# 	feature_names.append(feature_name)
# 	feature_start = int(row[1])
# 	feature_map_positions.append(feature_start)
# 	feature_length = int(row[2])- int(row[1])
# 	feature_lengths.append(feature_length)
# 	feature_color = str(row[3])
# 	feature_colors.append(feature_color)
#
# ax3 = fig.add_subplot(10,2,20)
#
# ax3.set_ylim(0, 1)
# ax3.set_xlim(-graph_pad,(genome_length+graph_pad))
# ax3.set_xticklabels([])
# ax3.set_yticklabels([])
# ax3.grid(False)
# ax3.get_yaxis().set_visible(False)
# for i in range(len(feature_map_positions)):
# 	ax3.broken_barh([(int(feature_map_positions[i]), int(feature_lengths[i]))],(0.2, 0.3), facecolors=(str(feature_colors[i])), edgecolors=('None'))
# 	ax3.annotate(feature_names[i], xy=(1,1), xytext=(feature_map_positions[i], 0.6), fontsize=4)


plt.tight_layout()
plt.savefig(outprefix+'_differences.pdf')
print "That's All Folks!"
