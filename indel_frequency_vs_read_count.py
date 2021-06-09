import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.stats import pearsonr
import os
import adding_percent_summary_file
import re
import argparse
import pprint

def determine_indel_frequency(args):
	determined_indel_freq = {}
	with open(args.summary_indels_percent,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			if line.strip()[0] == '@':
				oligo = eval(line.strip().split('@@@Cpf1_Oligo')[1])
				determined_indel_freq[oligo] = 100
				continue
			indel, oligo_indel, count, perc = line.strip().split('\t')
			if indel == '-':
				determined_indel_freq[oligo] -= eval(perc)*100
	return determined_indel_freq

def find_reported_indel_freq_and_count(args):
	reported_indel_freq = {}	
	guide_count = 	{}
	with open(args.all_guides_freq,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			guide_no, guide, count, perc = line.strip().split('\t')
			reported_indel_freq[eval(guide_no)] = eval(perc)
			guide_count[eval(guide_no)] = eval(count)
	return reported_indel_freq, guide_count
			
def generate_lists(determined_indel_freq, reported_indel_freq, guide_count):
	determined_frequencies = [determined_indel_freq[x] for x in range(1,1252)]
	reported_frequencies = [reported_indel_freq[x] for x in range(1,1252)]
	read_count = [guide_count[x] for x in range(1,1252)]
	return determined_frequencies, reported_frequencies, read_count
	
def plot_reported_and_determined_against_read_count(read_count_indel_freq, determined_frequencies, reported_frequencies, read_count):
	fig,ax = plt.subplots(figsize = (12,12))
	plt.ylim(0,100)
	plt.xscale('log')
	plt.scatter(read_count,reported_frequencies, color = 'cornflowerblue', label = 'Reported indel frequency')
	plt.scatter(read_count, determined_frequencies, color = 'slateblue', label = 'Determined indel frequency')
	plt.xscale('log')

	plt.xlabel('Read count')
	plt.ylabel('Indel frequency')
	for axis in [ax.xaxis]:
		axis.set_major_formatter(ScalarFormatter())
	plt.legend()
	plt.savefig(read_count_indel_freq)
	plt.close()
	
def plot_reported_against_determined(reported_against_determined,determined_frequencies, reported_frequencies):
	fig,ax = plt.subplots()
	plt.scatter(reported_frequencies, determined_frequencies)
	plt.xlabel('Reported indel frequency')
	plt.ylabel('Determined indel frequency')
	plt.axis((0,100,0,100))
	plt.savefig(reported_against_determined)
	print pearsonr(reported_frequencies, determined_frequencies)
	plt.close()

def split_by_activity(args, filter = 100, no_groups = 5):
	if filter is None:
		determined_indel_freq = determine_indel_frequency(args)
		sorted_by_activity = sorted([(key,val) for key,val in determined_indel_freq.iteritems()], key =lambda x:x[1])
	else:
		determined_indel_freq = determine_indel_frequency(args)
		guide_count = adding_percent_summary_file.construct_oligo_indel_dict(args)
		sorted_by_activity = sorted([(key,val) for key,val in determined_indel_freq.iteritems() if guide_count[key][0]>=filter], key =lambda x:x[1])
	
	activity_groups = np.array_split(sorted_by_activity, no_groups)
	return activity_groups

	
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('base_dir', help = 'Path to directory for output files')
	parser.add_argument('-in','--summary_indels_percent', help = 'Path to summary indels percent file',  default = 'C:\\Users\\jt20\\Documents\\Summer17\\Final_results\\summary_file_percent.txt')
	parser.add_argument('--all_guides_freq', help = 'Path to file giving the guide number, guide, read count and reported percentage of indel frequency', default = 'C:\Users\\jt20\\Documents\\Summer17\\Final_results\\core_files\\all_guides_freq.txt')
	args = parser.parse_args()
	
	
	reported_against_determined = os.path.join(args.base_dir, 'Figures','reported_against_determined.png')
	read_count_indel_freq = os.path.join(args.base_dir, 'Figures', 'scatterplot_indel_freq_read_count.png')
	
	determined_indel_freq = determine_indel_frequency(args)
	reported_indel_freq, guide_count = find_reported_indel_freq_and_count(args)
	determined_frequencies, reported_frequencies, read_count = generate_lists(determined_indel_freq, reported_indel_freq, guide_count)
	
	plot_reported_and_determined_against_read_count(read_count_indel_freq, determined_frequencies,reported_frequencies, read_count)
	plot_reported_against_determined(reported_against_determined, determined_frequencies, reported_frequencies)
