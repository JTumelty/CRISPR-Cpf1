import matplotlib.pyplot as plt
import numpy as np
import indel_frequency_vs_read_count
import mean_read_number
import argparse
import os
def compute_indel_length(indel):
	if indel == '-':
		return 0
	elif indel[0] == 'D':
		del_length = eval(indel.split('_')[0][1:None])
		if indel.find('I') == -1:
			return 0 - del_length
		else:
			for j in range(indel.find('I')+1,len(indel)):
				if indel[j] in ['C','R']:
					ins_length = eval(indel[indel.find('I')+1:j])
					break
			return ins_length - del_length
	elif indel[0] == 'I':
		ins_length = eval(indel.split('_')[0][1:None])
		if indel.find('D') == -1:
			return ins_length
		else:
			for j in range(indel.find('D')+1,len(indel)):
				if indel[j] in ['C','R']:
					del_length = eval(indel[indel.find('D')+1:j])
					break
			return ins_length - del_length

def find_all_indel_lengths(args, mean_reads, guides = range(1,1252)):
	normalised_indel_length = {}
	all_indel_length = {}
	with open(args.summary_indels_percent,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			elif line.strip()[0] == '@':
				oligo_no = eval(line.strip().split('@@@Cpf1_Oligo')[1])
				continue
			if oligo_no in guides:
				indel, oligo_indel, count, perc = line.strip().split('\t')
				indel_length = compute_indel_length(indel)
				if indel_length in normalised_indel_length:
					normalised_indel_length[indel_length] += round(eval(perc)*mean_reads)
					all_indel_length[indel_length] += eval(count)
				else:
					normalised_indel_length[indel_length] = round(eval(perc)*mean_reads)
					all_indel_length[indel_length] = eval(count)
	print normalised_indel_length[0]
	return normalised_indel_length, all_indel_length
	
def plot_indel_lengths(normalised_indel_length, fig_name_1, fig_name_2, y_label, y_lim, both = True):
	indel_length_keys = range(min(normalised_indel_length.keys()), max(normalised_indel_length.keys()))
	indel_length_counts = []
	ind = np.arange(len(indel_length_keys))
	for indel_length in indel_length_keys:
		if indel_length in normalised_indel_length:
			indel_length_counts.append(normalised_indel_length[indel_length])
		else:
			indel_length_counts.append(0)
			
	fig, ax = plt.subplots(figsize = (15,15))
	plt.bar(ind, indel_length_counts)
	plt.xticks(ind, indel_length_keys)
	plt.xlabel('Indel length')
	plt.ylabel(y_label)
	plt.ylim((0,y_lim))
	plt.savefig(fig_name_2)
	plt.close()
	
	if both:		
		fig, ax = plt.subplots(figsize = (15,15))
		plt.bar(ind, indel_length_counts)
		plt.xticks(ind, indel_length_keys)
		plt.xlabel('Indel length')
		plt.ylabel(y_label)
		plt.savefig(fig_name_1)
		plt.close()
			
def plot_percentages(normalised_indel_length, activity_groups):
	indel_length_keys = range(min(normalised_indel_length.keys()), max(normalised_indel_length.keys()))
	indel_length_counts = []
	ind = np.arange(len(indel_length_keys))
	for indel_length in indel_length_keys:
		if indel_length in normalised_indel_length:
			indel_length_counts.append(normalised_indel_length[indel_length])
		else:
			indel_length_counts.append(0)
			
	fig,ax = plt.subplots(figsize = (15,15))
	count = 0
	for activity_group in activity_groups:
		print activity_group[0], activity_group[-1]
		count += 1
		sep_normalised_indel_length = find_all_indel_lengths(args, mean_reads, guides = activity_group)[0]
		sep_indel_length_counts = []
		for indel_length in indel_length_keys:
			if indel_length in sep_normalised_indel_length:
				sep_indel_length_counts.append(sep_normalised_indel_length[indel_length])
			else:
				sep_indel_length_counts.append(0)
		
		perc = []		
		for i in range(len(sep_indel_length_counts)):
			try:
				perc.append(sep_indel_length_counts[i]/float(indel_length_counts[i])*100)
			except ZeroDivisionError:
				perc.append(0)
		plt.scatter(ind, perc, label = 'Activity group ' + str(count))
	
	plt.xticks(ind, indel_length_keys)
	plt.xlabel('Indel length')
	plt.ylim((0,100))
	plt.ylabel('Percentage of indels of that length')
	plt.legend()
	plt.savefig(os.path.join(args.base_dir, 'Figures', 'indel_length_split_activity.png'))

def make_dirs(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
		
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Creates a profile of indel lengths across all guides')
	parser.add_argument('base_dir', help = 'Base directory to store output files')
	parser.add_argument('summary_indels_percent', help = 'Path to summary indels percent file')
	parser.add_argument('--filter', help = 'Include only guides with a read count greater than this value. Default = 100.', default = 100, type = int)
	parser.add_argument('--no_groups', help = 'Split the guides up by activity into this number of groups to produce separate indel length plots for each. 0 (default) to not produce these plots.', default = 0, type = int)
	args = parser.parse_args()
	
	activity_groups = indel_frequency_vs_read_count.split_by_activity(args, filter = args.filter, no_groups = 1)[0]
	mean_reads,std = mean_read_number.mean_read_count(args)
	normalised_indel_length, all_indel_length = find_all_indel_lengths(args,mean_reads, guides = activity_groups)
	plot_indel_lengths(normalised_indel_length,os.path.join(args.base_dir,'Figures','indel_lengths_all_height.png'),os.path.join(args.base_dir,'Figures', 'indel_lengths_truncated.png'), 'Normalised number of reads', 25000)
	
	if args.no_groups != 0:
		fig_sub_dir = os.path.join(args.base_dir, 'Figures','indel_lengths_truncated')
		make_dirs(fig_sub_dir)
		count = 0
		activity_groups = indel_frequency_vs_read_count.split_by_activity(args,filter = args.filter, no_groups = args.no_groups)
		plot_percentages(normalised_indel_length, activity_groups)
		for activity_group in activity_groups:
			separated_normalised_indel_length = find_all_indel_lengths(args, mean_reads, guides = activity_group)[0]
			plot_indel_lengths(separated_normalised_indel_length,os.path.join(fig_sub_dir, ''.join(['indel_lengths_all_height_group_',str(count+1),'.png'])),os.path.join(fig_sub_dir, ''.join(['indel_lengths_truncated_group_',str(count+1),'.png'])), 'Normalised number of reads', 6000, both = False)
			count += 1