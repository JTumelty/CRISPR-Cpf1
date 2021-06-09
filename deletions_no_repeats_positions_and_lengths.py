import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare
import mean_read_number
import indel_frequency_vs_read_count
import argparse
import os

def count_guides(args):
	guide_count = {}
	with open(args.guides_file,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			guide_no, guiide, count = line.strip().split('\t')
			guide_count[guide_no] = eval(count)
	return guide_count

def parse_indels_associated_condensed(indels_associated_condensed_file,guide_count, mean_reads, guides = range(1,1252)):
	left_positions = {}
	right_positions = {}
	deletion_lengths  = {}
	with open(indels_associated_condensed_file ,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			guide_no,perc, indel, rep, del_len, start,stop = line.strip().split('\t')
			if indel[0] != 'D' or indel.find('I')!= -1 or indel.find('C') != -1 or rep != '-':
				continue
			if eval(guide_no) in guides:
				update_dicts(guide_count, guide_no, indel, deletion_lengths, left_positions, right_positions, perc, mean_reads)	
	return deletion_lengths, left_positions, right_positions
	
def update_dicts(guide_count, guide_no, indel, deletion_lengths, left_positions, right_positions, perc, mean_reads):
	perc = eval(perc)
	del_length = eval(indel.split('_')[0][1:None])
	positions = re.split('[LCIRD]', indel.split('_')[1])
	left_pos = eval(positions[1])
	right_pos = eval(positions[-1])
	
	deletion_lengths = add_to_dict(del_length, deletion_lengths,round(mean_reads*perc))
	left_positions = add_to_dict(left_pos, left_positions,round(mean_reads*perc))
	right_positions = add_to_dict(right_pos, right_positions,round(mean_reads*perc))

	return deletion_lengths, left_positions, right_positions

def add_to_dict(entry, my_dict, value):
	if entry in my_dict:
		my_dict[entry] += value
	else:
		my_dict[entry] = value
	return my_dict
	
def plot_positions(left_positions, right_positions, file_name_left, file_name_right):
	left_positions_keys = range(min(left_positions.keys()), max(left_positions.keys())+1)
	ind = np.arange(len(left_positions_keys))
	left_position_counts = []
	for left_pos in left_positions_keys:
		if left_pos in left_positions:
			left_position_counts.append(left_positions[left_pos])
		else:
			left_position_counts.append(0)
	ax,fig = plt.subplots(figsize = (15,15))
	plt.bar(ind,left_position_counts)
	plt.xlabel('Left position')
	plt.xticks(ind, [left_positions_key + 19 for left_positions_key in left_positions_keys])
	plt.savefig(file_name_left)
	plt.close()
	
	right_positions_keys = range(min(right_positions.keys()), max(right_positions.keys())+1)
	ind = np.arange(len(right_positions_keys))
	ind_cor = [ind[x]+ 19 for x in range(len(ind))]
	right_position_counts = []
	for right_pos in right_positions_keys:
		if right_pos in right_positions:
			right_position_counts.append(right_positions[right_pos])
		else:
			right_position_counts.append(0)
	ax,fig = plt.subplots(figsize = (15,15))
	plt.bar(ind,right_position_counts)
	plt.xlabel('Right position')
	plt.xticks(ind, [right_positions_key + 19 for right_positions_key in right_positions_keys])
	plt.savefig(file_name_right)
	plt.close()
	return left_positions_keys, right_positions_keys, left_position_counts, right_position_counts

def plot_deletion_lengths(deletion_lengths, file_name):
	deletion_length_keys = range(-max(deletion_lengths.keys())-1, -min(deletion_lengths.keys())+1)
	ind = np.arange(len(deletion_length_keys))
	deletion_lengths_counts = []
	for del_len in deletion_length_keys:
		if -del_len in deletion_lengths:
			deletion_lengths_counts.append(deletion_lengths[-del_len])
		else:
			deletion_lengths_counts.append(0)
	ax,fig = plt.subplots(figsize = (15,15))
	plt.bar(ind,deletion_lengths_counts)
	plt.xlabel('Deletion')
	plt.xticks(ind, deletion_length_keys)
	plt.savefig(file_name)
	plt.close()
	return deletion_length_keys, deletion_lengths_counts

def determine_left_right_inderpendence_data(left_positions_keys, right_positions_keys, left_position_counts, right_position_counts):
	cumulative_left = np.cumsum(left_position_counts)
	cumulative_right = np.cumsum(right_position_counts)
	exp_del_lengths = {}
	count = 0
	while count <= cumulative_left[-1]:
		left_rand = np.random.randint(round(cumulative_left[-1]))
		right_rand = np.random.randint(round(cumulative_right[-1]))
		for i in range(len(cumulative_left)):
			if left_rand <= cumulative_left[i]:
				left_pos = left_positions_keys[i]
				break
		for i in range(len(cumulative_right)):
			if right_rand <= cumulative_right[i]:
				right_pos = right_positions_keys[i]
				break
		if left_pos < right_pos-1:
			del_length = -left_pos + right_pos -1
			exp_del_lengths = add_to_dict(del_length, exp_del_lengths, 1)
			count += 1
	return exp_del_lengths
	
def plot_left_right_inderpendence_data(all_del_lengths, obs_deletion_length_keys,obs_deletion_length_counts, file_name):
	deletion_length_keys = [-x for x in sorted(all_del_lengths.keys(), reverse = True)]
	max_del_length = min(deletion_length_keys+ obs_deletion_length_keys)
	all_ind = range(max_del_length,0)
	exp_deletion_length_counts = []
	new_obs_deletion_length_counts = []
	for i in range(len(all_ind)):
		if -all_ind[i] in all_del_lengths:
			exp_deletion_length_counts.append(all_del_lengths[-all_ind[i]])
		else:
			exp_deletion_length_counts.append(0)
		if all_ind[i] in obs_deletion_length_keys:
			new_obs_deletion_length_counts.append(obs_deletion_length_counts[obs_deletion_length_keys.index(all_ind[i])])
		else:
			new_obs_deletion_length_counts.append(0)
	ind = np.arange(len(deletion_length_keys))
	fig,ax = plt.subplots(figsize = (15,15))
	plt.bar([all_ind[i] + 0.1 for i in range(len(all_ind))], exp_deletion_length_counts, width = 0.4, label = 'expected')
	plt.bar([all_ind[i] + 0.5 for i in range(len(all_ind))], new_obs_deletion_length_counts, width = 0.4, label = 'observed')
	plt.xlabel('Deletion')
	plt.legend()
	plt.savefig(file_name)
	plt.close()
	return deletion_length_keys, exp_deletion_length_counts


def make_dirs(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	return dir_name
		
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Plots the lengths of the deletions')
	parser.add_argument('base_dir', help = 'Base directory where the output files are stored')
	parser.add_argument('summary_indels_percent', help = 'Path to summary percent file')
	parser.add_argument('guides_file', help = 'Path to "all_guides.txt" a file with guide number, guide sequence and read count')
	parser.add_argument('--activity_groups', help = 'Number of activity groups to plot the indel profiles for, 0 is default', default = 0, type = int)
	args = parser.parse_args()
	
	indels_associated_condensed_file = os.path.join(args.base_dir, 'indels_associated_condensed.txt')
	fig_dir_name = os.path.join(args.base_dir,'Figures','deletion_positions')
	make_dirs(fig_dir_name)
	
	mean_reads, std = mean_read_number.mean_read_count(args)
	
	guide_count = count_guides(args)
	deletion_lengths, left_positions, right_positions = parse_indels_associated_condensed(indels_associated_condensed_file,guide_count, mean_reads)
	obs_deletion_length_keys, obs_deletion_length_counts = plot_deletion_lengths(deletion_lengths, os.path.join(args.base_dir, 'Figures','deletions_lengths_no_repeats_total_counts.png'))
	left_positions_keys, right_positions_keys, left_position_counts, right_position_counts = plot_positions(left_positions, right_positions, os.path.join(args.base_dir, 'Figures','deletions_left_pos_no_repeats_total_counts.png'), os.path.join(args.base_dir, 'Figures','deletions_right_pos_no_repeats_total_counts.png'))
	exp_del_lengths = determine_left_right_inderpendence_data(left_positions_keys, right_positions_keys, left_position_counts, right_position_counts)
	exp_deletion_length_keys, exp_deletion_length_counts = plot_left_right_inderpendence_data(exp_del_lengths, obs_deletion_length_keys, obs_deletion_length_counts, os.path.join(args.base_dir, 'Figures','exp_del_len_change_left_right_pos_counts.png'))
	
	if args.activity_groups != 0:
		count = 0
		activity_groups = indel_frequency_vs_read_count.split_by_activity(args,no_groups = args.activity_groups )
		for activity_group in activity_groups:
			guide_count = count_guides(args)
			deletion_lengths, left_positions, right_positions = parse_indels_associated_condensed(indels_associated_condensed_file, guide_count, mean_reads, guides= activity_group)
			obs_deletion_length_keys, obs_deletion_length_counts = plot_deletion_lengths(deletion_lengths, os.path.join(fig_dir_name,'deletions_lengths_no_repeats_group_'+str(count + 1)+'.png'))
			left_positions_keys, right_positions_keys, left_position_counts, right_position_counts = plot_positions(left_positions, right_positions, os.path.join(fig_dir_name,'deletions_left_pos_no_repeats_group_'+str(count+1)+'.png'), os.path.join(fig_dir_name,'deletions_right_pos_no_repeats_group_'+str(count+1)+'.png'))
			exp_del_lengths = determine_left_right_inderpendence_data(left_positions_keys, right_positions_keys, left_position_counts, right_position_counts)
			exp_deletion_length_keys, exp_deletion_length_counts = plot_left_right_inderpendence_data(exp_del_lengths, obs_deletion_length_keys, obs_deletion_length_counts, os.path.join(fig_dir_name,'exp_del_len_change_left_right_pos_counts_group_'+str(count+1)+'.png'))
			count += 1