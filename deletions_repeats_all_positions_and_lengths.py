import re
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
def parse_indels_associated_condensed(indels_associated_condensed_file):
	left_positions = {}
	right_positions = {}
	deletion_lengths  = {}
	deletion_length_and_pos = {}
	with open( indels_associated_condensed_file,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			guide_no,perc, indel, rep, del_len, start,stop = line.strip().split('\t')
			if indel[0] != 'D':
				continue
			if indel.find('I')!= -1:
				continue
			
			perc = eval(perc)
			del_length = eval(indel.split('_')[0][1:None])
			positions = re.split('[LCIRD]', indel.split('_')[1])
			left_pos = eval(positions[1])
			right_pos = eval(positions[-1])
			if indel.find('C')!= -1:
				amb_len = eval(positions[2])
			else:
				amb_len = 0
			
			deletion_lengths = update_deletion_lengths(deletion_lengths, del_length, rep, perc)			
			right_positions = add_pos(right_positions, right_pos, perc)
			left_positions = add_pos(left_positions, left_pos, perc)
			deletion_length_and_pos = update_deletion_length_and_pos(deletion_length_and_pos, del_length, perc, left_pos, right_pos, amb_len)
			
	return deletion_lengths, left_positions, right_positions, deletion_length_and_pos

def update_deletion_length_and_pos(deletion_length_and_pos, del_length, perc, left_pos, right_pos, amb_len):
	if del_length in deletion_length_and_pos:
		if (left_pos, amb_len, right_pos) in deletion_length_and_pos[del_length]:
			deletion_length_and_pos[del_length][(left_pos,amb_len, right_pos)] += perc
		else:
			deletion_length_and_pos[del_length][(left_pos,amb_len, right_pos)] = perc
	else:
		deletion_length_and_pos[del_length] = {(left_pos,amb_len, right_pos):perc}
	return deletion_length_and_pos
	
def update_deletion_lengths(deletion_lengths, del_length, rep, perc):
	if del_length in deletion_lengths:
		if rep == '-':
			if rep in deletion_lengths[del_length]:
				deletion_lengths[del_length][rep] += perc
			else:
				deletion_lengths[del_length][rep] = perc
		else:
			len_rep = len(rep)
			if len_rep in deletion_lengths[del_length]:
				deletion_lengths[del_length][len_rep] += perc
			else:
				deletion_lengths[del_length][len_rep] = perc
	else:
		if rep == '-':
			deletion_lengths[del_length] = {rep:perc}
		else:
			deletion_lengths[del_length] = {len(rep):perc}
	return deletion_lengths
	
def add_pos(pos_dict, pos, value):
	if pos in pos_dict:
		pos_dict[pos] += value
	else:
		pos_dict[pos] = value
	return pos_dict
	
# def plot_positions(left_positions, right_positions):
	# left_positions_keys = range(min(left_positions.keys()), max(left_positions.keys())+1)
	# ind = np.arange(len(left_positions_keys))
	# left_position_counts = []
	# for left_pos in left_positions_keys:
		# if left_pos in left_positions:
			# left_position_counts.append(left_positions[left_pos])
		# else:
			# left_position_counts.append(0)
	# ax,fig = plt.subplots(figsize = (12,12))
	# plt.bar(ind,left_position_counts)
	# plt.xlabel('Left position')
	# plt.xticks(ind, left_positions_keys)
	# plt.savefig('C:/Users/jt20/Documents/Summer17/test_22_07/SRR3710048_3/Results_2/mapped_against_null/figures/deletions_left_pos_repeats_len_1.png')
	
	# right_positions_keys = range(min(right_positions.keys()), max(right_positions.keys())+1)
	# ind = np.arange(len(right_positions_keys))
	# right_position_counts = []
	# for right_pos in right_positions_keys:
		# if right_pos in right_positions:
			# right_position_counts.append(right_positions[right_pos])
		# else:
			# right_position_counts.append(0)
	# ax,fig = plt.subplots(figsize = (12,12))
	# plt.bar(ind,right_position_counts)
	# plt.xlabel('Right position')
	# plt.xticks(ind, right_positions_keys)
	# plt.savefig('C:/Users/jt20/Documents/Summer17/test_22_07/SRR3710048_3/Results_2/mapped_against_null/figures/deletions_right_pos_repeats_len_1.png')
	# plt.close()
	
def plot_deletion_lengths(deletion_lengths):
	deletion_length_keys = range(min(deletion_lengths.keys()), max(deletion_lengths.keys()) + 1)
	ind = np.arange(len(deletion_length_keys))
	total_perc = []
	for del_length in deletion_length_keys:
		ind_perc = 0
		if del_length in deletion_lengths:
			for rep in deletion_lengths[del_length]:
				ind_perc += deletion_lengths[del_length][rep]
		else:
			ind_perc = 1
		total_perc.append(ind_perc)
	
	height_stacks = [[],[],[],[],[],[],[],[],[]]
	for i in range(9):
		if i == 0:
			for j in range(len(deletion_length_keys)):
				del_length = deletion_length_keys[j]
				if del_length in deletion_lengths:
					if '-' in deletion_lengths[del_length]:
						height_stacks[i].append(deletion_lengths[del_length]['-']/total_perc[j])
				
					else:
						height_stacks[i].append(0)
				else:
					height_stacks[i].append(0)
		else:
			for j in range(len(deletion_length_keys)):
				del_length = deletion_length_keys[j]
				if del_length in deletion_lengths:
					if i in deletion_lengths[del_length]:
						height_stacks[i].append(deletion_lengths[del_length][i]/total_perc[j])
						
					else:
						height_stacks[i].append(0)
				else:
					height_stacks[i].append(0)
	# bottom_heights = []
	# for j in range(30):
		# bottom_heights.append([])
		# bottom_heights[j] = np.cumsum([height_stacks[i][j] for i in range(9)])
	# ax, fig = plt.subplots(figsize= (15,15))
	# plt.bar(ind, height_stacks[0], label = '-')

	# for i in range(1,9):
		# plt.bar(ind, height_stacks[i], bottom = [bottom_heights[x][i-1] for x in range(30)], label = str(i))
		# bottom_heights.append([height_stacks])
	# plt.xlabel('Length of deletion')
	# plt.xticks(ind, deletion_length_keys)
	# plt.legend()
	# plt.ylabel('Proportion of deletions')
	# plt.savefig('C:/Users/jt20/Documents/Summer17/test_22_07/SRR3710048_3/Results_2/mapped_against_null/figures/proportion_deletion_lengths_repeats.png')
	# plt.close()
	return deletion_length_keys, height_stacks

def scatter_plot_deletion_length_proportions(deletion_length_keys, height_stacks, sub_fig_dir_name):
	ind = np.arange(len(deletion_length_keys))

	for i in range(9):
		ax, fig = plt.subplots(figsize = (12,12))
		fit = np.polyfit(deletion_length_keys[0:27], height_stacks[i][0:27],deg = 1)
		plt.plot(ind,[fit[0]*deletion_length_key + fit[1] for deletion_length_key in deletion_length_keys])
		plt.scatter(deletion_length_keys,height_stacks[i], label = str(i))
		plt.axis((0,29,0,1))
		plt.xlabel('Length of deletion')
		plt.ylabel('Proportion of deletions')
		plt.title('Deletions related to length ' + str(i) +' repeats')
		plt.savefig(os.path.join(sub_fig_dir_name, 'mh_len_'+str(i)+'_del_len.png'))
		plt.close()

def plot_just_regression_lines(deletion_length_keys, height_stacks, fig_dir_name,permuted):
	ind = np.arange(len(deletion_length_keys))
	ax, fig = plt.subplots(figsize = (12,12))
	for i in range(9):
		fit = np.polyfit(deletion_length_keys[0:27], height_stacks[i][0:27],deg = 1)
		plt.plot(ind,[fit[0]*deletion_length_key + fit[1] for deletion_length_key in deletion_length_keys], label = str(i))
		#plt.scatter(deletion_length_keys,height_stacks[i], label = str(i))
	plt.legend()
	plt.axis((0,28,0,1))
	plt.xlabel('Length of deletion')
	plt.ylabel('Proportion of deletions')
	plt.title('Regression lines for proportion of deletions due to repeats')
	if permuted == 1:
		plt.savefig(os.path.join(fig_dir_name, 'mh_proportions_del_len_permuted.png'))
	else:
		plt.savefig(os.path.join(fig_dir_name, 'mh_proportions_del_len.png'))
	plt.close()

def make_dirs(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'For each deletion length it finds the proportion of deletions of this length which may be associated with each  length of repeat or microhomology in the target sequence and plots a scatter plot with this information as well as regression lines showing how the relationship varies with the length of the repeat.')
	parser.add_argument('base_dir', help = 'Base directory where output files are stored')
	parser.add_argument('--permuted', help = '1 if considering the relationships with permuted guides and indels, 0 (default) to just produce the normal plots.', type = int, default = 0)
	args = parser.parse_args()
	
	fig_dir_name = os.path.join(args.base_dir, 'Figures')
	
	if args.permuted == 1:
		indels_associated_condensed_file = os.path.join(args.base_dir, 'indels_associated_condensed_permuted.txt')
		sub_fig_dir_name = os.path.join(fig_dir_name, 'mh_lengths_regression_plots_permuted')
	else:
		sub_fig_dir_name = os.path.join(fig_dir_name, 'mh_lengths_regression_plots')
		indels_associated_condensed_file = os.path.join(args.base_dir, 'indels_associated_condensed.txt')
	
	make_dirs(fig_dir_name)
	make_dirs(sub_fig_dir_name)
	deletion_lengths, left_positions, right_positions, deletion_length_and_pos = parse_indels_associated_condensed(indels_associated_condensed_file)

	deletion_length_keys, height_stacks = plot_deletion_lengths(deletion_lengths)
	scatter_plot_deletion_length_proportions(deletion_length_keys, height_stacks, sub_fig_dir_name)
	plot_just_regression_lines(deletion_length_keys, height_stacks, fig_dir_name, args.permuted)
	# plot_positions(left_positions, right_positions)