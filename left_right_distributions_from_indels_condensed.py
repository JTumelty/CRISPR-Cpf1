import matplotlib.pyplot as plt
import numpy as np
import mean_read_number
import indel_frequency_vs_read_count
import os
import argparse

""" One-sentence description
@param indels_associated_condensed_file string path to input file
@requires guides is not Null
@requires indels_associated_condensed_file has 8 columns that are ..
@effects creates a new file with the condensed results
@return map of (left,right) coordinates to number of indels observed
"""
def left_v_right_right_v_left(indels_associated_condensed_file, mean_reads, guides=range(1,1252)):
	right_left_dict = {}
	left_right_dict = {}
	with open(indels_associated_condensed_file,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			indel = line.strip().split('\t')[2]
			if indel.find('I') != -1:
				continue
			if indel == '-':
				continue
			rep = line.strip().split('\t')[3]
			if rep != '-':
				continue
			if indel.find('C') != -1:
				continue

			if eval(line.strip().split('\t')[0]) in guides:
				perc = eval(line.strip().split('\t')[1])
		
				right_pos = eval(indel[indel.find('R')+1:None])
				left_pos = eval(indel[indel.find('L')+1:indel.find('R')])
				
				right_left_dict = add_to_dict(right_left_dict, right_pos, left_pos, round(perc*mean_reads))
				left_right_dict = add_to_dict(left_right_dict, left_pos, right_pos, round(perc*mean_reads))

	return right_left_dict, left_right_dict

def add_to_dict(my_dict, entry1, entry2,value):
	if entry1 in my_dict:
		if entry2 in my_dict[entry1]:
			my_dict[entry1][entry2] += value
		else:
			my_dict[entry1][entry2] = value
	else:
		my_dict[entry1] = {entry2:value}
	return my_dict

def plot_heatmap(right_left_dict, left_right_dict, file_name):
	fig,ax = plt.subplots(figsize = (15,15))
	left_pos = list(reversed(range(-19,5)))
	right_pos = list(reversed(range(-6,15)))
	left_right = np.zeros((21,24))
	for i in range(24):
		for j in range(21):
			left = left_pos[23-i]
			right = right_pos[j]
			if left in left_right_dict:
				if right in left_right_dict[left]:
					left_right[j][i] = left_right_dict[left][right]
	
	plt.imshow(left_right, cmap= 'hot_r', extent = (-19,5,-6,17))
	plt.xticks([x+ 0.5 for x in range(-19,5)], [str(x+19) for x in range(-19,5)])
	plt.yticks([x+ 0.5 for x in range(-6,15)], [str(x+19) for x in range(-6,15)])
	plt.xlabel('Left (number of nucleotides after PAM)')
	plt.ylabel('Right (number of nucleotides after PAM)')
	plt.title('Position of unambiguous deletions')
	cbar = plt.colorbar(fraction=0.046, pad=0.04)
	cbar.set_label('Number of deletions across all guides')
	plt.savefig(file_name)
	plt.close()

def make_dirs(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)


def main():
	parser = argparse.ArgumentParser(description = 'Produces a heat map showing the left and right positions of each deletion, when there is no ambiguity or microhomology associated. Also produces heat maps for each of the separate activity groups as well as a cumulative option.')
	parser.add_argument('base_dir', help = 'Path to base directory where the heat maps will be stored')
	parser.add_argument('summary_indels_percent', help = 'Path to the summary indels file with percent')
	parser.add_argument('--all_guides', help = '0 to produce heat map with the distribution for all guides, 1 to not. Default = 1', type = int, default = 1)
	parser.add_argument('--separate_activity_groups', help = 'Produce heat map with the distribution for inputted number of separate activity groups, 0 to not. Default = 0', type = int, default = 0)
	parser.add_argument('--cumulative_activity_groups', help = 'Produce heat map for the cumulative distribution with increasing guide activity, 0 to not. Default = 10', type = int, default = 0)
	parser.add_argument('--filter', help = 'Lower limit on read number to consider the guide in the plots. Default = 100', type = int, default = 100)
	args = parser.parse_args()
	
	indels_associated_condensed_file = os.path.join(args.base_dir, 'indels_associated_condensed.txt')
	
	mean_reads, std = mean_read_number.mean_read_count(args)

	
	if args.all_guides == 0:
		right_left_dict, left_right_dict = left_v_right_right_v_left(indels_associated_condensed_file, mean_reads)
		plot_heatmap(right_left_dict, left_right_dict, os.path.join(args.base_dir,'Figures', 'heat_left_v_right_all_guides.png'))
		
	if args.cumulative_activity_groups != 0:
		cumulative_left_right_dir_name = os.path.join(args.base_dir, 'Figures','Cumulative_activity_groups')
		make_dirs(cumulative_left_right_dir_name)
		count = 0
		guides_considered = []
		activity_groups = indel_frequency_vs_read_count.split_by_activity(args,filter= args.filter, no_groups = args.cumulative_activity_groups )
		for activity_group in activity_groups:
			guides_considered.extend([guide[0] for guide in activity_group])
			right_left_dict, left_right_dict = left_v_right_right_v_left(indels_associated_condensed_file, mean_reads, guides = guides_considered)
			plot_heatmap(right_left_dict, left_right_dict, os.path.join(cumulative_left_right_dir_name,''.join(['heat_left_v_right_group_cumulative_',str(count+1),'_filtered_', str(args.filter), '.png'])))
			count += 1
	
	if args.separate_activity_groups != 0:
		separate_left_right_dir_name = os.path.join(args.base_dir, 'Figures','Separate_activity_groups')
		make_dirs(separate_left_right_dir_name)
		count = 0
		activity_groups = indel_frequency_vs_read_count.split_by_activity(args,filter= args.filter, no_groups = args.separate_activity_groups)
		for activity_group in activity_groups:
			right_left_dict, left_right_dict = left_v_right_right_v_left(indels_associated_condensed_file, mean_reads, guides = activity_group)
			plot_heatmap(right_left_dict, left_right_dict, os.path.join(separate_left_right_dir_name, ''.join(['heat_left_v_right_group_separate_',str(count+1),'_filtered_', str(args.filter), '.png'])))
			count += 1	
	
			
if __name__ == '__main__':
	main()
