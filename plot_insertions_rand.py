import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import pprint

def condense_insertions(all_insertions_with_insert, out_file):
	with open(all_insertions_with_insert,'r') as in_file:
		with open(out_file,'w') as out_file:
			for line in in_file:
				line_info = line.strip().split('\t')
				insertions = line_info[5]
				end_7_plus_A = line_info[4][16:24]
				out_file.write(line_info[1] + '\t' + line_info[2]+'\t' +insertions + '\t' + end_7_plus_A + '\t' + line_info[-2]+'\n')

def classify_insertions(out_file): 
	output_dict = {}
	pos_left_right_1 = {}
	pos_left_right_2 = {}
	with open(out_file, 'r') as in_file:
		for line in in_file:
			line_info = line.strip().split('\t')
			ins_len = eval(line_info[0].split('_')[0][1:])
			if ins_len >= 1:
				if ins_len not in output_dict:
					output_dict[ins_len] = [eval(line_info[1]),{}]
				else:
					output_dict[ins_len][0] += eval(line_info[1])
					
				pos_positions = [eval(line_info[2].split(':')[i].split('_')[1])  for i in range(len(line_info[2].split(':')))]
				tmp_info = []
				for (pos_insert, position) in [(line_info[2].split(':')[i].split('_')[0],eval(line_info[2].split(':')[i].split('_')[1])) for i in range(len(line_info[2].split(':')))]:
					if line_info[-1][position-1:position+ins_len-1] == pos_insert:
						tmp_info.append((position, line_info[2],line_info[3]))
			
				if len(tmp_info)>= 1:
					choice = np.random.randint(0, len(tmp_info))
					pos_in_guide = tmp_info[choice][0]
					if pos_in_guide not in output_dict[ins_len][1]:
						output_dict[ins_len][1][pos_in_guide] = [eval(line_info[1]), [(tmp_info[choice][1],tmp_info[choice][2])]]
					else:	
						output_dict[ins_len][1][pos_in_guide][0] += eval(line_info[1])
						output_dict[ins_len][1][pos_in_guide][1].append((tmp_info[choice][1],tmp_info[choice][2]))
						
					if pos_in_guide not in pos_left_right_1:
						pos_left_right_1[pos_in_guide] = {}
					if pos_in_guide + ins_len -1 not in pos_left_right_1[pos_in_guide]:
						pos_left_right_1[pos_in_guide][pos_in_guide + ins_len-1] = 1
					else:
						pos_left_right_1[pos_in_guide][pos_in_guide + ins_len-1] += 1
						
					if pos_in_guide +ins_len not in pos_left_right_2:
						pos_left_right_2[pos_in_guide+ ins_len] = {}
					if pos_in_guide + 2*ins_len -1 not in pos_left_right_2[pos_in_guide+ins_len]:
						pos_left_right_2[pos_in_guide+ins_len][pos_in_guide + 2*ins_len-1] = 1
					else:
						pos_left_right_2[pos_in_guide+ins_len][pos_in_guide + 2*ins_len-1] += 1
					
				else:
					if '-' not in output_dict[ins_len][1]:
						output_dict[ins_len][1]['-'] = [eval(line_info[1]), [(line_info[2],line_info[3])]]
					else:
						output_dict[ins_len][1]['-'][0] += eval(line_info[1])
						output_dict[ins_len][1]['-'][1].append((line_info[2],line_info[3]))
						
	return output_dict, pos_left_right_1, pos_left_right_2
			
def plot_positions(output_dict, fig_dir_name):
	found_pos = range(5,25)
	fig, ax = plt.subplots(figsize = (12,12))
	for ins_len in range(1,6):
		counts = []
		for position in found_pos:
			if position in output_dict[ins_len][1]:
				counts.append(100*output_dict[ins_len][1][position][0]/float(output_dict[ins_len][0]))
			else:
				counts.append(0)

		plt.plot(found_pos,counts, label = str(ins_len))
	plt.legend(title = 'Insertion size')
	plt.xlabel('Position of inserted bases in guide sequence')
	plt.ylabel('Percentage of insertions that length')
	plt.ylim((0,60))
	plt.xlim((5,25))
	plt.xticks(range(5,25))
	plt.savefig(os.path.join(fig_dir_name, 'insertions_match_guide_positions.png'))
	plt.close()
	
	fig, ax = plt.subplots(figsize = (12,12))
	for ins_len in range(6,10):
		counts = []
		for position in found_pos:
			if position in output_dict[ins_len][1]:
				counts.append(100*output_dict[ins_len][1][position][0]/float(output_dict[ins_len][0]))
			else:
				counts.append(0)

		plt.plot(found_pos,counts, label = str(ins_len))
	plt.legend(title = 'Insertion size')
	plt.xlabel('Position of inserted bases in guide sequence')
	plt.ylabel('Percentage of insertions that length')
	plt.ylim((0,100))
	plt.xlim((5,25))
	plt.xticks(range(5,25))
	plt.savefig(os.path.join(fig_dir_name, 'insertions_match_guide_positions_long_ins.png'))
	plt.close()
	
	fig,ax = plt.subplots(figsize = (12,12))
	height_counts = []
	for ins_len in range(1,12):
		height_counts.append(0)
		for pos in found_pos:
			if pos in output_dict[ins_len][1]:
				height_counts[ins_len-1]+= output_dict[ins_len][1][pos][0]
				
	unmatched_count = []
	for ins_len in range(1,12):
		unmatched_count.append(output_dict[ins_len][1]['-'][0])
	plt.bar(range(1,12), height_counts, label = 'Inserted sequence in target at position of insertion')
	plt.bar(range(1,12), unmatched_count, bottom = height_counts, label = 'Inserted sequence not in target at position of insertion')
	plt.xlabel('Insertion size')
	plt.xticks(range(1,12))
	plt.legend()
	plt.savefig(os.path.join(fig_dir_name, 'insertions_match_guide_not.png'))
	plt.close()
	

def make_dir(fig_dir_name):
	if not os.path.exists(fig_dir_name):
		os.makedir(fig_dir_name)

def plot_heatmap_insertions(pos_left_right_1, pos_left_right_2, fig_dir_name):
	fig,ax = plt.subplots(figsize = (12,12))
	array_1 = np.zeros([20,20])
	for left_pos in range(20):
		for right_pos in range(20):
			if left_pos+8 in pos_left_right_1:
				if right_pos+8 in pos_left_right_1[left_pos+8]:
					array_1[left_pos][right_pos] = pos_left_right_1[left_pos+8][right_pos+8]
	plt.imshow(array_1, cmap = 'hot_r')
	cbar = plt.colorbar(fraction=0.046, pad=0.04)
	plt.xticks([x for x in range(19)], range(8,27))
	plt.xlabel('End of insertion')
	plt.yticks([x for x in range(19)],range(8,27))
	plt.ylabel('Start of insertion')
	plt.savefig(os.path.join(fig_dir_name, 'insert_pos_1.png'))
	pprint.pprint(pos_left_right_1)
	
	fig, ax = plt.subplots(figsize = (12,12))
	array_1 = np.zeros([28,28])
	for left_pos in range(9,36):
		for right_pos in range(9,36):
			if left_pos in pos_left_right_2:
				if right_pos in pos_left_right_2[left_pos]:
					array_1[left_pos-9][right_pos-9] = pos_left_right_2[left_pos][right_pos]
	plt.imshow(array_1, cmap = 'hot_r')
	cbar = plt.colorbar(fraction=0.046, pad=0.04)
	plt.xticks([x for x in range(27)], range(9,36))
	plt.xlabel('End of insertion')
	plt.yticks([x for x in range(27)],range(9,36))
	plt.ylabel('Start of insertion')
	pprint.pprint(pos_left_right_2)
	plt.savefig(os.path.join(fig_dir_name, 'insert_pos_2.png'))
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Looks at the insertions and determines if the inserted region is within the target sequence, if it is then the position is plotted for each length of insertion.' \
													+ 'A second plot showing the count of whether it is or is not in the guide is also produced.')
	parser.add_argument('base_dir', help = 'Base directory for all output files. Requires "all_insertions_with_insert.txt" in directory.')
	args = parser.parse_args()
	
	fig_dir_name = os.path.join(args.base_dir, 'Figures')
	make_dir(fig_dir_name)
	
	all_insertions_with_insert = os.path.join(args.base_dir, 'all_insertions_with_insert.txt')
	out_file = os.path.join(args.base_dir, 'all_insertions_with_insert_condensed.txt')
	
	condense_insertions(all_insertions_with_insert,out_file)
	output_dict, pos_left_right_1, pos_left_right_2 = classify_insertions(out_file)
	plot_positions(output_dict, fig_dir_name)
	plot_heatmap_insertions(pos_left_right_1, pos_left_right_2, fig_dir_name)