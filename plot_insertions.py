import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

def condense_insertions(all_insertions_with_insert, out_file):
	with open(all_insertions_with_insert,'r') as in_file:
		with open(out_file,'w') as out_file:
			for line in in_file:
				line_info = line.strip().split('\t')
				insertions = line_info[5]
				end_7_plus_A = line_info[4][16:24]
				out_file.write(line_info[1] + '\t' +insertions + '\t' + end_7_plus_A + '\n')

def classify_insertions(out_file):
	output_dict = {}
	with open(out_file, 'r') as in_file:
		for line in in_file:
			line_info = line.strip().split('\t')
			ins_len = eval(line_info[0].split('_')[0][1:])
			if ins_len >= 1:
				if ins_len not in output_dict:
					output_dict[ins_len] = [1,{}]
				else:
					output_dict[ins_len][0] += 1
				pos_positions = [eval(line_info[1].split(':')[i].split('_')[1])  for i in range(len(line_info[1].split(':')))]
				for (pos_insert, position) in [(line_info[1].split(':')[i].split('_')[0],eval(line_info[1].split(':')[i].split('_')[1])) for i in range(len(line_info[1].split(':')))]:
					if line_info[2][position-17:].find(pos_insert) != -1:
						pos_in_guide = line_info[2][position-17:].find(pos_insert) +position
						if pos_in_guide in pos_positions:						
							if pos_in_guide not in output_dict[ins_len][1]:
								output_dict[ins_len][1][pos_in_guide] = [1, [(line_info[1],line_info[2])]]
							else:	
								output_dict[ins_len][1][pos_in_guide][0] += 1
								output_dict[ins_len][1][pos_in_guide][1].append((line_info[1],line_info[2]))
							break
				else:
					if '-' not in output_dict[ins_len][1]:
						output_dict[ins_len][1]['-'] = [1, [(line_info[1],line_info[2])]]
					else:
						output_dict[ins_len][1]['-'][0] += 1
						output_dict[ins_len][1]['-'][1].append((line_info[1],line_info[2]))

	return output_dict
			
def plot_positions(output_dict, fig_dir_name):
	found_pos = range(5,25)
	fig, ax = plt.subplots()
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
	plt.ylabel('Percentage of insertions ')
	plt.ylim((0,50))
	plt.xlim((5,25))
	plt.xticks(range(5,25))
	plt.savefig(os.path.join(fig_dir_name, 'insertions_match_guide_positions.png'))
	plt.close()
	
	fig,ax = plt.subplots()
	height_counts = []
	for ins_len in range(1,6):
		height_counts.append(0)
		for pos in found_pos:
			if pos in output_dict[ins_len][1]:
				height_counts[ins_len-1]+= output_dict[ins_len][1][pos][0]
				
	unmatched_count = []
	for ins_len in range(1,6):
		unmatched_count.append(output_dict[ins_len][1]['-'][0])
	plt.bar(range(1,6), height_counts, label = 'Insertion found in guide')
	plt.bar(range(1,6), unmatched_count, bottom = height_counts, label = 'Insertion not found in guide')
	plt.ylabel('Number of insertions')
	plt.xlabel('Insertion size')
	plt.legend()
	plt.savefig(os.path.join(fig_dir_name, 'insertions_match_guide_not.png'))
	plt.close()

def make_dir(fig_dir_name):
	if not os.path.exists(fig_dir_name):
		os.makedir(fig_dir_name)

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
	output_dict = classify_insertions(out_file)
	plot_positions(output_dict, fig_dir_name)