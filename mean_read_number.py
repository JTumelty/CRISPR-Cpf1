import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import os

def mean_read_count(args):
	counts = []
	with open(args.summary_indels_percent,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			if line.strip()[0] == '@':
				guide_no = eval(line.strip().split('@@@Cpf1_Oligo')[1])
				change = True
				continue
			indel, oligo_indel, count, perc = line.strip().split('\t')
			if change:
				counts.append(eval(count))
				change = False
			else:
				counts[-1]+= eval(count)
			mean_guide_count = np.mean(counts)
			read_std = np.std(counts)
	return mean_guide_count, read_std
	
def indels(args, mean_guide_count):
	indels_dict = {}
	with open(args.summary_indels_percent,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			if line.strip()[0] == '@':
				continue
			indel, oligo_indel, count, perc = line.strip().split('\t')
			if indel == '-':
				continue
			if indel in indels_dict:
				indels_dict[indel][0] += 1
				indels_dict[indel][1] += round(mean_guide_count*eval(perc))
			else:
				indels_dict[indel] = [1, round(mean_guide_count*eval(perc))]
	return indels_dict
	
def find_extreme_indels(indels_dict):
	for indel in indels_dict:
		if indels_dict[indel][0] >= 100 or indels_dict[indel][1] >= 1300:
			print indel, indels_dict[indel][0], indels_dict[indel][1]
			
def plot_indels(indel_count_and_guides_plot, indels_dict):
	indel_list = sorted(indels_dict.keys())
	indel_count = [indels_dict[indel][1] for indel in indel_list]
	indel_guide_count = [indels_dict[indel][0] for indel in indel_list]
	fig,ax = plt.subplots()
	plt.scatter(indel_count, indel_guide_count, s = 10)
	plt.xlabel('Number of indel occurrences')
	plt.ylabel('Number of guides')
	ax.set_yscale('log')
	plt.xlim((0,2000))
	ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	plt.savefig(indel_count_and_guides_plot)
	plt.close()
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description= 'Computes the mean number of reads mapped across all guides. Creates a scatter plot showing, for each indel, the number of times that it occurs and the number of guides it was found in.')
	parser.add_argument('base_dir', help = 'Path to directory to store output files')
	parser.add_argument('-in','--summary_indels_percent', help = 'Path to summary indels percent file', default = 'C:\\Users\\jt20\\Documents\\Summer17\\Final_results\\summary_file_percent.txt')
	args = parser.parse_args()

	indel_count_and_guides_plot = os.path.join(args.base_dir, 'Figures', 'indels_count_vs_guide.png')
	mean_guide_count, read_std = mean_read_count(args)
	print mean_guide_count
	indels_dict = indels(args, mean_guide_count)
	find_extreme_indels(indels_dict)
	plot_indels(indel_count_and_guides_plot, indels_dict)