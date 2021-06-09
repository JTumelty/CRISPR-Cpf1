import argparse
import os

def find_indel_dict(args):
	indel_dict = {}
	for i in range(1,1252):
		indel_dict[i] = {}
		total_count = 0
		in_file_name = os.path.join(args.indelmap_output_directory,'indelmap_output_'+str(i)+'.txt')
		with open(in_file_name,'r') as in_file:
			for line in in_file:
				if len(line.strip()) == 0:
					continue
				toks = line.split('\t')
				if toks[1] == '':
					continue
				if toks[2] in indel_dict[i]:
					if toks[3] in indel_dict[i][toks[2]]:
						indel_dict[i][toks[2]][toks[3]] += 1
					else:
						indel_dict[i][toks[2]][toks[3]] = 1
				else:
					indel_dict[i][toks[2]] = {toks[3]:1}
	return indel_dict

def output_indels_in_summary_file(summary_file):
	with open(summary_file,'w') as out_file:
		for i in range(1,1252):
			out_file.write('@@@Cpf1_Oligo'+str(i)+ '\n')
			for indel in indel_dict[i]:
				for mismatch in indel_dict[i][indel]:
					out_file.write(indel + '\t' + mismatch + '\t' + str(indel_dict[i][indel][mismatch]) +'\n')
def make_results_directory(base_dir):
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
		
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('base_dir', help= 'Base directory for output to be stored')
	parser.add_argument('indelmap_output_directory', help = 'Path to the driectory where the mapped files are stored')
	args = parser.parse_args()
	
	make_results_directory(args.base_dir)
	indel_dict = find_indel_dict(args)
	summary_file = os.path.join(args.base_dir, 'summary_file.txt')
	output_indels_in_summary_file(summary_file)