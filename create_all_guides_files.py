import argparse
import os

def parse_nmeth(nmeth):
	guide_dict = {}
	guide_no = 1
	with open(nmeth, 'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			guide, all, PAM, target, perc = line.strip().split('\t')
			guide_dict[guide_no] = [guide, perc, 0]
			guide_no += 1
	return guide_dict

def parse_summary(summary, guide_dict):
	with open(summary, 'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			if line.strip()[0] == '@':
				guide_no = eval(line.strip().split('@@@Cpf1_Oligo')[1])
				continue
			indel, oligo_indel, count = line.strip().split('\t')
			guide_dict[guide_no][2] += eval(count)
	return guide_dict
	
def produce_all_guides_freq_file(all_guides_freq, guide_dict):
	with open(all_guides_freq, 'w') as out_file:
		for guide_no in range(1,1252):
			out_file.write(str(guide_no) + '\t' + guide_dict[guide_no][0] + '\t' + str(guide_dict[guide_no][2]) + '\t' + guide_dict[guide_no][1] + '\n')
	
def produce_all_guides_file(all_guides, guide_dict):
	with open(all_guides, 'w') as out_file:
		for guide_no in range(1,1252):
			out_file.write(str(guide_no) + '\t' + guide_dict[guide_no][0] + '\t' + str(guide_dict[guide_no][2]) +'\n')
			
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Produces the files "all_guides.txt" and "all_guides_freq.txt" from "nmeth.4104-f3.txt" and "summary_file.txt" and outputs them in the core files folder in the base directory.')
	parser.add_argument('base_dir', help = 'Base directory for output files to be stored.')
	args = parser.parse_args()
	
	core_dir = os.path.join(args.base_dir, 'core_files')
	nmeth = os.path.join(core_dir, 'nmeth.4104-f3.txt')
	summary = os.path.join(args.base_dir, 'summary_file.txt')
	all_guides = os.path.join(core_dir, 'all_guides.txt')
	all_guides_freq = os.path.join(core_dir, 'all_guides_freq.txt')
	
	guide_dict = parse_nmeth(nmeth)
	guide_dict = parse_summary(summary, guide_dict)
	produce_all_guides_freq_file(all_guides_freq, guide_dict)
	produce_all_guides_file(all_guides, guide_dict)