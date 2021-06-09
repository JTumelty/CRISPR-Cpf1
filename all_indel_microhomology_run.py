import re
import numpy as np
import argparse
import os

def find_guide_list(guide_file):
	guide_list = []
	with open(guide_file,'r') as file_in:
		for line in file_in:
			if len(line.strip('\n')) > 0:
				guide_list.append(line.strip('\n').split('\t'))
	return guide_list

def find_microhomology_dict(guide_list):
	microhomology_dict = {}
	for guide in guide_list:
		sequence = 'TTTA' + guide[1] + 'AGCTTGGCGTAACTAGATCT'
		for length in range(1,10):
			for substring_start in range(0,47 - length):
				substring_to_match = sequence[substring_start:substring_start+length]
				for position in range(substring_start+1,47-length):
					if sequence[position:position+length] == substring_to_match:
						if (guide[0],guide[1], guide[2]) in microhomology_dict.keys():
							microhomology_dict[(guide[0],guide[1], guide[2])].append((substring_to_match, position-substring_start, substring_start, position))
						else:
							microhomology_dict[(guide[0],guide[1], guide[2])]= [(substring_to_match, position-substring_start, substring_start, position)]
	return microhomology_dict
	
def output_microhomologies(output_file, microhomology_dict):
	with open(output_file,'w') as file_out:
		for guide in microhomology_dict:
			for microhomology in microhomology_dict[guide]:
				file_out.write(str(guide[0]) + '\t' + guide[1] +'\t' + guide[2]+ '\t' + microhomology[0] + '\t' + str(microhomology[1]) +'\t' + str(microhomology[2]) + '\t' + str(microhomology[3])+ '\n')

def find_all_indels(summary_file, oligo_permutation):
	indels = {}
	with open(summary_file,'r') as f_in:
		for line in f_in:
			if line.strip('\n')[0] == '@':
				orig_no = eval(line.strip('\n').split('@@@Cpf1_Oligo')[1])
				oligo = oligo_permutation[orig_no -1]
				indels[oligo] = {}
			try:
				[indel, oligo_indel, count] = line.strip('\n').split('\t')
				if indel in indels[oligo]:
					indels[oligo][indel].append((oligo_indel,eval(count)))
				else:
					indels[oligo][indel] = [(oligo_indel, eval(count))]
			except ValueError:
				pass
	return indels

def indel_deletion_length(indels):
	indel_del_length = {}
	for i in range(1,1252):
		oligo = i
		indel_del_length[i] = {}

		for indel in indels[i].keys():
		
			if indel.find('I') != -1:
				continue
			if indel[0] == 'D':
				if indel.find('I') != -1:
					continue
				del_length = int(indel.split('_')[0][1:None])
				if del_length in indel_del_length[i].keys():
					indel_del_length[i][del_length].append(indel)
				else:
					indel_del_length[i][del_length] = [indel]
					
	return indel_del_length
	
def map_indel_microhomology_run(out_file, in_file, permuted_guides, indel_del_length):
	with open(out_file,'w') as file_out:
		with open(in_file,'r') as file_in:
			for line in file_in:
				found = False
				if len(line.strip('\n')) > 0:
					[guide_no,seq,guide_count, microhomology, length, start,stop] = line.strip('\n').split('\t')
					guide_no = int(guide_no)
					length = int(length)
					start = int(start)
				
					all_indels = []
					indels_right_length = []
					indels_to_write = []

					for k in indel_del_length[guide_no].keys():
						for indel in indel_del_length[guide_no][k]:
							all_indels.append(indel)

					for k in range(length,length+1):
						if k in indel_del_length[guide_no]:
							indels_right_length.extend(indel_del_length[guide_no][k])

					for indel in indels_right_length:
						del_length = int(re.split('[ILCRD]',indel.split('_')[1])[1])
						if del_length in range(start-22-1,start-22+2):
							indels_to_write.append(indel)
							found = True
							
					if found:
						file_out.write(line.strip('\n') + '\t' +','.join(indels_to_write) + '\t' + ','.join(all_indels) + '\n')
					else:
						file_out.write(line.strip('\n') + '\t' + '--' + '\t' + ','.join(all_indels) +'\n')

def run_analysing_microhomology_deletions(summary_indels,microhomologies_with_matches, indels_associated_file, permuted_guides):
	with open(microhomologies_with_matches,'r') as in_file:		
		guide_microhomologies_matched = {}
		for line in in_file:
			if len(line.strip('\n'))>0:
				line_info = line.strip('\n').split('\t')
				if line_info[7] != '--':
					my_indel = line_info[7].split(',')[0]
					if int(line_info[0]) in guide_microhomologies_matched:
						
						if my_indel in guide_microhomologies_matched[ int(line_info[0])]:
							guide_microhomologies_matched[int(line_info[0])][my_indel].append(line_info[3:7])
						else:
						
							guide_microhomologies_matched[int(line_info[0])][my_indel] =[line_info[3:7]]
					else:
						guide_microhomologies_matched[int(line_info[0])]= {my_indel:[line_info[3:7]]}
					

	with open(indels_associated_file,'w') as out_file:
		guide_indel = make_guide_indel_dict(summary_indels, permuted_guides)
		for guide in guide_indel.keys():
			
			for indel in guide_indel[guide]:
				if guide in guide_microhomologies_matched.keys():
					if indel[0] in guide_microhomologies_matched[guide]:
						for match in guide_microhomologies_matched[guide][indel[0]]:
							out_file.write(str(guide) + '\t' + str(indel[1]) + '\t' + indel[0] + '\t' + '\t'.join(match) + '\n')
					else:
						out_file.write(str(guide) + '\t' + str(indel[1]) + '\t' + indel[0]  + '\t' + '-' + '\n')
				else:
					out_file.write(str(guide) + '\t' + str(indel[1]) + '\t' + indel[0]  + '\t' + '-' + '\n')
					
def make_guide_indel_dict(summary_indels, permuted_guides):
	guide_indel = {}
	oligo_count = {}
	with open(summary_indels,'r') as file_in:
		for line in file_in:
			if line.strip('\n')[0] == '@':
				orig_no = eval(line.strip('\n').split('@@@Cpf1_Oligo')[1])
				oligo = permuted_guides[orig_no -1]
				oligo_count[oligo] = 0
				guide_indel[oligo] = []
				continue
			try:
				indel, oligo_indel, count = line.strip('\n').split('\t')
				guide_indel[oligo].append((indel, eval(count)))
				oligo_count[oligo] += eval(count)
			except ValueError:
				continue
	
	for oligo in guide_indel.keys():
		guide_indel[oligo] = [(x,y/float(oligo_count[oligo])) for (x,y) in guide_indel[oligo]]
	return guide_indel

def run_summarise_indels_associated(indels_associated_file, indels_associated_condensed_file):
	guide_dict = {}
	with open(indels_associated_file,'r') as in_file:
		for line in in_file:
			if len(line.strip('\n')) > 0:
				line_info = line.strip('\n').split('\t')
				guide = line_info[0]
				if guide not in guide_dict:
					guide_dict[guide] = {}
				indel = (line_info[1],line_info[2])
				if line_info[3] != '-':
					if indel in guide_dict[guide]:
						if guide_dict[guide][indel][1] != '-':
							if int(line_info[4]) >= int(guide_dict[guide][indel][1]):
								guide_dict[guide][indel] = line_info[3:]
						else:
							guide_dict[guide][indel] = line_info[3:]
					else:
						guide_dict[guide][indel] = line_info[3:]
				else:
					guide_dict[guide][indel] = ['-','-','-','-']

	with open(indels_associated_condensed_file,'w') as out_file:
		for i in range(1,1252):
			try:
				for indel in guide_dict[str(i)]:
					out_file.write(str(i) + '\t' +'\t'.join(indel) + '\t' + guide_dict[str(i)][indel][0] + '\t' +guide_dict[str(i)][indel][1] + '\t' + guide_dict[str(i)][indel][2] + '\t'+guide_dict[str(i)][indel][3] +  '\n')		
			except KeyError:
				import pprint
				pprint.pprint(guide_dict)
				exit()
def generate_random_guides_file(args, k):
	base_dir = args.base_dir
	out_file_name = os.path.join([base_dir, 'random_guides_'+str(k)+'.txt'])
	for i in range(1251):
		guide_char = []
		char = ['A','C','G','T']
		for j in range(23):
			rand_no = random.randint(0,3)
			guide_char.append(char[rand_no])
		out_file.write(str(i+1) + '\t' + ''.join(guide_char) + '\t' + '0' + '\n')

def run_normal(args):
	base_dir = args.base_dir
	all_microhomology_file = os.path.join(base_dir, 'all_microhomologies.txt')
	microhomologies_with_matches =  os.path.join(base_dir, 'all_microhomologies_matched.txt')
	indels_associated_file =   os.path.join(base_dir, 'indels_associated.txt') 
	summary_indels = args.summary_indels_file
	indels_associated_condensed_file =  os.path.join(base_dir, 'indels_associated_condensed.txt')
	guide_file = args.guides_file
	
	guide_list = find_guide_list(guide_file)
	microhomology_dict = find_microhomology_dict(guide_list)
	output_microhomologies(all_microhomology_file, microhomology_dict)
	indels = find_all_indels(summary_indels, range(1,1252))
	permuted_guides = range(1,1252)
	indel_del_length = indel_deletion_length(indels)
	map_indel_microhomology_run(microhomologies_with_matches,all_microhomology_file, permuted_guides, indel_del_length)
	run_analysing_microhomology_deletions(summary_indels, microhomologies_with_matches, indels_associated_file, permuted_guides)
	run_summarise_indels_associated(indels_associated_file, indels_associated_condensed_file)

def run_permuted(args):
	base_dir = args.base_dir
	all_microhomology_file = os.path.join(base_dir, 'all_microhomologies_permuted.txt')
	microhomologies_with_matches =  os.path.join(base_dir, 'all_microhomologies_matched_permuted.txt')
	indels_associated_file =   os.path.join(base_dir, 'indels_associated_permuted.txt')
	summary_indels = args.summary_indels_file
	indels_associated_condensed_file = os.path.join(base_dir, 'indels_associated_condensed_permuted.txt')
	guide_file = args.guides_file
	
	
	guide_list = find_guide_list(guide_file)
	microhomology_dict = find_microhomology_dict(guide_list)
	output_microhomologies(all_microhomology_file, microhomology_dict)
	permuted_guides = np.random.permutation(range(1,1252))
	indels = find_all_indels(summary_indels, permuted_guides)
	indel_del_length = indel_deletion_length(indels)
	map_indel_microhomology_run(microhomologies_with_matches,all_microhomology_file, permuted_guides, indel_del_length)
	run_analysing_microhomology_deletions(summary_indels, microhomologies_with_matches, indels_associated_file, permuted_guides)
	run_summarise_indels_associated(indels_associated_file, indels_associated_condensed_file)
		
def run_random_sequences(args):
	base_dir = args.base_dir
	all_microhomology_file = os.path.join(base_dir, 'all_microhomologies_random.txt')
	microhomologies_with_matches =  os.path.join(base_dir, 'all_microhomologies_matched_random..txt')
	indels_associated_file =   os.path.join(base_dir, 'indels_associated_random.txt')
	summary_indels = args.summary_indels_file
	indels_associated_condensed_file = os.path.join(base_dir, 'indels_associated_condensed_random.txt')
	guide_file = generate_random_guides_file(args, k)
	
	guide_list = find_guide_list(guide_file)
	microhomology_dict = find_microhomology_dict(guide_list)
	output_microhomologies(all_microhomology_file, microhomology_dict)
	indels = find_all_indels(summary_indels, range(1,1252))
	permuted_guides = range(1,1252)
	indel_del_length = indel_deletion_length(indels)
	map_indel_microhomology_run(microhomologies_with_matches,all_microhomology_file, permuted_guides, indel_del_length)
	run_analysing_microhomology_deletions(summary_indels, microhomologies_with_matches, indels_associated_file, permuted_guides)
	run_summarise_indels_associated(indels_associated_file, indels_associated_condensed_file)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Finds microhomologies of all lengths for each of the guide sequences (normal, random or permuted) and matches these with the deletions. ')
	parser.add_argument('base_dir', help = 'Directory for output files')
	parser.add_argument('summary_indels_file', help = 'All indels partitioned by guide')
	parser.add_argument('guides_file',help = 'Path to "all_guides.txt"')
	parser.add_argument('--normal_random_permuted', type=int,help = '0: Run with actual guides with matched indels, 1: Run with actual guides with indels swapped between guides, 2: Run with 1251 random guides with sets of indels associated with them. Default = 0', default = 0)
	args = parser.parse_args()
	
	if args.normal_random_permuted == 0:
		run_normal(args)
	elif args.normal_random_permuted == 1:
		run_permuted(args)
	elif args.normal_random_permuted == 2:
		run_random_sequences(args)
	
	
	
	
