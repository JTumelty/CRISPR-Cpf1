import io
import argparse
import os

def insert_repeat(target, indel, final):
	indel_length = eval(indel.split('_')[0][1:None])
	right_det = eval(indel[indel.find('R')+1:None])
	if indel.find('C')!= -1:
		left_det = eval(indel[indel.find('L')+1:indel.find('C')])
	else:
		left_det = eval(indel[indel.find('L')+1:indel.find('R')])
	last_temp = target[18+right_det:None]
	beginning = target[0: 19+left_det]
	initial = target[19+left_det:18+right_det]
	insertion = final[23+left_det: len(final)- len(last_temp)]
	
	if indel_length == 1 and indel.find('C') == -1:
		for i in range(len(target)):
			if final[4:4+i] == target[0:i] and final[4+i+1:None] == target[i:None]:
				insertion = final[4+i]
				break
	return initial, insertion

def find_guide_dict(all_guides):	
	guide_dict = {}
	with open(all_guides,'r') as in_file:
		for line in in_file:
			try:
				guide_no,guide,count = line.strip('\n').split('\t')
				guide_dict[eval(guide_no)] = guide
			except ValueError:
				continue
	return guide_dict

def find_insertions(indels_associated_condensed):
	insertions = {}
	with open(indels_associated_condensed, 'r') as in_file:
		for line in in_file:
			if len(line.strip('\n')) == 0:
				continue
			guide, perc, indel, repeat, length, start, stop = line.strip().split('\t')
			if indel[0] == 'I' and indel.find('D') == -1:
				if eval(guide) in insertions:
					if indel not in insertions[eval(guide)]:
						insertions[eval(guide)][indel] = eval(perc)
					else:
						insertions[eval(guide)][indel] += eval(perc)
				else:
					insertions[eval(guide)] = {indel:eval(perc)}
	return insertions


def find_all_insertions(insertions, example_reads_indels):
	all_insertions = []
	for guide in insertions:
		target = guide_dict[guide] + 'AGCTTGGCGTAACTAGATCT'
		for insertion in insertions[guide]:
			if insertion.find('C')!= -1:
				left_det = eval(insertion.split('_')[1][1:min(insertion.split('_')[1].find('C'), insertion.split('_')[1].find('R'))])
				right_det = eval(insertion.split('R')[1])
				amb_det = eval(insertion.split('_')[1][insertion.split('_')[1].find('C')+1:insertion.split('_')[1].find('R')])
			else:
				left_det = eval(insertion.split('_')[1][1:insertion.split('_')[1].find('R')])
				right_det = eval(insertion.split('R')[1])
				amb_det = None
			insertion_length = eval(insertion.split('_')[0][1:None])
			with open(example_reads_indels,'r') as in_file:
				for line in in_file:
					if line.strip('\n')[0] == '@':
						oligo = eval(line.strip().split('@@@Cpf1_Oligo')[1])
						continue
					if oligo == guide:
						if line.startswith(insertion):
							template_init = target[19+left_det:18+right_det]
							seq = line.strip().split('\t')[2]
							guide_start = seq.find('TTTTT')
			
							if seq[23+guide_start:27+guide_start] == 'TTTA':
								inserted= seq[guide_start+23:None]
							else:
								for i in [22,24,21,25, 20, 26, 19, 27, 18, 28, 17, 29]:
									if seq[i+guide_start:i+4+guide_start] == 'TTTA':
										inserted = seq[guide_start+i:None]
										break
								else:
										continue
							initial, inserted_2 = insert_repeat(target,insertion, inserted)
							all_insertions.append((str(guide), insertion, str(insertions[guide][insertion]), inserted))
							break

	return all_insertions
				
def output_insertions(all_insertions_file, all_insertions, guide_dict):
	with open(all_insertions_file,'w') as out_file:
		for insertion in all_insertions:
			target = guide_dict[eval(insertion[0])]+'AGCTTGGCGTAACTAGATCT'
			if insertion[-1][4:] == target:
				continue
			out_file.write('\t'.join(insertion) + '\t' + target + '\n')
	

if __name__ == '__main__':
		parser = argparse.ArgumentParser(description = 'Extracts all the insertions and examples of relevant reads and outputs them into a file called "all_insertions.txt".')
		parser.add_argument('base_dir', help = 'Base directory where output will be stored. This directory should contain "all_guides.txt" in a subdirectory "core_files", "indels_associated_condensed.txt" and "example_reads_indels.txt"')
		args = parser.parse_args()
		
		guide_dict = find_guide_dict(os.path.join(args.base_dir, 'core_files', 'all_guides.txt'))
		insertions = find_insertions(os.path.join(args.base_dir, 'indels_associated_condensed.txt'))
		all_insertions = find_all_insertions(insertions, os.path.join(args.base_dir, 'example_reads_indels.txt'))
		output_insertions(os.path.join(args.base_dir, 'all_insertions.txt'), all_insertions, guide_dict)
