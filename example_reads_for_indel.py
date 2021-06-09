import argparse
import os

def find_indels_mismatch(summary):
	oligo_indel_dict = {}
	with open(summary,'r') as in_file:
		for line in in_file:
			if len(line.strip()) == 0:
				continue
			if line.strip()[0] == '@':
				oligo = line.strip().split('@@@Cpf1_Oligo')[1]
				oligo_indel_dict[oligo] = {}
				continue
			indel = line.strip().split('\t')[0]
			if indel not in oligo_indel_dict[oligo]:
				oligo_indel_dict[oligo][indel] = {line.strip().split('\t')[1]:{}}
			else:	
				oligo_indel_dict[oligo][indel][line.strip().split('\t')[1]] = {}
	return oligo_indel_dict
	
def find_representative_reads(indelmap_output_dir, indelmap_output_name, reads_dir, oligo_indel_dict):
	for i in range(1,1252):
		oligo = str(i)
		reads_to_find = {}
		with open(os.path.join(indelmap_output_dir, indelmap_output_name + '_' + oligo + '.txt'),'r') as in_file:
			for line in in_file:
				if len(line.strip('\n')) == 0:
					continue
				read_id, oligo_id, indel, mismatch, time, count_1, count_2 = line.strip().split('\t')
				if oligo_id == '':
					continue
				indels = indel.split(',')
				oligo_ids = oligo_id.split(',')
				
				for j in range(len(indels)):
					indel = indels[j]
					oligo_id = oligo_ids[j]
					if len(oligo_id.split(':')[0].split('_')) == 3:
						oligo_indel = oligo_id.split(':')[0].split('_')[-1]+'_'+oligo_id.split(':')[1]
					elif len(oligo_id.split(':')[0].split('_')) == 4:
						oligo_indel = oligo_id.split(':')[0].split('_')[2]+'_'+oligo_id.split(':')[0].split('_')[3]+'_'+oligo_id.split(':')[1]

					if indel in oligo_indel_dict[oligo].keys():
						if oligo_indel in oligo_indel_dict[oligo][indel]:
							if len(oligo_indel_dict[oligo][indel][oligo_indel].keys()) == 0:
								oligo_indel_dict[oligo][indel][oligo_indel] = {read_id:None}
								if read_id in reads_to_find:
									reads_to_find[read_id].extend([indel, oligo_indel])
								else:
									reads_to_find[read_id] = [indel, oligo_indel]
		
		with open(os.path.join(reads_dir, 'Cpf1_Oligo'+oligo, 'all_reads.fastq'),'r') as in_file:
			for line in in_file:
				if len(line.strip()) == 0:
					continue
				if line.strip('\n')[1:None] in reads_to_find.keys():
					seq = in_file.next().strip()
					reads_to_find[line.strip('\n')[1:None]].append(seq)
					
		for read_id in reads_to_find.keys():
			no_indels = int(len(reads_to_find[read_id]))
			for j in range(int(no_indels/2)):
				indel = reads_to_find[read_id][2*j]
				oligo_indel = reads_to_find[read_id][2*j+1]
				seq = reads_to_find[read_id][-1]
				oligo_indel_dict[oligo][indel][oligo_indel][read_id] = seq
	return oligo_indel_dict
	
def output_example_reads(oligo_indel_dict,example_reads_indels):	
	with open(example_reads_indels,'w') as out_file:
		for i in range(1,1252):
			oligo = str(i)
			out_file.write('@@@Cpf1_Oligo'+str(i) + '\n')
			for indel in oligo_indel_dict[oligo]:
				for oligo_indel in oligo_indel_dict[oligo][indel]:
					for read in oligo_indel_dict[oligo][indel][oligo_indel]:
						out_file.write(indel + '\t' + oligo_indel + '\t' + oligo_indel_dict[oligo][indel][oligo_indel][read] + '\n')
					
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Finds an example of a read for each indel and oligo indel pair in the summary file. ')
	parser.add_argument('base_dir', help = 'Base directory where all output will be stored. Should already contain a folder for mapped files, summary files')
	parser.add_argument('indelmap_output_dir', help = 'Directory containing all of the mapped test files')
	parser.add_argument('indelmap_output_name', help = 'Example of base region in the name of the mapped files e.g. for "4_indelmap_output_32.txt", give "4_indelmap_output".')
	args = parser.parse_args()
	
	reads_dir = 'C:\\Users\\jt20\\Documents\\Summer17\\test_22_07\\SRR3710048_3\\Results_2\\Oligos'
	summary = os.path.join(args.base_dir, 'summary_file.txt')
	example_reads_indels = os.path.join(args.base_dir, 'example_reads_indels.txt')
	
	oligo_indel_dict = find_indels_mismatch(summary)
	oligo_indel_dict = find_representative_reads(args.indelmap_output_dir, args.indelmap_output_name, reads_dir, oligo_indel_dict)
	output_example_reads(oligo_indel_dict, example_reads_indels)