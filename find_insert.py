import os
import argparse

def parse_all_insertions(all_insertions, all_insertions_with_insert, all_insertions_with_insert_unmatched, all_insertions_with_insert_wrong_len):
	with open(all_insertions,'r') as in_file:
		with open(all_insertions_with_insert,'w') as out_file:
			with open(all_insertions_with_insert_unmatched,'w') as out_file_2:
				with open(all_insertions_with_insert_wrong_len,'w') as out_file_3:
					for line in in_file:
						if len(line.strip()) == 0:
							continue
						guide_no, indel, perc, end, target = line.strip().split('\t')
						end = end[4:]
		
						ins_len = eval(indel.split('_')[0][1:])
						inserted = []
						score = []
						if len(end) == len(target) + ins_len:
							for i in range(len(target)):
								ind_score = 0
								for j in range(i):
									if end[j] != target[j]:
										ind_score += 1
								for j in range(i, len(target)):
									if end[j+ins_len] != target[j]:
										ind_score += 1
								score.append(ind_score)
							min_score = min(score)
							for i in range(len(target)):
								if score[i] == min_score:
									inserted.append( end[i:i+ins_len] + '_'+str(i+1))
							out_file.write(line.strip() + '\t' + ':'.join(inserted)+ '\n')
						else:
							out_file_3.write(line.strip() + '\n')
						
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Finds potential insertions')
	parser.add_argument('base_dir', help = 'Base directory to store output files. Should also containt the file "all_insertions.txt".')
	args = parser.parse_args()
	
	all_insertions = os.path.join(args.base_dir, 'all_insertions.txt')
	all_insertions_with_insert = os.path.join(args.base_dir, 'all_insertions_with_insert.txt')
	all_insertions_with_insert_unmatched = os.path.join(args.base_dir, 'all_insertions_with_insert_unmatched.txt')
	all_insertions_with_insert_wrong_len = os.path.join(args.base_dir, 'all_insertions_with_insert_wrong_len.txt')
	parse_all_insertions(all_insertions,all_insertions_with_insert, all_insertions_with_insert_unmatched, all_insertions_with_insert_wrong_len)