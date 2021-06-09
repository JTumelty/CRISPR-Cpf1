import argparse
import re
import numpy as np
import os
import pprint

def add_read(indel, indel_dict, oligo, oligo_null):
	if indel in indel_dict[oligo]:
		if oligo_null in indel_dict[oligo][indel]:
			indel_dict[oligo][indel][oligo_null] += 1
		else:
			indel_dict[oligo][indel][oligo_null] = 1
	else:
		indel_dict[oligo][indel] = {oligo_null: 1}
	return indel_dict
	
def parse_indelmap_results(args):
	indel_dict = {}
	for i in range(1, 1252):
		oligo = 'Cpf1_Oligo' + str(i)
		indel_dict[oligo] = {}
		indelmap_output = os.path.join(args.indelmap_output, '4_indelmap_output_'+str(i)+'.txt')
		with open(indelmap_output, 'r') as in_file:
			for line in in_file:
					if len(line.strip('\n'))==0:
						continue
					toks = line.strip('\n').split('\t')
					
					if toks[1] == 'None' or toks[1] == '':
						continue
					try:
						nulls = [(indel, null.split('Oligo_'+str(i)+'_')[1].split(':')[0],null.split(':')[1],eval(null.split(':')[2]), muts) for (null,indel,muts) in zip(toks[1].split(','), toks[2].split(','),toks[3].split(','))]
					except IndexError:
						print line
						exit()
					nulls = [(x[0], x[1]+'_'+x[2], x[3],x[4]) for x in nulls]

					if len(nulls) == 1:
						indel, oligo_null, perc_oligo_indel,muts = nulls[0]
						indel_dict = add_read(indel, indel_dict, oligo, oligo_null)

					else:
						#If there is a null oligo that results in a null indel, select that
						for (indel, oligo_null, perc_oligo_indel, muts) in nulls:

							if indel == '-':
								if oligo_null == '-_-':
									indel_dict = add_read(indel, indel_dict, oligo, oligo_null)
									break
						if indel == '-' and oligo_null == '-_-':
							continue
						
						#Otherwise, if 18 closest region in gRNA or PAM of oligo is affected, exclude from consideration,
						#otherwise assign randomly according to oligo percentages in plasmids
						
						ok_nulls = [x for x in nulls if indelOutofGuideSeedPAM(x[1])]
						if len(ok_nulls) == 0:
							ok_nulls = nulls
						percs = [x[2] for x in ok_nulls]
						total_perc = sum(percs)
						idx = np.where(np.random.multinomial(1,np.array(percs)/total_perc)>0)[0][0]
						selected_null = ok_nulls[idx]
						
				
						indel = selected_null[0]
						oligo_null = selected_null[1]
						muts = selected_null[2]
						indel_dict = add_read(indel, indel_dict, oligo, oligo_null)

	return indel_dict

	
def tokFullIndel(indel):
	indel_toks = indel.split('_')
	indel_type, indel_details = indel_toks[0], ''
	if len(indel_toks) > 1:
		indel_details =  indel_toks[1]
	cigar_toks = re.findall(r'([CLRDI]+)(-?\d+)', indel_details)
	details, muts = {'I':0,'D':0,'C':0}, []
	for (letter,val) in cigar_toks:
		details[letter] = eval(val)
	if len(indel_toks) > 2 or (indel_type == '-' and len(indel_toks) > 1):
		mut_toks = re.findall(r'([MNDSI]+)(-?\d+)(\[[ATGC]+\])?', indel_toks[-1])
		for (letter,val,nucl) in mut_toks:
			if nucl == '':
				nucl = '[]'
			muts.append((letter, eval(val), nucl[1:-1]))
	try:	
		if indel_type[0] == '-':
			isize = 0
		else:
			isize = eval(indel_type[1:])
	except IndexError:
		return False

	return indel_type[0],isize,details, muts
	
def indelOutofGuideSeedPAM(indel):
	if indel[0] == '-':
		return True

	try:
		itype, isize, details, muts = tokFullIndel(indel)
		if itype != '-':
			left, right = details['L'], details['R']
			if (left > -23 and left < 5) or (right > -22 and right < 6):
				return False
		for (letter, pos, nucl) in muts:
			if letter == 'M' or letter == 'S':
				if pos > -22 and pos < 0 and not (letter == 'M' and pos == 0):
					return False
		return True
	except TypeError:
		return False

	
def create_summary(summary_file, indel_dict):
	with open(summary_file, 'w') as out_file:
		for i in range(1, 1252):
			oligo = 'Cpf1_Oligo'+str(i)
			out_file.write('@@@'+oligo+'\n')
			for indel in indel_dict[oligo]:
				for initial_oligo in indel_dict[oligo][indel]:
					out_file.write(indel + '\t' + initial_oligo + '\t' + str(indel_dict[oligo][indel][initial_oligo]) + '\n')
					
def make_results_directory(base_dir):
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
		
def run_all(args):
	make_results_directory(args.base_dir)
	make_results_directory(os.path.join(args.base_dir,'Figures'))
	summary_file = os.path.join(args.base_dir, 'summary_file.txt')
	indel_dict = parse_indelmap_results(args)
	create_summary(summary_file,  indel_dict)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Generates the summary file from the 1251 indelmap output filea (path to these may need to be changed). When a read is mapped to multiple oligo files, one is selected randomly out of suitable matches.')
	parser.add_argument('base_dir', help= 'Base directory for output files to be stored')
	parser.add_argument('indelmap_output', help = 'Path to the driectory where the mapped files are stored')
	args = parser.parse_args()

	run_all(args)