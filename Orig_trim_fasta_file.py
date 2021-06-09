# Script to create a fastq file with trimmed reads 
import pickle
import pprint
mismatch_guides = ['CTGATGGTCCATGTCTGTTACTC','CTCCGCCCCACCTCCCCCAGCCC', 'CTCTCCTAGGACCCTCCCTCTCG', 'AAGAAGAAGACGAGCCTGAGTCC']

def write_to_fastq(fastq_out, line_id, trimmed_seq, first_line):
	if first_line:
		fastq_out.write(''.join(['@', line_id]))
		first_line = False
	else:
		fastq_out.write('\n' + ''.join(['@', line_id]))
	fastq_out.write('\n' + trimmed_seq)
	fastq_out.write('\n' + '+')
	fastq_out.write('\n' + 'I'*len(trimmed_seq))
	return first_line


def trim_sequence(line_id, line, count_codes, reads_used,   guide_list, first_line_fastq_mismatch):
	start_T_parameters = [66, 67, 65,68,64, 69, 63, 70, 62,71, 61,60]
	for start_T in start_T_parameters:
		if line[start_T:start_T+5] == 'TTTTT':
			break
		elif start_T == 60:
			start_index = line.find('AATTTCTACTCTTGTAGAT')
			start_end_flanking = line.find('TAACTAGATCT') - 9
			if start_index != -1 and start_end_flanking != -10:
				if -6 <= start_end_flanking - start_index - 19 <= 6:
					count_codes[3] += 1 
					return False, None, count_codes, first_line_fastq_mismatch
				elif 7 <= start_end_flanking - start_index -19  <= 27:
					count_codes[4] += 1
					return False, None, count_codes, first_line_fastq_mismatch
				else:
					count_codes[6] += 1
					return False, None, count_codes, first_line_fastq_mismatch
			else:
				count_codes[6] += 1
				return False, None, count_codes, first_line_fastq_mismatch
	
	start_index_parameters = [24, 23, 25,22, 26]
	for start_index in start_index_parameters:
		if line[start_index:start_index + 19] == 'AATTTCTACTCTTGTAGAT': 
			break
		elif start_index == 26: 
			for start_index in start_index_parameters:
				if line[start_index: start_index+20] == 'AATTTCTACTAAGTGTAGAT':
					count_codes[2] += 1
					return False, None, count_codes, first_line_fastq_mismatch
			count_codes[1] += 1
			start_index = 24

	end_index = line[start_index:None].find('TAACTAGATCT') + 11	+ start_index
	if end_index - (start_index+12) <= 70 or end_index - (start_index + 12) >= 150:
		count_codes[7] += 1
		return False, None, count_codes, first_line_fastq_mismatch
	else:	
		if end_index != 10 + start_index:
			if line[start_T+6:start_T+30].find('TTTA') != -1:
				if line[start_index + 19:start_index+42] in guide_list:
					count_codes[0] += 1
					guide = line[start_index + 19:start_index+42]
					#n if guide not in mismatch_guides:
					reads_used.write(line_id + '\t' + guide + '\t' + line[start_T+6:start_T+21]+'\n')
					return True, line[start_index+12:end_index], count_codes, first_line_fastq_mismatch
					# else:
						# mismatch_reads.write(line_id + '\t' + guide + '\t' + line[start_T + 6:start_T + 21] + '\n')
						# first_line_fastq_mismatch = write_to_fastq(mismatch_fastq_out, line_id, line[start_index+12:end_index], first_line_fastq_mismatch)
						# return False, None, count_codes, first_line_fastq_mismatch
				else:
				
					# Mixture of important/not important code to consider if the guide was used in the tests with the mismatch target sequences or one of those deemed to have a 
					# low read count or high mismatch frequency as well as attempting to cope with mutations in the scaffold or length of T's or barcode that result in the indices 
					# being misaligned. 
					
					for j in [20,18,21,17,22,16]:
						if line[start_index+j:start_index+j+23] in guide_list and line[start_index + j+23:start_index + j+ 28] == 'TTTTT':
							guide = line[start_index + j:start_index+j+23]
							reads_used.write(line_id + '\t' + guide + '\t' + line[start_T+6:start_T+21]+'\n')
							return True, line[start_index+12:end_index], count_codes, first_line_fastq_mismatch
					else:	
						return False, None, count_codes, first_line_fastq_mismatch
					
					if start_T - 19- start_index < 22 or start_T - 19- start_index > 24:
						count_codes[10] += 1 
						wrong_guide.write(line_id +'\t' +  line[start_index + 19:start_T] +'\t' + line[start_T+6:start_T+21]+'\n')
						return False, None, count_codes, first_line_fastq_mismatch
					elif start_T - 19- start_index  == 24 and line[start_index + 19] == 'T':
						count_codes[0] += 1
						guide = line[start_index+ 19:start_T]
						if guide not in mismatch_guides:
							reads_used.write(line_id +'\t' +  guide +'\t' + line[start_T+6:start_T+21]+'\n')
							return True, line[start_index + 12:end_index], count_codes,first_line_fastq_mismatch
						else:
							mismatch_reads.write(line_id + '\t' + guide + '\t' + line[start_T + 6: start_T + 21] + '\n')
							first_line_fastq_mismatch = write_to_fastq(mismatch_fastq_out, line_id, line[start_index+12:end_index], first_line_fastq_mismatch)
							return False, None, count_codes, first_line_fastq_mismatch
					else:
						count_codes[9] += 1
						guide = line[start_index+ 19:start_T]
						if guide not in mismatch_guides:
							reads_used.write(line_id +'\t' +  guide+'\t' + line[start_T+6:start_T+21]+'\n')
							return True, line[start_index + 12:end_index], count_codes, first_line_fastq_mismatch
						else:
							mismatch_reads.write(line_id + '\t' + guide + '\t' + line[start_T + 6: start_T + 21] + '\n')
							first_line_fastq_mismatch = write_to_fastq(mismatch_fastq_out, line_id, line[start_index+12:end_index], first_line_fastq_mismatch)
							return False, None, count_codes, first_line_fastq_mismatch
			else:
				count_codes[8] += 1
				return False, None, count_codes, first_line_fastq_mismatch
		else:
			count_codes[5] += 1
			return False, None, count_codes, first_line_fastq_mismatch

	
	
	
def	run(input_fasta, output_fastq, count_codes, reads_used,  guide_list):
	with open(input_fasta, 'r') as fasta_in:
		with open(output_fastq, 'w') as fastq_out:
			first_line_fastq = True
			first_line_unused = True
			first_line_fastq_mismatch = True
			for line in fasta_in:
				if line[0] == '>':
					line_id = line.strip()[1:None]
				else: 
					line = line.strip()
					trimmed, trimmed_seq, count_codes, first_line_fastq_mismatch = trim_sequence(line_id, line, count_codes,  reads_used,  guide_list, first_line_fastq_mismatch)
					if trimmed:		
						first_line_fastq = write_to_fastq(fastq_out, line_id, trimmed_seq, first_line_fastq)
	return count_codes

def run_all(number):
	count_codes = [0,0,0,0,0,0,0,0,0,0,0] 
	
	with open(r'C:\\Users\\jt20\\Documents\\Summer17\\data\\guide_list.pkl','r') as f:
		with open('C:\\Users\\jt20\\Documents\\Summer17\\data\\guide_list.txt','w') as g:
			guide_list = pickle.load(f)
			pprint.pprint(guide_list, stream =g)
	
	input_fasta = 'C:\\Users\\jt20\\Documents\\Summer17\\SRR3710048_all\\SRR3710048.' + str(number) + '.fasta'
	output_fastq = 'C:\\Users\\jt20\\Documents\\Summer17\\test_22_07\\SRR3710048_3\\updated\\' + str(number) + '.fastq'

	reads_used = open(''.join(['C:\\Users\\jt20\\Documents\\Summer17\\test_22_07\\SRR3710048_3\\updated\\reads_used_', str(number),'.txt']), 'w')
	
	count_codes = run(input_fasta, output_fastq, count_codes, reads_used,  guide_list)

if __name__ == '__main__':
	for i in range(1,22):
		print('--------------------------------')
		print i
		print('--------------------------------')
		run_all(i)
		