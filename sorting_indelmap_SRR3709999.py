from update_cut_site_number import update_cut_site_number, consider_mismatch, consider_indel
import argparse
import os

def create_dicts(args):
	oligo_dict = {}
	read_id_seq_dict = {}

	for i in range(1,21):
		oligo_dict = parse_initial_mapped(oligo_dict,i,args.mapped_files)
		read_id_seq_dict = parse_fastq(read_id_seq_dict, i, args.fastq_files)					
	return read_id_seq_dict, oligo_dict				

def parse_initial_mapped(oligo_dict, i, mapped_files)
	with open(os.path.join(mapped_files, str(i)+'_mapped.txt'),'r') as in_file:
		for line in in_file:
			if len(line.strip('\n')) > 0:
				if line.strip('\n').split('\t')[1] == 'None':
					continue
				oligo = line.strip('\n').split('\t')[1].split('_')[0]
				if oligo in oligo_dict:
					oligo_dict[oligo][0] += 1
					if (line.strip('\n').split('\t')[2], line.strip('\n').split('\t')[3]) not in oligo_dict[oligo][1]:
						oligo_dict[oligo][1][(line.strip('\n').split('\t')[2], line.strip('\n').split('\t')[3])] = [1,line.strip('\n').split('\t')[0]]
					else:
						oligo_dict[oligo][1][(line.strip('\n').split('\t')[2], line.strip('\n').split('\t')[3])][0] += 1
				else:
					oligo_dict[oligo]=[1,{(line.strip('\n').split('\t')[2], line.strip('\n').split('\t')[3]):[0, line.strip('\n').split('\t')[0]]}]
	return oligo_dict
	
def parse_fastq(read_id_seq_dict, i, fastq):
	with open('C:/Users/jt20/Documents/Summer17/SRR3709999/filtered_fastq/trimmed/SRR3709999.'+str(i)+'.fastq','r') as fastq_file:
		for line in fastq_file:
			if len(line.strip('\n'))>0:
				line_id = line.strip('\n')
				if line_id[0] == '@':
					read = fastq_file.next().strip('\n')
					read_id_seq_dict[line_id[1:None]] = read
	return read_id_seq_dict
								
def output_template_files(read_id_seq_dict, oligo_dict, cut_site):
	for i in range(1,1252):
		with open('C:/Users/jt20/Documents/Summer17/SRR3709999/filtered_fastq/mapped_11_09/'+str(i)+'.fasta','w') as oligo_file:
			read_count = oligo_dict['Oligo'+str(i)][0]
			for (indel, mut) in oligo_dict['Oligo'+str(i)][1]:
				indel_count = oligo_dict['Oligo'+str(i)][1][(indel,mut)][0]
				perc = indel_count/float(read_count)*100
				if perc < 0.5:
					continue
				cut_site_number = update_cut_site_number(cut_site+60,indel, mut)
				oligo_file.write('>Oligo_'+str(i)+'_'+indel+':'+mut+':'+str(perc) + ' ' + str(cut_site_number) + ' FORWARD\n')
				oligo_file.write(read_id_seq_dict[oligo_dict['Oligo'+str(i)][1][(indel,mut)][1]]+'\n')
			
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Produces the template fasta files considering the initial mutations and indels at the start. ')
	parser.add_argument('mapped_files', help = 'Path to the directory containing the text files with the indelmap output from SRR3709999.')
	parser.add_argument('fastq_files', help = 'Path to the directory containing the SRR3709999 fastq files')
	parser.add_argument('--cut_site', help = 'Position of the cut site in a null oligo. In the template file, the second column in the description line would be the cut site plus 60. Default = 18', default = 18, type = int)
	args = parser.pare_args()
	
	read_id_seq_dict, oligo_dict = create_dicts(args)
	output_template_files(read_id_seq_dict, oligo_dict, args.cut_site)