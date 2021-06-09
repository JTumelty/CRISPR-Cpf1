import argparse

def construct_oligo_indel_dict(args):
	oligo_indel_dict = {}
	try:
		in_file =  open(args.summary_indels,'r')
	except AttributeError:
		in_file =  open(args.summary_indels_percent,'r') 
			
	for line in in_file:
		if len(line.strip()) == 0:
			continue
		if line.strip()[0] == '@':
			oligo = eval(line.strip().split('Cpf1_Oligo')[1])
			oligo_indel_dict[oligo] = [0,{}]
			continue
		info = line.strip().split('\t')
		indel= info[0]
		oligo_indel = info[1]
		count = info[2]
		oligo_indel_dict[oligo][0] += eval(count)
		if indel in oligo_indel_dict[oligo][1]:
			oligo_indel_dict[oligo][1][indel][oligo_indel] = eval(count)
		else:
			oligo_indel_dict[oligo][1][indel] = {oligo_indel:eval(count)}
			
	in_file.close()
	return oligo_indel_dict

def output_summary_indels_percent(args, oligo_indel_dict):			
	with open(args.summary_indels_percent,'w') as out_file:
		for i in range(1,1252):
			oligo = i
			out_file.write('@@@Cpf1_Oligo'+str(oligo) +'\n')
			for indel in oligo_indel_dict[oligo][1]:
				for oligo_indel in oligo_indel_dict[oligo][1][indel]:
					out_file.write(indel + '\t' + oligo_indel + '\t' + str(oligo_indel_dict[oligo][1][indel][oligo_indel]) + '\t' + str(oligo_indel_dict[oligo][1][indel][oligo_indel]/float(oligo_indel_dict[oligo][0])) + '\n')
				
def run_all(args):	
	oligo_indel_dict = construct_oligo_indel_dict(args)
	output_summary_indels_percent(args, oligo_indel_dict)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Outputs the summary indels file with an extra column with the percentage of times that indel occurred for each guide')
	parser.add_argument('-in','--summary_indels', help = 'Path to input summary file', default = 'C:/Users/jt20/Documents/Summer17/Final_results/summary_file.txt', type = str)
	parser.add_argument('-out','--summary_indels_percent', help = 'Path for the output summary file with percents to be stored', default = 'C:/Users/jt20/Documents/Summer17/Final_results/summary_file_percent.txt', type =str)
	args = parser.parse_args()

	run_all(args)