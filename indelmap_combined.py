import subprocess

indelmap_path = 'C:\\Users\\jt20\\Documents\\Summer17\\indelmap_6.exe'

guide_index = 'C:\\Users\\jt20\\Documents\\Summer17\\guide_index.txt'

for i in range(1,1252):
	try:
		fastq_file = 'C:\\Users\\jt20\\Documents\\Summer17\\test_22_07\\SRR3710048_3\\Results_2\\Oligos\\Cpf1_Oligo'+str(i) + '\\all_reads.fastq'
		output_file = 'C:\\Users\\jt20\\Documents\\Summer17\\test_22_07\\SRR3710048_3\\Results_2\\mapped_against_null_11_09\\4_indelmap_output_'+str(i)+'.txt'
		oligo_file = 'C:\\Users\\jt20\\Documents\\Summer17\\SRR3709999\\filtered_fastq\\mapped_11_09\\'+str(i) + '.fasta'
		subprocess.call(indelmap_path +' '+ fastq_file +' '+ oligo_file+' ' +output_file +' '+ str(0) + ' 5 '+ guide_index, shell = True)
		print('-------------------')
		print str(i)
		print('-------------------')
	except IOError:
		print 'Error: '+ str(i)
		continue