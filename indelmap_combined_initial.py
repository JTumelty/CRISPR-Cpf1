import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description = 'Maps the day 0 filtered reads aginst the original template oligo file with cut site at 18th nucleotide after the PAM.')
parser.add_argument('indelmap_path', help = 'Path to indelmap.exe') # 'C:\\Users\\jt20\\Documents\\Summer17\\indelmap_6.exe'
parser.add_argument('guide_index', help = 'Path to file containing the start and end of the guides in the read') # 'C:\\Users\\jt20\\Documents\\Summer17\\guide_index.txt'
parser.add_argument('oligo_file', help = 'Path to the template oligo file') #'C:\\Users\\jt20\\Documents\\Summer17\\test_22_07\\orig_oligo_template.fasta'
parser.add_argument('fastq_dir', help = 'Path to the directory containing the trimmed and filtered SRR3709999  fastq files')
parser.add_argument('output_dir', help = 'Path to directory where mapped files should be stored')
args = parser.parse_args()

for i in range(1,21):
	try:
		fastq_file = os.path.join(args.fastq_dir, 'SRR3709999.'+str(i)+'.fastq')
		output_file = os.path.join(args.output_dir,str(i)+'_mapped.txt')
		
		subprocess.call(args.indelmap_path +' '+ fastq_file +' '+ args.oligo_file+' ' +output_file +' '+ str(1) + ' 5 '+ args.guide_index, shell = True)
		print('-------------------')
		print str(i)
		print('-------------------')
	except IOError:
		print 'Error: '+ str(i)
		continue
		
		