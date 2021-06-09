import re
import argparse
def update_cut_site_number(cut_site_number,indel,mismatch):
	cut_site_number = consider_mismatch(mismatch, cut_site_number)
	cut_site_number = consider_indel(indel, cut_site_number)
	return cut_site_number
	
def consider_mismatch(mismatch, cut_site_number):
	for i in range(len(mismatch)):
		if mismatch[i] == 'D':
			for j in range(i,len(mismatch)):
				if mismatch[j] == 'S':
					deletion_size = int(mismatch[i+1:j])
					break
			for k in range(j+1,len(mismatch)):
				if mismatch[k] in ['M','I','D']:
					del_pos = int(mismatch[j+1:k])
					break
				elif k == len(mismatch) -1:
					del_pos = int(mismatch[j+1:])
					break
			if del_pos < 0:
				cut_site_number -= deletion_size
		if mismatch[i] == 'I':
			for j in range(i,len(mismatch)):
				if mismatch[j] == 'S':
					insertion_size = int(mismatch[i+1:j])
					break
			for k in range(j+1,len(mismatch)):
				if mismatch[k] in ['M','I','D']:
					ins_pos = int(mismatch[j+1:k])
					break
				elif k == len(mismatch) -1:
					ins_pos = int(mismatch[j+1:])
					break
			if ins_pos < 0:
				cut_site_number += insertion_size	
	return cut_site_number

def consider_indel(indel, cut_site_number):
	if indel == '-':
		return cut_site_number
	else:
		positions = re.split('[LRCID]', indel.split('_')[1])
		left_pos = int(positions[1])
		right_pos = int(positions[2])
		first_indel_length = int(indel.split('_')[0][1:])
		
		if left_pos >= 0:
			return cut_site_number
		else:
			if indel[0] == 'D':
				if indel.split('_')[1].find('I') == -1 and indel.split('_')[1].find('C') == -1:
					return cut_site_number - first_indel_length
				elif indel.split('_')[1].find('I') == -1:
					return cut_site_number - first_indel_length # may need to edit
				else:
					i_index = indel.split('_')[1].find('I')
					for j in range(i_index + 1,len(indel.split('_')[1])):
						if indel.split('_')[1][j] in ['D','C','R','I']:
							break
					insertion_length = int(indel.split('_')[1][i_index+1:j])
					if insertion_length == first_indel_length:
						return cut_site_number
					else:
						net_size = insertion_length - first_indel_length
						return cut_site_number + net_size # may want to edit
							
			else:
				if indel.split('_')[1].find('D') == -1 and indel.split('_')[1].find('C') == -1:
					return cut_site_number + first_indel_length
				elif indel.split('_')[1].find('D') == -1:
					return cut_site_number + first_indel_length # may need to edit
				else:
					i_index = indel.split('_')[1].find('D')
					for j in range(i_index + 1,len(indel.split('_')[1])):
						if indel.split('_')[1][j] in ['D','C','R','I']:
							break
					deletion_length = int(indel.split('_')[1][i_index+1:j])
					if deletion_length == first_indel_length:
						return cut_site_number
					else:
						net_size = first_indel_length - deletion_length
						return cut_site_number + net_size # may want to edit
				

		
		
