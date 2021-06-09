# CRISPR-Cpf1

Summer project undertaken July--September 2017 in the Parts group, Wellcome Sanger Institute. The aim was to analyse a dataset of genome sequences after CRISPR-Cpf1 activity. 

The original dataset was filtered (orig_trim_fasta_file.py) to ensure that each read had the correct structure for CRISPR-Cpf1 activity and each read was trimmed to reduce redundant information.
The cut reads were compared to the original genome sequence to determine which insertions and deletions took place.

The frequency and locations of insertions and deletions were determined in a probabilistic way.
