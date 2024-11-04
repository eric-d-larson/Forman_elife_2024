#!/usr/bin/env python
import sys, re, gzip
from Bio import SeqIO

def split_to_chr_and_replace_nonACGTN(infasta, path):

    file_path = str(path)
    file_path_final = ''
	
    if not file_path.endswith('/'):
	    file_path_final = file_path + '/'
    else:
	    file_path_final = file_path
	
    for rec in SeqIO.parse(open(infasta), 'fasta'): 
	    with gzip.open((file_path_final + str(rec.id) +'.fa.gz'), 'wt') as f:
		    f.write('>'+rec.id+'\n')
		    f.write(re.sub('[^ACGTN]', 'N', str(rec.seq).upper())+'\n')
	    f.close()

if __name__ == "__main__":
	infasta, path = sys.argv[1:]
	split_to_chr_and_replace_nonACGTN(infasta, path)
