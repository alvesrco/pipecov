'''
Author: Renato Oliveira
Date: 24-01-2022
Version: 1.0


1.Accepts a alignment fasta file
2.Extracts the first sequence
3.Deletes initials and finals gaps '-' and replaces interior gaps by 'N'

Usage: python correct_frameshift.py alignment.fasta output
'''

import sys
import os

def read_correct_extract(arq, output_name):
	novoFasta=""
	novoFasta+=arq.readline()

	for line in arq:
		
		if(line[0]!=">"):
			begin=0
			while(line[begin]=="-"):
				begin+=1

			end=len(line)-2
			while(line[end]=="-"):
				end-=1
			end+=1

			for i in range(begin, end):
				if(line[i] == "-"):
					novoFasta+="N"
				else:
					novoFasta+=line[i].upper()

		else:
			break
				

	novoArq = open(output_name+"_corrected.fasta", "w")
	novoArq.write(novoFasta)

	novoArq.close()


input_file = open(sys.argv[1], "r")
output_name = sys.argv[2]

read_correct_extract(input_file, output_name)

input_file.close()
