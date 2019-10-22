""""""

import argparse

if __name__ == "__main__":
	parser=argparse.ArgumentParser(prog='debruij.py',
description='Assembleur de séquence basé sur la méthode de Debruij.')
	parser.add_argument('-i',
action='store_true',
help='fichier fastq, single end')
	parser.add_argument('-k',
type=int,
default=21,
help='taille des kmer (optionnel - par defaut:21)')
	parser.add_argument('-r',
action='store_true',
default='',
	args = parser.parse_args()

def read_fastq(nom_fichier_fastq):
	

