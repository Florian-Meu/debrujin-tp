""""""

###Import des modules
import argparse
import os

if __name__ == "__main__":
	parser=argparse.ArgumentParser(prog='debruij.py',
description='Assembleur de séquence basé sur la méthode de Debruij.')
	parser.add_argument('--i',
type=str,
help='fichier fastq, single end')
	parser.add_argument('--k',
type=int,
default=21,
help='taille des kmer (optionnel - par defaut:21)')
	parser.add_argument('--r',
type=str,
default='')
	args = parser.parse_args()

def read_fastq(nom_fichier_fastq):
	lire=False
	chemin="".join(["../data/",nom_fichier_fastq])
	with open(chemin,"r") as fastq :
		for ligne in fastq :
			if ligne.startswith("@"):
				lire=True
			elif ligne.startswith("+"):
				lire=False
			else :
				if lire :
					yield ligne

if __name__ == "__main__":
	nom_fichier_fastq=args.i
	for sequence in read_fastq(nom_fichier_fastq):
		print (sequence)
