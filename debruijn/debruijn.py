"""Ce programme va permettre d'assembler des reads.

Pour se faire, la méthode des graphs de Debruijn sera utilisée.
"""

###Import des modules
import argparse
#import os

def main():
    """La fonction main() correspond à ce qui sera executé si le code source est lancé"""
    parser = argparse.ArgumentParser(prog='debruij.py',\
    description='Assembleur de séquence basé sur la méthode de Debruij.')
    parser.add_argument('--i', type = str, help='fichier fastq, single end')
    parser.add_argument('--k', type = int, default = 21, help='taille des kmer (optionnel - par defaut:21)')
    parser.add_argument('--r', type = str, default='')
    args = parser.parse_args()

    liste_reads = []
    for sequence in read_fastq(args.i):
        print(sequence)
        liste_reads.append(sequence)
    
    liste_kmers=[]
    for read in liste_reads:
        for kmers in cut_kmer(read,args.k):
            liste_kmers.append(kmers)
    print(liste_kmers)
        



def read_fastq(nom_fichier_fastq):
    """Cette fonction a pour but de récupérer les reads contenus dans un fichier Fastq."""
    lire = False
    chemin = "".join(["../data/", nom_fichier_fastq])
    with open(chemin, "r") as fastq:
        for ligne in fastq:
            if ligne.startswith("@"):
                lire = True
            elif ligne.startswith("+"):
                lire = False
            else:
                if lire:
                    yield ligne

def cut_kmer(sequence,taille):
    """Cette fonction va permettre d'obtenir des k-mers de taille désirée à
    partir d'un read renseigné."""
    for indice in range(len(sequence[:-taille])) :
        yield sequence[indice:indice+taille]

if __name__ == "__main__":
    main()
