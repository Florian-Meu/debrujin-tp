"""Ce programme va permettre d'assembler des reads.

Pour se faire, la méthode des graphs de Debruijn sera utilisée.
"""

### Import des modules
import argparse
#import os
import networkx as nx

### Liste des fonctions.
def read_fastq(fichier_fastq):
    """Cette fonction a pour but de récupérer les reads contenus
    dans un fichier Fastq."""
    lire = False
    with open(fichier_fastq, "r") as fastq:
        for ligne in fastq:
            if ligne.startswith("@"):
                lire = True
            elif ligne.startswith("+"):
                lire = False
            else:
                if lire:
                    yield ligne[:-1]

def cut_kmer(sequence, taille_kmer):
    """Cette fonction va permettre d'obtenir des k-mers de taille désirée à
    partir d'un read renseigné."""
    for indice in range(len(sequence[:-taille_kmer+1])):
        yield sequence[indice:indice+taille_kmer]

def build_kmer_dict(fichier_fastq, taille_kmer):
    """Cette fonction va permettre de calculer les occurrences de
    chaque Kmers contenus au sein des reads issus du fastq"""
    liste_reads = []
    for sequence in read_fastq(fichier_fastq):
        liste_reads.append(sequence)
    occurrence_kmers = {}
    for read in liste_reads:
        for kmer in cut_kmer(read, taille_kmer):
            if kmer in occurrence_kmers.keys():
                occurrence_kmers[kmer] += 1
            else:
                occurrence_kmers[kmer] = 1
    return occurrence_kmers

def build_graph(dico_kmers):
    """Cette fonction va permettre de créer un digraph qui permettra,
    à terme, d'aligner les reads"""
    graphique=nx.DiGraph()
    for kmer,poids in dico_kmers.items() :
        graphique.add_edge(kmer[:-1], kmer[1:], weight=poids)
    return graphique       
    

def get_starting_nodes():
    pass

def get_sink_nodes():
    pass

def get_contigs():
    pass

def save_contigs():
    pass

def std():
    pass

def path_average_weight():
    pass

def remove_paths():
    pass

def select_best_path():
    pass

def solve_bubble():
    pass

def simplify_bubbles():
    pass

def solve_entry_tips():
    pass

def solve_out_tips():
    pass



###Définition de la fonction Main
def main():
    """La fonction main() correspond à ce qui sera executé si le code
    source est lancé"""
    parser = argparse.ArgumentParser(prog='debruij.py',\
    description='Assembleur de séquence basé sur la méthode de Debruij.')
    parser.add_argument('--i', type = str, help='fichier fastq, single end')
    parser.add_argument('--k', type = int, default = 21,\
    help='taille des kmer (optionnel - par defaut:21)')
    parser.add_argument('--r', type = str, default='')
    args = parser.parse_args()

    liste_reads = []
    for sequence in read_fastq(args.i):
        #print(sequence)
        liste_reads.append(sequence)

    liste_kmers = []
    for read in liste_reads:
        for kmers in cut_kmer(read, args.k):
            liste_kmers.append(kmers)
    #print(liste_kmers)

    occurrence_kmers = build_kmer_dict(args.i, args.k)
    #print(occurrence_kmers)

    graphique = build_graph(occurrence_kmers)

###Si fichier lancé on execute la boucle main.
if __name__ == "__main__":
    main()
