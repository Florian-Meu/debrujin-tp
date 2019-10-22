"""Ce programme va permettre d'assembler des reads.

Pour se faire, la méthode des graphs de De Bruijn sera utilisée.
"""

### Import des modules
import argparse
import os
import statistics
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
    graphique = nx.DiGraph()
    for kmer, poids in dico_kmers.items():
        graphique.add_edge(kmer[:-1], kmer[1:], weight=poids)
    return graphique


def get_starting_nodes(graphique):
    """Fonction qui permet de relever les noeuds d'entrée"""
    noeuds_entree = []
    for noeud in graphique.nodes:
        if len(list(graphique.predecessors(noeud)))==0:
            noeuds_entree.append(noeud)
    return noeuds_entree

def get_sink_nodes(graphique):
    """Fonction qui permet de relever les noeuds de sortie"""
    noeuds_sortie = []
    for noeud in graphique.nodes:
        if len(list(graphique.successors(noeud))) == 0:
            noeuds_sortie.append(noeud)
    return noeuds_sortie

def get_contigs(graphique, debuts, fins):
    """Fonction permettant de générer une liste de tulpes
    contenant les contigs associés à leur taille."""
    contigs = []
    for noeud_depart in debuts:
        for noeud_fin in fins:
            for path in nx.all_simple_paths(graphique,\
            source=noeud_depart, target = noeud_fin):
                prep_contig = path
                contig_ecrit = []
                contig_ecrit.append(prep_contig[0])
                for i in range(1, len(prep_contig)):
                    contig_ecrit.append(prep_contig[i][-1:])
                contig_ecrit = "".join(contig_ecrit)
                contigs.append((contig_ecrit, len(contig_ecrit)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(liste_contigs, nom_fichier):
    with open(nom_fichier, "w") as fichier_sortie:
        numero = 0
        for contigs in liste_contigs:
            fichier_sortie.write(">contig_{0} len={1}\n".format(numero, contigs[1]))
            fichier_sortie.write("{0}\n".format(fill(contigs[0])))
            numero += 1

def std(liste_valeurs):
    """Calcul l'écart-type de la liste de valeurs"""
    return statistics.stdev(liste_valeurs)


def path_average_weight():
    pass

def remove_paths(graph, delete_entry_node=False, delete_sink_node=False):
    pass

def select_best_path():
    pass

def solve_bubble(graphique):
    """Fonction permettant de nettoyer le graphique de tous
    les éventuels noeuds présents dans le graphique."""
    debuts = get_starting_nodes(graphique)
    fins = get_sink_nodes(graphique)
    ensemble_chemins = get_contigs(graphique)
    poids_moyen=[]

    for chemin in ensemble_chemins:
        poids_moyen.append(path_average_weight(graphique,chemin))

    ecarttype = std(poids_moyen)

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
    description='Assembleur de séquence basé sur la méthode de De Bruij.')
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

    debuts = get_starting_nodes(graphique)
    #print(debuts)

    fins = get_sink_nodes(graphique)
    #print(fins)

    liste_contigs = get_contigs(graphique, debuts, fins)
    #print(liste_contigs)

    save_contigs(liste_contigs, "Export_contigs.fna")


###Si fichier lancé on execute la boucle main.
if __name__ == "__main__":
    main()
