
# -*- coding: utf-8 -*-
"""Ce programme va permettre d'assembler des reads.

Pour se faire, la méthode des graphs de De Bruijn sera utilisée.
"""

### Import des modules
import argparse
import os
import statistics
import matplotlib
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


def path_average_weight(graph, chemin):
    """Cette fonction permet de retourner le poids moyen d'un
    chemin."""
    poids=0
    nbre_edges=0
    edges=graph.subgraph(chemin).edges(data=True)
    for u, v, e in edges:
        poids+=e['weight']
        nbre_edges+=1
    poids=poids/nbre_edges
    #print(poids)
    return poids

def remove_paths(graph, liste_chemins, delete_entry_node=False, delete_sink_node=False):
    """Fonction qui va permettre de nettoyer le graphique de chemins
    indésirables."""
    for chemin in liste_chemins:
        if delete_entry_node==False and delete_sink_node==False :
            for noeud in chemin[1:-1]:
                if noeud in graph.nodes:
                    graph.remove_node(noeud)
        elif delete_entry_node==True and delete_sink_node==False :
            for noeud in chemin[:-1]:
                if noeud in graph.nodes:
                    graph.remove_node(noeud)
        elif delete_entry_node==False and delete_sink_node==True :
            for noeud in chemin[1:]:
                if noeud in graph.nodes:
                    graph.remove_node(noeud)
        elif delete_entry_node==True and delete_sink_node==True :
            for noeud in chemin[:]:
                if noeud in graph.nodes:
                    graph.remove_node(noeud)
    return graph

def select_best_path(graph, ensemble_chemins, ensemble_longueurs, poids_moyen, delete_entry_node=False, delete_sink_node=False):
    """Fonction qui va permettre d'identifier le meilleur chemin."""
    while len(ensemble_chemins) > 1 :
        #print("ensemble_chemins = {}".format(ensemble_chemins))
        a_retirer = []
        taille_max=max(poids_moyen)
        chemin_fort_poids = []
        for i in range(len(ensemble_chemins)):
            if poids_moyen[i]==taille_max:
                chemin_fort_poids.append(ensemble_chemins[i])
        for chemin in chemin_fort_poids:
            ensemble_chemins.remove(chemin)
        a_retirer = ensemble_chemins + a_retirer
        if len(chemin_fort_poids)==1:
            #print("chemin_fort_poids = {}".format(chemin_fort_poids))
            graph = remove_paths(graph, a_retirer, delete_entry_node, delete_sink_node)
        else :
            longueur_max = max(ensemble_longueurs)
            grands_chemins = []
            for i in range(len(chemin_fort_poids)):
                if ensemble_longueurs[i]==longueur_max:
                    grands_chemins.append(chemin_fort_poids[i])
            for chemin in grands_chemins:
                chemin_fort_poids.remove(chemin)
            a_retirer = chemin_fort_poids + a_retirer
            if len(grands_chemins)==1:
                #print("grands_chemins = {}".format(grands_chemins))
                graph = remove_paths(graph, a_retirer, delete_entry_node, delete_sink_node)
            else :
                random.seed(9001)
                choix=random.randint(0,len(grands_chemins))
                ensemble_chemins=grands_chemins[choix]
                grands_chemins.remove(ensemble_chemins)
                a_retirer = grands_chemins + a_retirer
                #print("ensemble_chemins = {}".format(ensemble_chemins))
                graph = remove_paths(graph, a_retirer, delete_entry_node, delete_sink_node)
    return graph

def find_bubbles(graphique):
    """Fonction qui va permettre de trouver les bulles au sein
    de l'arbre. Elle retournera les origines et les fins de ces
    bulles."""
    bulles = []
    ensemble_noeuds = list(graphique.nodes)
    #print(ensemble_noeuds)
    for i in range(len(ensemble_noeuds)):
        #print(ensemble_noeuds)
        if len(list(graphique.successors(ensemble_noeuds[i])))>1:
            debut = ensemble_noeuds[i]
            j = 0
            fin = ""
            while len(list(graphique.predecessors(ensemble_noeuds[i+j])))==1 and j<len(ensemble_noeuds)-i:
                j+=1
                if len(list(graphique.predecessors(ensemble_noeuds[i+j])))>1:
                    fin = ensemble_noeuds[i+j]
            if fin!="":
                bulles.append([debut, fin])
    return bulles      

def solve_bubble(graphique,debut,fin):
    """Fonction permettant de nettoyer le graphique de résoudre
    la bulle contenue entre les deux bornes."""
    ensemble_chemins = []
    for path in nx.all_simple_paths(graphique,\
            source = debut, target = fin):
                ensemble_chemins.append(path)
    poids_moyen = []
    for chemin in ensemble_chemins:
        poids_moyen.append(path_average_weight(graphique,chemin))
    #ecarttype = std(poids_moyen)
    ensemble_longueurs = []
    for chemin in ensemble_chemins:
        ensemble_longueurs.append(len(chemin))
    graphique = select_best_path(graphique, ensemble_chemins, ensemble_longueurs, poids_moyen, delete_entry_node=False, delete_sink_node=False)
    return graphique

def simplify_bubbles(graphique):
    """Fonction permettant de nettoyer le graphique de toutes
    les éventuelles bulles présentes dans le graphique."""
    liste_bulles = find_bubbles(graphique)
    #print(liste_bulles)

    for bulle in liste_bulles:
        graphique = solve_bubble(graphique,bulle[0],bulle[1])

    return graphique

def solve_entry_tips(graphique, entrees):
    """Fonction qui permet d'enlever les entrées indésirables"""
    #établissement des bornes de chemins d'entrée.
    bornes_initiales = []
    for noeuds_entree in entrees:
        ensemble_noeuds = list(graphique.nodes)
        for i in range(len(ensemble_noeuds)):
            if len(list(graphique.predecessors(ensemble_noeuds[i])))>1:
                fin = ensemble_noeuds[i]
        bornes_initiales.append([noeuds_entree,fin])
    #Définitions des chemins et valeurs associées
    ensemble_chemins = []
    poids_moyen = []
    ensemble_longueurs = []
    for borne in bornes_initiales:
        for path in nx.all_simple_paths(graphique,\
                source = borne[0], target = borne[1]):
            ensemble_chemins.append(path)
            poids_moyen.append(path_average_weight(graphique,path))
            ensemble_longueurs.append(len(path))
    #Nettoyage du graphique
    graphique = select_best_path(graphique, ensemble_chemins, ensemble_longueurs, poids_moyen, delete_entry_node=True, delete_sink_node=False)
    return graphique

def solve_out_tips(graphique, sorties):
    """Fonction qui permet d'enlever les entrées indésirables"""
    #établissement des bornes de chemins de sortie.
    bornes_initiales = []
    for noeuds_sortie in sorties:
        ensemble_noeuds = list(graphique.nodes)
        for i in range(len(ensemble_noeuds)):
            if len(list(graphique.successors(ensemble_noeuds[i])))>1:
                debut = ensemble_noeuds[i]
        bornes_initiales.append([debut,noeuds_sortie])
    #Définitions des chemins et valeurs associées
    ensemble_chemins = []
    poids_moyen = []
    ensemble_longueurs = []
    for borne in bornes_initiales:
        for path in nx.all_simple_paths(graphique,\
                source = borne[0], target = borne[1]):
            ensemble_chemins.append(path)
            poids_moyen.append(path_average_weight(graphique,path))
            ensemble_longueurs.append(len(path))
    #Nettoyage du graphique
    graphique = select_best_path(graphique, ensemble_chemins, ensemble_longueurs, poids_moyen, delete_entry_node = False, delete_sink_node = True)
    return graphique

###Définition de la fonction Main
def main():
    """La fonction main() correspond à ce qui sera executé si le code
    source est lancé"""
    parser = argparse.ArgumentParser(prog='debruij.py',\
    description='Assembleur de séquence basé sur la méthode de De Bruij.')
    parser.add_argument('--i', type = str, help='fichier fastq, single end')
    parser.add_argument('--k', type = int, default = 21,\
    help='taille des kmer (optionnel - par defaut : 21)')
    parser.add_argument('--r', type = str, default='',\
    help='genome de reference (optionnel)')
    parser.add_argument('--o', type = str, default='',\
    help='fichier contig (optionnel)')
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

    graphique = simplify_bubbles(graphique)

    #recherche des entrées indésirables.
    noeuds_entree = get_starting_nodes(graphique)
    if len(noeuds_entree)>0:
        graphique = solve_entry_tips(graphique, noeuds_entree)
    
    noeuds_terminaux = get_sink_nodes(graphique)
    if len(noeuds_terminaux)>0:
        graphique = solve_out_tips(graphique, noeuds_terminaux)

    #nx.draw_networkx(graphique)

###Si fichier lancé on execute la boucle main.
if __name__ == "__main__":
    main()
