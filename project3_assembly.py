import random 
import resource
import argparse
import time
from Bio import SeqIO
import re

def extract_sequences(filename):
    """
    Cette fonction prend en entrée le nom d'un fichier (au format FASTA ou FASTQ) et renvoie une liste de séquences 
    nucléotidiques en format texte.
    
    Args:
    - filename (str): Le nom du fichier à lire
    
    Returns:
    - sequences (list): Une liste de séquences nucléotidiques en format texte
    
    """
    
    # Initialisation de la liste de séquences
    sequences = []
    
    # Ouverture du fichier en mode lecture
    with open(filename, "r") as f:
        
        # Lecture des enregistrements (records) et extraction des séquences
        for record in SeqIO.parse(f, "fasta" if filename.endswith(".fa") else "fastq"):
            sequences.append(str(record.seq))
    
    # Retourne la liste de séquences extraites
    return sequences


def index_kmers(sequences, k):
    """
    Cette fonction prend une liste de séquences d'ADN et un entier k en entrée, et renvoie un dictionnaire
    où les clés sont les k-mers présents dans les séquences et les valeurs sont le nombre d'occurrences de chaque k-mer.

    Entré:
        sequences (list): une liste de séquences d'ADN (chaque séquence est une chaîne de caractères)
        k (int): la longueur de chaque k-mer à extraire des séquences

    Returns:
        dict: un dictionnaire où les clés sont les k-mers et les valeurs sont le nombre d'occurrences de chaque k-mer
    """

    kmers = {}  # Initialisation d'un dictionnaire vide pour stocker les k-mers et leur nombre d'occurrences
    #Compile une expression régulière pour valider que la séquence ne contient que les nucléotides A, C, G, et T
    pattern = re.compile('^[ACGT]+$')
    # Parcours de chaque séquence dans la liste des séquences
    for sequence in sequences:
        # Vérifie que la séquence ne contient que les nucléotides A, C, G, et T
        if not pattern.match(sequence):
            continue 
        # Parcours de chaque k-mer dans la séquence courante
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]  # Extraction du k-mer courant
            if kmer in kmers:  # Si le k-mer est déjà présent dans le dictionnaire, on incrémente son nombre d'occurrences
                kmers[kmer] += 1
            else:  # Sinon, on ajoute le k-mer au dictionnaire avec une occurrence initiale de 1
                kmers[kmer] = 1

    return kmers  # Renvoie le dictionnaire des k-mers et leur nombre d'occurrences

def filter_kmers(kmers, threshold):
    """Filtre les k-mers dont le seuil est inférieure à la valeur minimale spécifiée.
    
    Args:
        kmers (dict): un dictionnaire contenant des k-mers en clé et leur occurence en valeur.
        threshold (int): la valeur minimale d'occurence à conserver pour chaque k-mer.
        
    Returns:
        dict: un nouveau dictionnaire contenant uniquement les k-mers dont l'occurence est supérieure
              ou égale à la valeur minimale spécifiée.
    """
    filtered_kmers = {}
    # Parcourt tous les k-mers et leur abondance dans le dictionnaire 'kmers'.
    for kmer, abundance in kmers.items():
        # Vérifie si l'abondance du k-mer est supérieure ou égale à la valeur minimale spécifiée.
        if abundance >= threshold:
            # Ajoute le k-mer et son abondance correspondante au nouveau dictionnaire 'filtered_kmers'.
            filtered_kmers[kmer] = abundance
    # Retourne le nouveau dictionnaire 'filtered_kmers'.
    return filtered_kmers


def predecessors_and_successors(filtered_kmers):
    """
    Cette fonction prend en entrée un dictionnaire `kmer_index` qui contient les k-mers comme clés et les positions où 
    ils se produisent comme valeurs. La fonction renvoie un dictionnaire `index` qui contient pour chaque k-mer dans le
    `kmer_index`, une liste de ses prédécesseurs et de ses successeurs.

    :param kmer_index: un dictionnaire contenant les k-mers et leurs positions.
    :return: un dictionnaire contenant les k-mers avec leurs prédécesseurs et successeurs.
    """

    # Initialisation d'un dictionnaire vide pour contenir les prédécesseurs et successeurs de chaque k-mer
    index = {}

    # Boucle à travers chaque k-mer dans le dictionnaire `kmer_index`
    for w in filtered_kmers.keys():

        # Cherche les successeurs pour le k-mer actuel
        successors_l = []
        kmer = w[1:]
        posibility = [kmer + 'A', kmer + 'C', kmer + 'G', kmer +  'T']
        for i in posibility:
            if i in filtered_kmers.keys():
                successors_l.append(i)
        successors_index = successors_l

        # Cherche les prédécesseurs pour le k-mer actuel
        predecessors_l = []
        kmer = w[:-1]
        posibility = ['A' + kmer, 'C' + kmer,'G' + kmer,'T' + kmer]
        for i in posibility:
            if i in filtered_kmers.keys():
                predecessors_l.append(i)
        predecessors_index = predecessors_l

        # Ajoute les prédécesseurs et successeurs du k-mer actuel au dictionnaire `index`
        index[w] = {"successors": successors_index, "predecessors": predecessors_index}

    # Retourne le dictionnaire `index` contenant les prédécesseurs et successeurs de chaque k-mer
    return index


def simple_path(index):
    """
    Cette fonction prend en entrée un index et retourne une liste de contigs. Cette fonction construit les contigs à partir
    de l'index fourni en suivant les chemins simples dans le graphe de De Bruijn représenté par l'index.

    Args:
        index (dict): Un dictionnaire représentant l'index de De Bruijn. Chaque clé est un k-mer (une séquence
                      de longueur k) et chaque valeur est un dictionnaire avec les clés "predecessors" et "successors".
                      Les valeurs de ces clés sont des listes de k-mers prédécesseurs et successeurs du k-mer correspondant.

    Returns:
        list: Une liste de contigs construits à partir de l'index fourni.

    """
    
    # Créer un dictionnaire vide pour stocker les contigs construits et un ensemble vide pour les noeuds visités.
    contig_dict = {}
    visited = set()
    
    def build_contig(w):
        """
        Cette fonction récursive prend en entrée un k-mer w et retourne un contig construit à partir de ce k-mer
        en suivant les chemins simples dans le graphe de De Bruijn représenté par l'index fourni à la fonction.
        Le contig est stocké dans un dictionnaire pour éviter de recalculer les contigs déjà construits.

        Args:
            w (str): Un k-mer.

        Returns:
            str: Un contig construit à partir du k-mer w.

        """
        
        # Si le contig correspondant au k-mer w a déjà été construit, le retourner.
        if w in contig_dict:
            return contig_dict[w]
        
        # Sinon, construire un nouveau contig à partir du k-mer w.
        contig = w
        visited.add(w)
        successor = index[w]["successors"]
        
        # Tant que le k-mer a un seul successeur et que ce successeur n'a pas été visité,
        # ajouter la dernière lettre du successeur au contig et marquer le successeur comme visité.
        while len(successor) == 1 and successor[0] not in visited:
            contig += successor[0][-1]
            visited.add(successor[0])
            successor = index[successor[0]]["successors"]
        
        # Si le k-mer a plusieurs successeurs ou que le seul successeur restant a été visité,
        # explorer récursivement chaque successeur non visité et construire un contig à partir de chaque successeur.
        # Garder le contig le plus long construit à partir des successeurs.
        if len(successor) > 1:
            best_contig = contig
            for s in successor:
                if s not in visited:
                    new_contig = build_contig(s)
                    if len(new_contig) > len(best_contig):
                        best_contig = new_contig
            
            contig = best_contig
        
        # Faire de même avec les prédécesseurs du k-mer w, en ajoutant la première lettre de chaque prédécesseur au début du contig.

        pre = index[w]["predecessors"]
        while len(pre) == 1 and pre[0] not in visited:
            contig = pre[0][0] + contig
            visited.add(pre[0])
            pre = index[pre[0]]["predecessors"]
        
        if len(pre) > 1:
            best_contig = contig
            for p in pre:
                if p not in visited:
                    new_contig = build_contig(p) + contig
                    if len(new_contig) > len(best_contig):
                        best_contig = new_contig
            # Ajouter le contig construit pour le k-mer w au dictionnaire de contigs.
            contig = best_contig
        
        contig_dict[w] = contig
        return contig
    # Construire des contigs à partir de chaque k-mer dans l'index qui n'a pas été visité.
    list_contig = []
    for w in index.keys():
        if w not in visited:
            list_contig.append(build_contig(w))
    # Retourner la liste de contigs construits.
    return list_contig

def write_contigs_to_fasta(contigs, filename):
    """
    Cette fonction écrit une liste de séquences de contigs dans un fichier au format FASTA.

    Args:
        contigs (list): Une liste de séquences de contigs.
        filename (str): Le nom du fichier de sortie.

    Returns:
        None
    """
    with open(filename, 'w') as f:
        for i, contig in enumerate(contigs):
            f.write(f'>contig_{i+1}\n')
            f.write(f'{contig}\n')
    print(f'Les contigs ont été enregistrés dans le fichier {filename}.')
    
def main():
    parser = argparse.ArgumentParser(description="Assemble des contigs à partir d'un fichier FASTA ou FASTQ de séquences avec gestion greedy .")
    parser.add_argument("filename", help="Le nom du fichier contenant les séquences.")
    parser.add_argument("-k", "--kmer_length", type=int, default=21, help="La longueur des k-mers à indexer (par défaut : 21).")
    parser.add_argument("-a", "--threshold", type=int, default=2, help="La fréquence minimale à laquelle un k-mer doit apparaître pour être conservé (par défaut : 2).")
    parser.add_argument("-o", "--output_file", type=str, default="contigs.fasta", help="Le nom du fichier de sortie contenant les contigs assemblés (par défaut : 'contigs_greedy.fasta').")
    args = parser.parse_args()

    # Extraction des séquences à partir du fichier
    try:
        sequences = extract_sequences(args.filename)
    except FileNotFoundError:
        print(f"Erreur : le fichier '{args.filename}' n'existe pas.")
        return
    
    # Indexation des k-mers
    kmers = index_kmers(sequences, args.kmer_length)
    
    # Filtrage des k-mers à faible abondance
    filtered_kmers = filter_kmers(kmers, args.threshold)
    
    # Recherche des prédécesseurs et successeurs de chaque k-mer
    pred_et_succ = predecessors_and_successors(filtered_kmers)
    
    # Assemblage des contigs
    contigs = simple_path(pred_et_succ)
    
    # Écriture des contigs dans un fichier de sortie
    with open(args.output_file, 'w') as outfile:
        for i, contig in enumerate(contigs):
            outfile.write(f">Contig_{i+1}\n{contig}\n")

    elapsed_time = time.process_time()
    mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    
    
    
    # Affichage des résultats
    print(f"Nombre de séquences extraites : {len(sequences)}")
    print(f"Nombre de k-mers indexés : {len(kmers)}")
    print(f"Nombre de k-mers filtrés : {len(kmers) - len(filtered_kmers)}")
    print(f"Nombre de k-mers conservés : {len(filtered_kmers)}")
    print(f"Nombre de contigs assemblés : {len(contigs)}")
    print(f"Les contigs ont été ajoutés dans le fichier '{args.output_file}'.")
    print(f"{len(contigs)} contigs obtenu en  {elapsed_time:.2f} seconds, en utilisant {mem_usage:.2f}")

if __name__ == "__main__":
    main()
    start_time = time.time()