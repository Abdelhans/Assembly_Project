# Projet d'Assemblage basé sur un Graphe de De Bruijn

## Description
Ce projet a pour objectif de développer un assembleur basé sur un graphe de De Bruijn pour assembler plusieurs jeux de données, tout en abordant divers défis tels que la gestion des erreurs de séquençage et des répétitions. L'assemblage produit est ensuite évalué pour sa qualité à l'aide de l'outil QUAST.

## Objectifs
1. Écrire un assembleur utilisant un graphe de De Bruijn.
2. Assembler différents ensembles de données tout en gérant les erreurs de séquençage.
3. Évaluer la qualité des assemblages produits à l'aide de QUAST.

## Données disponibles
Les données de séquençage et les séquences de référence nécessaires pour ce projet peuvent être téléchargées à partir du lien suivant :
- [Données de séquençage](https://nextcloud.univ-lille.fr/index.php/s/Ao8kz6iXkA83mJR)

## Méthodologie
1. **Indexation des k-mers** : Les k-mers sont extraits des séquences et indexés à l'aide d'une structure appropriée.
2. **Construction du Graphe de De Bruijn** : Le graphe est construit à partir des k-mers indexés pour permettre la navigation et l'assemblage.
3. **Assemblage des Contigs** : Les contigs sont assemblés en suivant des chemins simples dans le graphe, en gérant les conflits entre les différents chemins possibles.
4. **Évaluation de l'Assemblage** : La qualité des assemblages est évaluée en utilisant QUAST pour comparer les contigs assemblés aux séquences de référence.

## Installation
Pour exécuter le projet, assurez-vous d'avoir Python et les bibliothèques suivantes installées :
- Biopython
- tqdm

Vous pouvez installer les dépendances avec pip :
```bash
pip install biopython tqdm
```
## Utilisation 

Pour exécuter le script d'assemblage, utilisez la commande suivante :
```bash
python project3_assembly.py --filename <chemin/vers/fichier.fasta> -k <taille_kmer> -o <fichier_sortie.fasta>
```

## Paramètres

--filename : Chemin vers le fichier contenant les séquences (au format FASTA ou FASTQ).

-k : (optionnel) Taille des k-mers à utiliser (par défaut : 21).

-o : (optionnel) Nom du fichier de sortie pour les contigs assemblés (par défaut : 'contigs.fasta').


