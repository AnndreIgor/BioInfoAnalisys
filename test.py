import os
from main import construir_arvores
from pathlib import Path

from collections import Counter
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

with open(os.path.join('data/out/tmp/PLASMODIUM8.aln'), "r") as handle:
    alignment = AlignIO.read(handle, "clustal")

# Supondo que 'alignment' seja um objeto Bio.Align.MultipleSeqAlignment
sequence_names = [record.id for record in alignment]
duplicates = [item for item, count in Counter(sequence_names).items() if count > 1]

if duplicates:
    print("Nomes duplicados encontrados:", duplicates)

    for i, record in enumerate(alignment):
        record.id = f"seq_{i}"

calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment) 

constructor = DistanceTreeConstructor()

tree = constructor.nj(distance_matrix)

# Salva a árvore
path_o_tree = os.path.join("data",
                           "out", 
                           "Trees", 
                           f'tree{2222}')
# 'tree_{Path('PLASMODIUM8.aln').stem}.{'nexus'}' 
Phylo.write(tree, path_o_tree, 'nexus')
