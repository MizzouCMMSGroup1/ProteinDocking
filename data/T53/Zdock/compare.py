from Bio.PDB import *
import sys

p=PDBParser()
s=p.get_structure('B', 'complex.1.pdb')
atom_list=Selection.unfold_entities(s,'A')

ns=NeighborSearch(atom_list)

for atom in atom_list:
    center=atom.get_coord()
    R=10.0

    neighbor_l=ns.search(center,R)
    residue_list=Selection.unfold_entities(neighbor_l,'R')

    for r in residue_list:
        print r



# center=atom_list[0].get_coord()
# R=5.0
# 
# neighbor_l=ns.search(center,R)
# residue_list=Selection.unfold_entities(neighbor_l,'R')
# 
# for r in residue_list:
#     print r

#print residue_list