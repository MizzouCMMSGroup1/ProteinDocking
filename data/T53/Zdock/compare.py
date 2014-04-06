from Bio.PDB import *
import sys

parser=PDBParser()
structure=parser.get_structure('model', 'complex.1.pdb')

model=structure[0]
chainA = model['A']
chainB = model['B']

atom_list1=Selection.unfold_entities(chainA,'A')
atom_list2=Selection.unfold_entities(chainB,'A')

ns = NeighborSearch(atom_list2)

for atom in atom_list1:
    center=atom.get_coord()
    R=5.0

    neighbor_l = ns.search(center,R)
    residue_list = Selection.unfold_entities(neighbor_l,'R')

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