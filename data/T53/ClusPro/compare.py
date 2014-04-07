from Bio.PDB import *
import sys

def parse_pdb_surfaces(name='input.pdb', search_radius=5.0):

    parser=PDBParser()
    structure=parser.get_structure('model', name)

    model=structure[0]
    chainA = model['A']
    chainB = model['B']

    atom_list1=Selection.unfold_entities(chainA,'A') #A is for atom, not chain id
    atom_list2=Selection.unfold_entities(chainB,'A')

    ns = NeighborSearch(atom_list2)

    surface_residues = []

    for atom in atom_list1:
        center=atom.get_coord()
        R=search_radius

        neighbor_l = ns.search(center,R)
        residue_list = Selection.unfold_entities(neighbor_l,'R')

        for r in residue_list:
            surface_residues += [(r.get_resname(), r.id[1])]

    unique_surface_residues = set(surface_residues)
    return unique_surface_residues



def main():
    for i in range(1,20+1):
        print "testing protein: ", i
        residues = parse_pdb_surfaces('model.000.'+("%02d" % i)+'.pdb', 5.0)
        
        print "surface residues: ", len(residues)
        for usr in sorted(residues, key=lambda x: x[1]):
            print usr

if __name__ == '__main__':
	main()

