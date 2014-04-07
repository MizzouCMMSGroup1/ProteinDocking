from Bio.PDB import *
import sys
import argparse

def parse_pdb_surfaces(folder='folder', name='input.pdb', search_radius=5.0):

    parser=PDBParser()
    structure=parser.get_structure('model', folder + '/' + name)

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


def parse_zdock(input_radius=5.0):
    for i in range(1,10+1):
        print "testing protein: ", i
        residues = parse_pdb_surfaces('Zdock', 'complex.'+str(i)+'.pdb', 5.0)
        
        print "surface residues: ", len(residues)
        for usr in sorted(residues, key=lambda x: x[1]):
            print usr
            
            
def parse_hex(input_radius=5.0):
    for i in range(1,20+1):
        print "testing protein: ", i
        residues = parse_pdb_surfaces('Hex', 'hex-result'+("%04d" % i)+'.pdb', 5.0)
        
        print "surface residues: ", len(residues)
        for usr in sorted(residues, key=lambda x: x[1]):
            print usr
            
            
def parse_cluspro(input_radius=5.0):
    for i in range(1,20+1):
        print "testing protein: ", i
        residues = parse_pdb_surfaces('ClusPro', 'model.000.'+("%02d" % i)+'.pdb', 5.0)
        
        print "surface residues: ", len(residues)
        for usr in sorted(residues, key=lambda x: x[1]):
            print usr


def main():
    
    parser = argparse.ArgumentParser(description="Runner for batch_compare_script")
    t_group = parser.add_mutually_exclusive_group()
    t_group.add_argument('-z','--zdock',action='store_true')
    t_group.add_argument('-x','--hex',action='store_true')
    t_group.add_argument('-c','--cluspro',action='store_true')
    
    args = parser.parse_args()
    
    if args.zdock:
        parse_zdock()
    if args.hex:
        parse_hex()
    if args.cluspro:
        parse_cluspro()

if __name__ == '__main__':
	main()

