from Bio.PDB import *
import sys
import argparse

surface_residues = []

class PDBSelect(Select):
    def accept_residue(self, residue):
        test_residue = (residue.get_parent().get_id(),residue.get_resname(), residue.get_id()[1])
        if test_residue in surface_residues:
            return 1
        else:
            return 0

def parse_pdb_surfaces(folder='name', name='input.pdb', search_radius=5.0, out_pdb='pdb_out'):

    global surface_residues
    
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
            surface_residues += [('B', r.get_resname(), r.id[1])]
            surface_residues += [('A', atom.get_parent().get_resname(), atom.get_parent().id[1])]

    surface_residues = set(surface_residues)
    #print surface_residues
    
    io=PDBIO()
    io.set_structure(structure)
    io.save(folder + '/' + out_pdb, PDBSelect())


def parse_zdock(input_radius=5.0):
    for i in range(1,10+1):
        parse_pdb_surfaces('Zdock', 'complex.'+str(i)+'.pdb', 5.0, 'interface.'+("%02d" % i)+'.pdb')

def parse_hex(input_radius=5.0):
    for i in range(1,20+1):
        parse_pdb_surfaces('Hex', 'hex-result'+("%04d" % i)+'.pdb', 5.0, 'interface.'+("%02d" % i)+'.pdb')

def parse_cluspro(input_radius=5.0):
    for i in range(1,20+1):
        parse_pdb_surfaces('ClusPro', 'model.000.'+("%02d" % i)+'.pdb', 5.0, 'interface.'+("%02d" % i)+'.pdb')

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
