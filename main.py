from Bio.PDB import PDBParser
import numpy as np
from Bio import PDB

def pobierz_strukture_z_pdb(pdb_id):
    pdbl = PDB.PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

    if pdb_file_path:
        print(f'Pobrano strukturę o ID {pdb_id} i zapisano w pliku {pdb_file_path}')
        return pdb_file_path
    else:
        print(f'Nie udało się pobrać struktury o ID {pdb_id}')
        return None


def calculate_dihedral_angle(atom1, atom2, atom3, atom4):
    #Oblicza kąt dihedryczny alpha
    if any(x is None for x in [atom1, atom2, atom3, atom4]):
        return None

    b1 = atom2 - atom1
    b2 = atom3 - atom2
    b3 = atom4 - atom3

    v1 = np.cross(b1, b2)
    v1 = v1 / np.linalg.norm(v1)

    v2 = np.cross(b2, b3)
    v2 = v2 / np.linalg.norm(v2)

    m1 = np.cross(v1, v2)
    m2 = np.cross(v1 / np.linalg.norm(v1), b2 / np.linalg.norm(b2))

    x = np.dot(m1, m2)
    y = np.dot(np.cross(v1, v2), m2)

    return np.degrees(np.arctan2(y, x))


def extract_coordinates(structure, atom_name, residue_number, chain_id='A'):
    #Pobiera współrzędne atomu w danej strukturze
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[1] == residue_number:
                        for atom in residue:
                            if atom.name == atom_name:
                                return atom.get_coord()
    return None


if __name__ == "__main__":
    pdb_id = '1FCW'
    pdb_file_path = pobierz_strukture_z_pdb(pdb_id)

    if pdb_file_path:
        # Wczytanie struktury z pliku PDB
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_file_path)

        # Obliczenia dla kąta alpha dla każdej reszty
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == ' ' and residue.id[2] == ' ':
                        residue_number = residue.id[1]
                        n = residue_number
                        alpha_atoms = [
                            extract_coordinates(structure, 'O3\'', n - 1),
                            extract_coordinates(structure, 'P', n - 1),
                            extract_coordinates(structure, 'O5\'', n),
                            extract_coordinates(structure, 'C5\'', n)
                        ]
                        alpha_angle = calculate_dihedral_angle(*alpha_atoms)

                        if alpha_angle is not None:
                            print(f'Reszta {residue_number}: Kąt dihedryczny alpha: {alpha_angle} stopni')
                        else:
                            print(f'Nie można obliczyć kąta dla reszty {residue_number}')
