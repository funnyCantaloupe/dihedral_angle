from Bio.PDB import *
import numpy as np

# Ścieżka do pliku PDB
pdb_file = "C:/Users/karol/PycharmProjects/Katy/1ehz.pdb"  # Zmień na właściwą ścieżkę

# Inicjalizacja parsera PDB
parser = PDBParser()

# Wczytanie struktury z pliku PDB
structure = parser.get_structure("struktura", pdb_file)

# Definicja atomów dla kątów torsyjnych
torsion_atoms = {
    'alpha': ('O3\'', 'P', 'O5\'', 'C5\''),
    'beta': ('P', 'O5\'', 'C5\'', 'C4\''),
    'gamma': ('O5\'', 'C5\'', 'C4\'', 'C3\''),
    'delta': ('C5\'', 'C4\'', 'C3\'', 'O3\''),
    'epsilon': ('C4\'', 'C3\'', 'O3\'', 'P'),
    'zeta': ('C3\'', 'O3\'', 'P', 'O5\''),
    'chi_purines': ('O4\'', 'C1\'', 'N9', 'C4'),
    'chi_pyrimidines': ('O4\'', 'C1\'', 'N1', 'C2')
}

# Lista kątów torsyjnych
torsion_angles = []

# Iteracja po łańcuchach i resztach w strukturze
for model in structure:
    for chain in model:
        residues = chain.get_residues()
        for residue in residues:
            residue_number = residue.get_id()[1]  # Numer reszty

            torsion_values = [f"{residue_number}"]  # Numer reszty jako pierwsza wartość

            for angle, atoms in torsion_atoms.items():
                try:
                    atom1, atom2, atom3, atom4 = atoms
                    # Pobieranie atomów dla kąta torsyjnego
                    atom1 = residue[atom1]
                    atom2 = residue[atom2]
                    atom3 = residue[atom3]
                    atom4 = residue[atom4]

                    # Obliczanie kąta torsyjnego
                    torsion_angle = calc_dihedral(atom1.get_vector(), atom2.get_vector(),
                                                  atom3.get_vector(), atom4.get_vector())

                    torsion_values.append(np.degrees(torsion_angle))
                except KeyError:
                    torsion_values.append(None)
            torsion_angles.append(torsion_values)

# Zapisanie wyników do pliku CSV jako macierz n x m
with open('kąty_torsyjne.csv', 'w') as file:
    header = ",".join([angle for angle in torsion_atoms.keys()])
    file.write(f"Reszta,{header}\n")

    for angle_row in torsion_angles:
        line = ",".join([str(value) for value in angle_row])
        file.write(f"{line}\n")
