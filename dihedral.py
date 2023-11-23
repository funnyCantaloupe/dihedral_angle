from Bio.PDB import *
import numpy as np

#Sciezka do pliku, wczytanie pliku
pdb_file = "C:/Users/karol/PycharmProjects/Katy/1ehz.pdb"
parser = PDBParser()
structure = parser.get_structure("struktura", pdb_file)

#Definicja atomow dla katow torsyjnych - slownik, lista katow
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

torsion_angles = []

#Przejscie po resztach i lancuchach w strukturze
for model in structure:
    for chain in model:
        residues = list(chain.get_residues())
        for i in range(len(residues)):
            residue = residues[i]
            residue_number = residue.get_id()[1]

            torsion_values = [f"{residue_number}"]

            for angle, atoms in torsion_atoms.items():
                try:
                    atom1, atom2, atom3, atom4 = atoms

                    #Dla alpha atom1 powinien byc z poprzedniej reszty, pierwszy alpha jest None
                    if angle == 'alpha' and i == 0:
                        torsion_values.append(None)
                        continue

                    if angle == 'alpha' and i > 0:
                        atom1 = residues[i-1][atom1]
                    else:
                        atom1 = residue[atom1]

                    atom2 = residue[atom2]

                    #Dla zeta atom3 powinien byc z nastepnej reszty
                    if angle == 'zeta' and i < len(residues) - 1:
                        atom3 = residues[i + 1][atom3]
                    else:
                        atom3 = residue[atom3]

                    #Dla epsilon i zeta atom4 powinien byc z nastepnej reszty
                    if angle in ['epsilon', 'zeta'] and i < len(residues) - 1:
                        atom4 = residues[i + 1][atom4]
                    else:
                        atom4 = residue[atom4]

                    #Obliczanie kata torsyjnego
                    torsion_angle = calc_dihedral(atom1.get_vector(), atom2.get_vector(),
                                                  atom3.get_vector(), atom4.get_vector())

                    torsion_values.append(np.degrees(torsion_angle))

                except KeyError:
                    torsion_values.append(None)
            torsion_angles.append(torsion_values)

            #Kat chi, jesli puryna to nie pirydymina
            for values in torsion_values:
                if torsion_values[7] is not None:
                    torsion_values[8] = None

#Zapisanie wynikow
with open('katy_torsyjne.csv', 'w') as file:
    header = ",".join([angle for angle in torsion_atoms.keys()])
    file.write(f"Reszta,{header}\n")

    for angle_row in torsion_angles:
        line = ",".join([f"{value:.1f}" if isinstance(value, float) else str(value) for value in angle_row])
        file.write(f"{line}\n")


pdb_file2 = "C:/Users/karol/PycharmProjects/Katy/katy_torsyjne.csv"
with open(pdb_file2, 'r+') as file:
    lines = file.readlines()
    file.seek(0)

    for line in lines[:77]:
        file.write(line)

    file.truncate()

