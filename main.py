from Bio.PDB import *
import matplotlib.pyplot as plt

#Sciezka do pliku, wczytanie pliku
pdb_file = "C:/Users/karol/PycharmProjects/Katy/101m.pdb"
parser = PDBParser()
structure = parser.get_structure("struktura", pdb_file)

#Lista katow phi i psi, lista kolorow dla struktur drugorzedowych
phi_values = []
psi_values = []
structure_colors = []

#Przejscie po lancuchach i resztach
for model in structure:
    for chain in model:
        polypeptides = PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides):
            phi_psi = poly.get_phi_psi_list()
            for phi, psi in phi_psi:
                if phi is not None and psi is not None:
                    #Radiany na stopnie
                    phi_deg = phi * (180.0 / 3.14159)
                    psi_deg = psi * (180.0 / 3.14159)
                    phi_values.append(phi_deg)
                    psi_values.append(psi_deg)

                    #Kolor dla alpha helix, beta lub inne
                    if phi_deg < 0 and psi_deg < 0:
                        structure_colors.append('blue')
                    elif phi_deg < 0 and psi_deg > 0:
                        structure_colors.append('red')
                    else:
                        structure_colors.append('green')

#Wykres Ramachandrana
plt.figure(figsize=(8, 6))
for color, phi, psi in zip(structure_colors, phi_values, psi_values):
    plt.scatter(phi, psi, s=5, alpha=0.5, c=color)

plt.title('Mapa Ramachandrana z kolorem dla struktur drugorzedowych')
plt.xlabel('Kat Phi (°)')
plt.ylabel('Kat Psi (°)')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.grid(True)

#Legenda
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label='Alpha Helix',
                              markerfacecolor='blue', markersize=8),
                   plt.Line2D([0], [0], marker='o', color='w', label='Beta Sheet',
                              markerfacecolor='red', markersize=8),
                   plt.Line2D([0], [0], marker='o', color='w', label='Inne',
                              markerfacecolor='green', markersize=8)]
plt.legend(handles=legend_elements)
plt.show()
