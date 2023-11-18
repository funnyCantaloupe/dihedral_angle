from Bio.PDB import *
import matplotlib.pyplot as plt

# Ścieżka do pliku PDB
pdb_file = "C:/Users/karol/PycharmProjects/Katy/1mbo.pdb"  # Zmień na właściwą ścieżkę

# Inicjalizacja parsera PDB
parser = PDBParser()

# Wczytanie struktury z pliku PDB
structure = parser.get_structure("struktura", pdb_file)

# Lista kątów phi i psi
phi_values = []
psi_values = []
structure_colors = []  # Lista przechowująca kolory dla poszczególnych struktur drugorzędowych

# Iteracja po łańcuchach i resztach w strukturze
for model in structure:
    for chain in model:
        polypeptides = PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides):
            phi_psi = poly.get_phi_psi_list()
            for phi, psi in phi_psi:
                if phi is not None and psi is not None:
                    # Konwersja kątów z radianów na stopnie
                    phi_deg = phi * (180.0 / 3.14159)
                    psi_deg = psi * (180.0 / 3.14159)
                    phi_values.append(phi_deg)
                    psi_values.append(psi_deg)

                    # Przypisanie koloru na podstawie kryteriów dla helisy alfa i arkusza beta
                    if phi_deg < 0 and psi_deg < 0:
                        structure_colors.append('blue')  # Helisa alfa - niebieski kolor
                    elif phi_deg < 0 and psi_deg > 0:
                        structure_colors.append('red')  # Arkusz beta - czerwony kolor
                    else:
                        structure_colors.append('green')  # Pozostałe - zielony kolor

# Stworzenie wykresu Ramachandrana z kolorami dla struktur drugorzędowych
plt.figure(figsize=(8, 6))
for color, phi, psi in zip(structure_colors, phi_values, psi_values):
    plt.scatter(phi, psi, s=5, alpha=0.5, c=color)

plt.title('Mapa Ramachandrana z kolorem dla struktur drugorzędowych')
plt.xlabel('Kąt Phi (°)')
plt.ylabel('Kąt Psi (°)')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.grid(True)

# Legenda dla typów struktury drugorzędowej
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label='Helisa alfa',
                              markerfacecolor='blue', markersize=8),
                   plt.Line2D([0], [0], marker='o', color='w', label='Arkusze beta',
                              markerfacecolor='red', markersize=8),
                   plt.Line2D([0], [0], marker='o', color='w', label='Inne',
                              markerfacecolor='green', markersize=8)]
plt.legend(handles=legend_elements)

plt.show()
