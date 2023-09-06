import numpy as np
import requests
from io import StringIO  # Importez StringIO depuis le module io
from Bio import PDB
import pandas as pd


# création d'un dictionnaire avec les rayons de Van der Waals moyens des 10 atomes les plus fréquemment retrouvés
vdw_radii = {
    'H': 1.20,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'S': 1.80,
    'P': 1.80,
    'MG': 1.50,
    'FE': 1.80,
    'ZN': 1.39,
    'CA': 2.00,
}


# Function to generate points on a sphere uniformly
def generate_points_on_sphere(radius, num_points):
    # Méthode avec golden angle pour une répartition uniforme
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(num_points)
    z = np.linspace(1 - 1.0 / num_points, 1.0 / num_points - 1, num_points)
    radius_scaled = radius * np.sqrt(1 - z * z)
    x = radius_scaled * np.cos(theta)
    y = radius_scaled * np.sin(theta)
    return np.column_stack((x, y, z))

# Function to fetch the data from the RCSB web service
def fetch_pdb_data(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_data = response.text
        return pdb_data
    else:
        raise Exception(f"Impossible de télécharger les données PDB pour {pdb_id}")

def atomic_coordinates(pdb_data):
    atom_data = {"res_name": [], "res_num": [], "atom_name": [], "x": [], "y": [], "z": []}
    atoms_names, res_names, res_nums, x, y, z = [], [], [], [], [], []
    
    current_model = None  # Gardez une trace du numéro de modèle actuel
    
    for line in pdb_data.splitlines():
        if line.startswith("MODEL"):
            current_model = int(line[10:14])  # Extrait le numéro de modèle
        elif line.startswith("ATOM") and (current_model is None or current_model == 1):
            atom_name = line[13:16].strip()
            if atom_name != "H":
                atoms_names.append(line[13].strip())
                res_names.append(line[17:20].strip())
                res_nums.append(line[22:26])
                x.append(float(line[30:38]))
                y.append(float(line[38:46]))
                z.append(float(line[46:54]))
    
    atom_data["res_name"] = res_names
    atom_data["res_num"] = res_nums
    atom_data["atom_name"] = atoms_names
    atom_data["x"], atom_data["y"], atom_data["z"] = x, y, z
    
    df = pd.DataFrame(atom_data)
    print(df)
    return df

# function to calculate surface area 
def calculate_surface_area(pdb_id, probe_radius):
    pdb_data = fetch_pdb_data(pdb_id)
    print(pdb_data)
    # Analyse de la structure PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', StringIO(pdb_data))

    atom_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_coords.append(atom.get_coord())
                    if atom.element != "H":  # exclure les atomes d'hydrogène

   
                    # Calcul des points accessibles et de la surface accessible
                        accessible_points = []
                        num_points = 100  # Nombre de points sur la sphère
                        points_on_sphere = generate_points_on_sphere(probe_radius, num_points)
                        atom_areas = []

                        for atom_coord in atom_coords:
                                sphere = points_on_sphere + atom_coord
                                for point in sphere:
                                    is_accessible = True
                                    for other_atom_coord in atom_coords:
                                        if np.linalg.norm(point - other_atom_coord) - vdw_radii[atom.element] < probe_radius:
                                            is_accessible = False
                                            break
                                    if is_accessible:
                                        accessible_points.append(point)
                                        print(f"Point {point} est accessible.")

                                atom_accessible_area = len(accessible_points) * (4 * np.pi * (probe_radius)**2) / num_points
                                atom_areas.append(atom_accessible_area)

                        accessible_surface_area = sum(atom_areas)
                        return accessible_surface_area     

pdb_id = "1CRN"  # Remplacez par le PDB ID de votre protéine
probe_radius = 1.4  # Rayon de la sonde (par exemple, pour l'oxygène)

accessible_area = calculate_surface_area(pdb_id, probe_radius)

print(f"Surface accessible absolue : {accessible_area:.2f} Å²")

