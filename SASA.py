import numpy as np
import requests
from io import StringIO  # Importez StringIO depuis le module io
from Bio import PDB
import pandas as pd
import math



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

def calculate_surface_area(i, atom1, df, probe_radius=1.4):

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
    'CA': 2.00,}
    name_atom1, coord_atom1, name_res, num_res = atom1["atom_name"], [atom1["x"], atom1["y"], atom1["z"]], atom1["res_name"], atom1["res_num"]
    r_vdw_atom1 = vdw_radii.get(name_atom1, 1.5)
    num_points = 10
    points_on_sphere = generate_points_on_sphere(r_vdw_atom1, num_points)
    points_on_sphere = points_on_sphere + coord_atom1
    unaccessible_points = []
    for j, atom2 in df.iterrows():
        if i == j: continue
        name_atom2, coord_atom2 = atom2["atom_name"], [atom2["x"], atom2["y"], atom2["z"]]
        r_vdw_atom2 = vdw_radii.get(name_atom2, 1.5)
        #distance = np.linalg.norm(np.array(coord_atom1) - np.array(coord_atom2))
        distance = math.dist(coord_atom1, coord_atom2)
        if distance < 5.0:
            for k, point in enumerate(points_on_sphere):
                if k not in unaccessible_points:
                    coord_point_atom1 = [point[0], point[1], point[2]]
                    #print ("distance entre atom1 et atom2 ", distance)
                    #print (coord_point_atom1)
                    #distance_point_atom2 = np.linalg.norm(np.array(coord_point_atom1) - np.array(coord_atom2))
                    distance_point_atom2 = math.dist(coord_point_atom1, coord_atom2)
                    print("la distance entre atom 1 et atom 2 :", distance_point_atom2 )
                    compare = 2 * probe_radius + r_vdw_atom2
                    print(compare)
                    if distance_point_atom2 < compare:
                        unaccessible_points.append(k)  

    print(len(unaccessible_points))
    ratio_exp = (num_points - len(unaccessible_points)) / num_points
    atom_accessible_area = (4 * np.pi * (r_vdw_atom1 + probe_radius) ** 2) * ratio_exp


    liste = [name_res, num_res, round(ratio_exp, 2), round(atom_accessible_area, 2)]
    return liste

def calculate_sasa(pdb_id):
    pdb_data = fetch_pdb_data(pdb_id)
    df = atomic_coordinates(pdb_data)
    print (df)
    results = []

    for i, atom1 in df.iterrows():
        sasa_data = calculate_surface_area(i, atom1, df)
        results.append(sasa_data)

    sasa_df = pd.DataFrame(results, columns=["Residue Name", "Residue Number", "Accessibility Ratio", "Atom Accessible Area"])

    # Calculez la somme totale des surfaces accessibles
    total_surface_area = sasa_df["Atom Accessible Area"].sum()
    
    return sasa_df, total_surface_area  # Renvoie le DataFrame et la surface totale