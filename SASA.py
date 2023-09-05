import numpy as np
import requests

# Function to generate points on a sphere uniformly
def generate_points_on_sphere_golden_angle(radius, num_points):
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
        print(pdb_data)
        return pdb_data
    else:
        raise Exception(f"Impossible de télécharger les données PDB pour {pdb_id}")

pdb_id = '1C75'
fetch_pdb_data = fetch_pdb_data(pdb_id)