import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def generate_points_on_sphere_random(radius, num_points):
    # Méthode avec des angles aléatoires
    theta = 2 * np.pi * np.random.rand(num_points)
    phi = np.arccos(2 * np.random.rand(num_points) - 1)
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)
    return np.column_stack((x, y, z))

radius = 1.0  # Rayon de la sphère
num_points = 100  # Nombre de points à générer

points_random = generate_points_on_sphere_random(radius, num_points)

# Créez une figure 3D pour afficher les points
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Tracez les points générés avec la méthode aléatoire en rouge
ax.scatter(points_random[:, 0], points_random[:, 1], points_random[:, 2], s=10, c='r', marker='x', label='Aléatoire')

# Paramètres d'axe
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Titre
ax.set_title('Points générés sur une sphère')

# Légende
ax.legend()

plt.show()