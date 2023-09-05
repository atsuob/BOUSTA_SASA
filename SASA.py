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

def generate_points_on_sphere_golden_angle(radius, num_points):
    # Méthode avec golden angle pour une répartition uniforme
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(num_points)
    z = np.linspace(1 - 1.0 / num_points, 1.0 / num_points - 1, num_points)
    radius_scaled = radius * np.sqrt(1 - z * z)
    x = radius_scaled * np.cos(theta)
    y = radius_scaled * np.sin(theta)
    return np.column_stack((x, y, z))

radius = 1.0  # Rayon de la sphère
num_points = 100  # Nombre de points à générer

points_random = generate_points_on_sphere_random(radius, num_points)
points_uniform = generate_points_on_sphere_golden_angle(radius, num_points)

# Créez une figure 3D pour afficher les points
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Tracez les points générés avec la méthode aléatoire en rouge
ax.scatter(points_random[:, 0], points_random[:, 1], points_random[:, 2], s=10, c='r', marker='x', label='Aléatoire')
# Tracez les points générés avec la méthode golden angle en bleu
ax.scatter(points_uniform[:, 0], points_uniform[:, 1], points_uniform[:, 2], s=10, c='b', marker='o', label='Golden Angle')


# Paramètres d'axe
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Titre
ax.set_title('Points générés sur une sphère')

# Légende
ax.legend()

plt.show()