import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def generate_points_on_sphere_golden_angle(radius, num_points):
    # Méthode avec golden angle pour une répartition uniforme
    golden_angle = np.pi * (3 - np.sqrt(5))
    theta = golden_angle * np.arange(num_points)
    z = np.linspace(1 - 1.0 / num_points, 1.0 / num_points - 1, num_points)
    radius_scaled = radius * np.sqrt(1 - z * z)
    x = radius_scaled * np.cos(theta)
    y = radius_scaled * np.sin(theta)
    return np.column_stack((x, y, z))

points_uniform = generate_points_on_sphere_golden_angle(radius, num_points)

