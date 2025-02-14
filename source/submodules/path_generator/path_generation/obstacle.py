import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle

@dataclass
class Obstacle:
    center: np.ndarray
    radius: np.double
    height: np.double = 0

def plot_2D_obstacle(obstacle: Obstacle, ax):
    circle = plt.Circle((obstacle.center.item(0), obstacle.center.item(1)), 
                        obstacle.radius, color='r')
    ax.add_patch(circle)

def plot_2D_obstacles(obstacles: list, ax):
    for i in range(len(obstacles)):
        plot_2D_obstacle(obstacles[i], ax)

def plot_3D_obstacle(obstacle: Obstacle, ax):
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = obstacle.radius * np.cos(u)*np.sin(v) + obstacle.center.item(0)
    y = obstacle.radius * np.sin(u)*np.sin(v) + obstacle.center.item(1)
    z = obstacle.radius * np.cos(v) + obstacle.center.item(2)
    ax.plot_surface(x, y, z, color="r")

def plot_3D_obstacles(obstacles: list, ax):
    for i in range(len(obstacles)):
        plot_3D_obstacle(obstacles[i], ax)

def plot_cylinder(obstacle: Obstacle, ax):
    z = np.linspace(0, -obstacle.height, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = obstacle.radius*np.cos(theta_grid) + obstacle.center.item(0)
    y_grid = obstacle.radius*np.sin(theta_grid) + obstacle.center.item(1)
    ax.plot_surface(x_grid,y_grid,z_grid, color="r")

def plot_cylinders(obstacles: list, ax):
    for i in range(len(obstacles)):
        plot_cylinder(obstacles[i], ax)



def plot_buildings(obstacles: list, ax):
    for i in range(len(obstacles)):
        obstacle = obstacles[i]
        center = obstacle.center
        plot_rectangular_prism(ax, center.item(0), center.item(1), -obstacle.height/2, obstacle.radius, obstacle.radius, -obstacle.height)

def plot_rectangular_prism(ax, north, east, down, length, width, height):

    # Calculate half-lengths in each dimension
    half_length = length / 2
    half_width = width / 2
    half_height = height / 2

    # Define vertices of the rectangular prism
    vertices = np.array([
        [north - half_length, east - half_width, down - half_height],
        [north + half_length, east - half_width, down - half_height],
        [north + half_length, east + half_width, down - half_height],
        [north - half_length, east + half_width, down - half_height],
        [north - half_length, east - half_width, down + half_height],
        [north + half_length, east - half_width, down + half_height],
        [north + half_length, east + half_width, down + half_height],
        [north - half_length, east + half_width, down + half_height]
    ])

    # Define faces of the rectangular prism
    faces = [
        [vertices[0], vertices[1], vertices[2], vertices[3]],
        [vertices[4], vertices[5], vertices[6], vertices[7]],
        [vertices[0], vertices[1], vertices[5], vertices[4]],
        [vertices[2], vertices[3], vertices[7], vertices[6]],
        [vertices[0], vertices[3], vertices[7], vertices[4]],
        [vertices[1], vertices[2], vertices[6], vertices[5]]
    ]

    # Plot faces
    ax.add_collection3d(Poly3DCollection(faces, facecolors='lime', linewidths=1, edgecolors='g', alpha=0.5))


def plot_2d_buildings(obstacles: list, ax):
    for i in range(len(obstacles)):
        obstacle = obstacles[i]
        center = obstacle.center
        plot_square(ax, center.item(0), center.item(1), obstacle.radius)

def plot_square(ax, north, east, width):
    # Calculate the coordinates of the square's corners
    x = north - width / 2
    y = east - width / 2
    
    # Create a rectangle patch
    rect = Rectangle((x, y), width, width, linewidth=1, edgecolor='r', facecolor='none')
    
    # Add the rectangle to the plot
    ax.add_patch(rect)