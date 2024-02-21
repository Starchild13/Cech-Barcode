import dionysus as d
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay

# Function to compute the radius of simplices for Delaunay triangulation
def compute_radius(points, simplex):
    A = points[simplex]
    cov = np.cov(A.T)
    barycenter = np.mean(A, axis=0)
    radius = np.sqrt(np.max(np.sum((A - barycenter)**2, axis=1)))
    return radius

# Create a random point cloud
points = np.random.random((200, 2))

# Delaunay triangulation
tri = Delaunay(points)

# Initialize a filtration
f = d.Filtration()

# Add vertices to the filtration
for i in range(len(points)):
    f.append(d.Simplex([i], 0))

# Add edges and triangles from the Delaunay triangulation
for simplex in tri.simplices:
    simplex_radius = compute_radius(points, simplex)  # Use compute_radius to avoid name conflict
    f.append(d.Simplex(list(simplex), simplex_radius))
    for edge in [[simplex[0], simplex[1]], [simplex[0], simplex[2]], [simplex[1], simplex[2]]]:
        f.append(d.Simplex(edge, simplex_radius))

# Sort the filtration
f.sort()

# Compute persistent homology
p = d.homology_persistence(f)
dgms = d.init_diagrams(p, f)

# Function to plot the persistence barcode
def plot_barcode(diagrams):
    for i, dgm in enumerate(diagrams):
        for p in dgm:
            plt.plot([p.birth, p.death], [i, i], 'k')
    plt.title("Persistence Barcode")
    plt.xlabel("Time")
    plt.ylabel("Homology dimension")
    plt.show()

# Call the plotting function with the calculated diagrams
plot_barcode(dgms)
