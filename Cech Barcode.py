import dionysus as d
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay

# Function to compute the radius of simplices for Delaunay triangulation
def compute_radius(points, simplex):
    A = points[simplex]
    # Note: 'cov' is calculated but not used; consider removing if unnecessary
    cov = np.cov(A.T)
    barycenter = np.mean(A, axis=0)
    radius = np.sqrt(np.max(np.sum((A - barycenter)**2, axis=1)))
    return radius

# Create a random point cloud
points = np.random.random((200, 3))

# Delaunay triangulation
tri = Delaunay(points)

# Initialize an empty filtration
f = d.Filtration()

# Track added simplices to avoid duplicates (especially important for edges)
added_simplices = set()

# Add vertices to the filtration
for i in range(len(points)):
    f.append(d.Simplex([i], 0))
    added_simplices.add(tuple([i]))

# Function to add all faces of a simplex
def add_faces(simplex, simplex_radius):
    # For a 3-simplex (tetrahedron) in 3D, add all faces (triangles)
    for i in range(4):
        face = tuple(sorted(set(simplex) - {simplex[i]}))  # Remove one vertex to get a face
        if face not in added_simplices:
            f.append(d.Simplex(face, simplex_radius))
            added_simplices.add(face)

# Add edges, faces, and tetrahedra from the Delaunay triangulation
for simplex in tri.simplices:
    simplex_radius = compute_radius(points, simplex)
    simplex_tuple = tuple(sorted(simplex))
    
    # Add faces of the simplex first
    add_faces(simplex, simplex_radius)
    
    # Then add the simplex itself (tetrahedron)
    if simplex_tuple not in added_simplices:
        f.append(d.Simplex(simplex, simplex_radius))
        added_simplices.add(simplex_tuple)
    
    # Iterate over edges of the simplex
    for i in range(4):
        for j in range(i + 1, 4):
            edge = tuple(sorted([simplex[i], simplex[j]]))
            if edge not in added_simplices:
                f.append(d.Simplex(edge, simplex_radius))
                added_simplices.add(edge)


# Sort the filtration
f.sort()

# Compute persistent homology, passing the filtration as an argument
p = d.homology_persistence(f)
dgms = d.init_diagrams(p, f)

# Function to plot the persistence barcode
def plot_barcode(diagrams):
    plt.figure(figsize=(10, 5))  # Adjust figure size for better visibility
    for i, dgm in enumerate(diagrams):
        for p in dgm:
            plt.plot([p.birth, p.death], [i, i], 'k')
    plt.title("Persistence Barcode")
    plt.xlabel("Time")
    plt.ylabel("Homology dimension")
    plt.show()

# Call the plotting function with the calculated diagrams
plot_barcode(dgms)

