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

# Helper functions for formatting tick labels
def H_formatter(x, pos):
    return f'H{x:.0f}'

def t_formatter(x, pos):
    return f'{x:.0f}s'
    
# Sort the filtration
f.sort()

def plot_barcode(diagrams):
# Create a single figure and axes for the plot
    fig, ax = plt.subplots()

    # Calculate the number of diagrams to set y-ticks accordingly
    num_diagrams = len(diagrams)
    y_ticks = np.arange(num_diagrams)
    y_tick_labels = [f'$H_{{{i}}}$' for i in range(num_diagrams)]

    base_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    colors = base_colors * (len(diagrams) // len(base_colors) + 1)

    # Iterate through each dimension's diagram
    for i, dgm in enumerate(diagrams):
        # Plot each bar in the barcode
        for p in dgm:
            ax.plot([p.birth, p.death], [i, i], colors[i])

    # Setting custom formatter for the y-axis
    formatter = FuncFormatter(H_formatter)
    ax.yaxis.set_major_formatter(formatter)

    # Setting custom formatter for the x-axis
    formatter_t = FuncFormatter(t_formatter)
    ax.xaxis.set_major_formatter(formatter_t)

    # Set the y-ticks and labels according to homology dimensions
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels)

    # Positions and labels for the x-axis
    tick_positions = [1, 2, 3, 4] # Assuming these positions are relevant for your data
    tick_labels = ['$t_1$', '$t_2$', '$t_3$', '$t_4$'] # Labels for the x-axis
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    # Setting the plot title and labels
    plt.title("Persistence Barcode")
    plt.xlabel("Filtration value")
    plt.ylabel("Homology dimension")

    # Show the complete plot
    plt.show()
    plt.savefig('persistence_barcode.png')



# Call the plotting function with the calculated diagrams
plot_barcode(dgms)

