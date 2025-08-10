import numpy as np
from stl import mesh  # Install with: pip install numpy-stl


def calculate_offset_vertices(vertices, offset):
    """Calculates outward-offset vertices in XY plane."""
    n = len(vertices)
    offset_vertices = []

    for i in range(n):
        prev = vertices[i - 1]
        curr = vertices[i]
        next = vertices[(i + 1) % n]

        # Edge vectors
        edge1 = curr - prev
        edge2 = next - curr

        # Perpendicular normals (outward-facing)
        norm1 = np.array([-edge1[1], edge1[0]])
        norm2 = np.array([-edge2[1], edge2[0]])

        # Normalize
        norm1 = norm1 / np.linalg.norm(norm1)
        norm2 = norm2 / np.linalg.norm(norm2)

        # Average normal (smooth corners)
        avg_norm = (norm1 + norm2) / 2
        avg_norm = avg_norm / np.linalg.norm(avg_norm)

        # Offset vertex
        offset_vertex = curr + avg_norm * offset
        offset_vertices.append(offset_vertex)

    return np.array(offset_vertices)


def generate_sealed_membrane(original_vertices, membrane_thickness, total_height):
    """Generates a sealed 3D membrane covering the pattern."""
    outer_vertices = calculate_offset_vertices(original_vertices, membrane_thickness)
    n = len(original_vertices)

    # All vertices (order: inner-top, outer-top, inner-bottom, outer-bottom)
    vertices = []

    # 1. Inner vertices (top)
    for v in original_vertices:
        vertices.append([v[0], v[1], total_height])

    # 2. Outer vertices (top)
    for v in outer_vertices:
        vertices.append([v[0], v[1], total_height])

    # 3. Inner vertices (bottom)
    for v in original_vertices:
        vertices.append([v[0], v[1], total_height - membrane_thickness])

    # 4. Outer vertices (bottom)
    for v in outer_vertices:
        vertices.append([v[0], v[1], total_height - membrane_thickness])

    vertices = np.array(vertices)

    # Define triangles for the mesh
    faces = []

    # Helper to add a triangle
    def add_face(a, b, c):
        faces.append([a, b, c])

    # --- Top cover (sealed surface) ---
    for i in range(n):
        next_i = (i + 1) % n
        add_face(i, next_i, n + i)  # Inner to outer
        add_face(next_i, n + next_i, n + i)  # Complete quad

    # --- Bottom cover ---
    for i in range(n):
        next_i = (i + 1) % n
        add_face(2 * n + i, 2 * n + next_i, 3 * n + i)  # Inner to outer
        add_face(2 * n + next_i, 3 * n + next_i, 3 * n + i)  # Complete quad

    # --- Side walls ---
    for i in range(n):
        next_i = (i + 1) % n
        # Inner wall
        add_face(i, 2 * n + i, 2 * n + next_i)
        add_face(i, 2 * n + next_i, next_i)
        # Outer wall
        add_face(n + i, 3 * n + i, 3 * n + next_i)
        add_face(n + i, 3 * n + next_i, n + next_i)

    # Create the mesh
    mesh_data = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            mesh_data.vectors[i][j] = vertices[f[j]]

    return mesh_data, vertices


# Input vertices (your 2D pattern)
original_vertices = np.array([
    [844, 1038],
    [837, 1043],
    [837, 1233],
    [870, 1233],
    [870, 1041],
    [866, 1038]
])

# Parameters
membrane_thickness = 20  # 20-unit thick film
total_height = 50  # Original pattern height

# Generate the sealed membrane
mesh_data, all_vertices = generate_sealed_membrane(original_vertices, membrane_thickness, total_height)

# Save to STL (3D printable/viewable)
mesh_data.save('sealed_membrane.stl')

# Save vertices to TXT (for verification)
with open('sealed_vertices.txt', 'w') as f:
    current_z = None
    for v in all_vertices:
        if v[2] != current_z:
            current_z = v[2]
            f.write(f"z = [{current_z}]\n")
        f.write(f"vertex{{ x = {v[0]}, y = {v[1]}, z = {v[2]} }}\n")

print("Files saved:")
print("- 3D model: 'sealed_membrane.stl' (open in any 3D viewer)")
print("- Vertex list: 'sealed_vertices.txt'")