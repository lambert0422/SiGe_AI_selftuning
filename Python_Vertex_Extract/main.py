import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, signal
# # ===== Parameters for figure2 =====
# VERTEX_DENSITY = 0.001  # Control vertex count (0.001=more vertices, 0.05=fewer)
# MIN_AREA = 100  # Minimum contour area to consider
# BLUR_SIZE = 5  # Gaussian blur kernel size
# THRESHOLD = 150  # Binary threshold value
# ===== Parameters for figure1 =====
# VERTEX_DENSITY = 0.004  # Control vertex count (0.001=more vertices, 0.05=fewer)
# MIN_AREA = 100  # Minimum contour area to consider
# BLUR_SIZE = 5  # Gaussian blur kernel size
# THRESHOLD = 200  # Binary threshold value
# scale = 400/215 # Correct the physical size to figure resolution
# image = cv2.imread('figure.png')
# gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# ===== Parameters for figure3 =====
VERTEX_DENSITY = 0.0001  # Control vertex count (0.001=more vertices, 0.05=fewer)
MIN_AREA = 50  # Minimum contour area to consider
BLUR_SIZE = 5  # Gaussian blur kernel size
THRESHOLD = 70  # Binary threshold value
scale = 3 # Correct the physical size to figure resolution
image = cv2.imread('figure3.png')
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
gray = 255 - gray  # Invert if needed


##################################################
# num_levels = 4
#
# # Calculate the bin size (256 / 4 = 64)
# bin_size = 256 // num_levels
#
# quantized = (gray // bin_size) * bin_size  # Alternatively: np.floor(gray / bin_size) * bin_size
#
# # quantized = np.floor(gray / 64) * 85  # Adjust as needed
#
# # Convert back to uint8 (if needed)
# quantized = np.uint8(quantized)
#
# # Display results
# cv2.imshow('Original', image)
# cv2.imshow('Grayscale', gray)
# cv2.imshow('Quantized (4 levels)', quantized)
# cv2.waitKey(0)
# cv2.destroyAllWindows()

##################################################
# Show processing steps
plt.figure(figsize=(15, 5))
plt.subplot(131), plt.imshow(gray, cmap='gray'), plt.title("1. Original")
blurred = cv2.GaussianBlur(gray, (BLUR_SIZE, BLUR_SIZE), 0)
plt.subplot(132), plt.imshow(blurred, cmap='gray'), plt.title("2. Blurred")
_, binary = cv2.threshold(blurred, THRESHOLD, 255, cv2.THRESH_BINARY_INV)
plt.subplot(133), plt.imshow(binary, cmap='gray'), plt.title("3. Binary")
plt.tight_layout()
plt.show()

contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
filtered_contours = [cnt for cnt in contours if cv2.contourArea(cnt) > MIN_AREA]
print(f"Found {len(filtered_contours)} patterns")


def generate_grid_with_density(vertices, axis='x', base_density=1.0, max_numPoints = 5,min_spacing = 2):
    coords = sorted({v[0] if axis == 'x' else v[1] for v in vertices})

    if len(coords) < 2:
        return coords

    distances = np.diff(coords)
    neighbor_dists = []
    for i in range(len(coords)):
        if i == 0:
            neighbor_dists.append(distances[0])
        elif i == len(coords) - 1:
            neighbor_dists.append(distances[-1])
        else:
            neighbor_dists.append(min(distances[i - 1], distances[i]))


    grid_positions = [coords[0]]  # Always include first point

    for i in range(len(coords) - 1):
        start = coords[i]
        end = coords[i + 1]
        segment_length = end - start


        local_density = base_density * (1 - neighbor_dists[i] / max(neighbor_dists))



        num_additional = min(max_numPoints,
                             int((segment_length / min_spacing)))


        for j in range(1, num_additional + 1):
            position = start + j * (segment_length / (num_additional ))
            # grid_positions.append(round(position, 1))
            grid_positions.append(position)

        grid_positions.append(end)  # Add endpoint

    return sorted(list(set(grid_positions)))


def generate_grid_definitions(vertices, axis='x', var_prefix='GRID_x', density=1.0):

    positions = generate_grid_with_density(vertices, axis, density)

    if len(positions) < 2:
        return []

    definitions = []
    current_start = positions[0]
    current_spacing = None
    prev_spacing = None
    spacing_count = 0


    definitions.append(
        f"#IF $M3D         line{{ pos = ${axis}_min             spacing = ${var_prefix}_Junc }}"
    )

    # Process middle points
    for i in range(1, len(positions) - 1):
        spacing = round(min(positions[i] - positions[i - 1], positions[i + 1] - positions[i]), 1)

        if spacing == prev_spacing:
            spacing_count += 1
            current_end = positions[i]


            if spacing_count >= 2:
                continue
        else:

            if prev_spacing is not None and spacing_count >= 1:
                definitions.append(
                    f"#IF $M3D         line{{ pos = {current_start:<10} spacing = {prev_spacing}}}"
                )
                if spacing_count >= 2:
                    definitions.append(
                        f"#IF $M3D         line{{ pos = {current_end:<10} spacing = {prev_spacing}}}"
                    )


            current_start = positions[i]
            prev_spacing = spacing
            spacing_count = 1


    if prev_spacing is not None and spacing_count >= 1:
        definitions.append(
            f"#IF $M3D         line{{ pos = {current_start:<10} spacing = {prev_spacing}}}"
        )
        if spacing_count >= 2:
            definitions.append(
                f"#IF $M3D         line{{ pos = {positions[-2]:<10} spacing = {prev_spacing}}}"
            )


    definitions.append(
        f"#IF $M3D         line{{ pos = ${axis}_max             spacing = ${var_prefix}_Junc }}"
    )

    return definitions
def min_distance(coords, axis):

    axis_coords = sorted({v[0] if axis == 'x' else v[1] for v in coords})
    if len(axis_coords) < 2:
        return 0
    return min(axis_coords[i + 1] - axis_coords[i] for i in range(len(axis_coords) - 1))


def get_vertices(contour, density=VERTEX_DENSITY):
    """Returns vertices with controllable density"""
    perimeter = cv2.arcLength(contour, True)
    epsilon = density * perimeter  # Key control parameter
    return cv2.approxPolyDP(contour, epsilon, True).reshape(-1, 2)



contour_img = image.copy()
vertex_img = image.copy()
patterns = []
for i, cnt in enumerate(filtered_contours):

    vertices = get_vertices(cnt)


    cv2.drawContours(contour_img, [cnt], -1, (0, 255, 0), 2)

    M = cv2.moments(cnt)
    if M["m00"] != 0:
        cX = int(M["m10"] / M["m00"])
        cY = int(M["m01"] / M["m00"])
    else:

        cX, cY = cnt[0][0]


    label = str(i+1)
    font = cv2.FONT_HERSHEY_SIMPLEX
    font_scale = 0.7
    font_thickness = 2
    text_color = (255, 255, 255)
    bg_color = (0, 0, 0)


    (text_width, text_height), _ = cv2.getTextSize(label, font, font_scale, font_thickness)


    cv2.rectangle(contour_img,
                  (cX - 5, cY - text_height - 5),
                  (cX + text_width + 5, cY + 5),
                  bg_color, -1)


    cv2.putText(contour_img, label,
                (cX, cY),
                font, font_scale, text_color, font_thickness)

    for (x, y) in vertices:
        cv2.circle(vertex_img, (x, y), 5, (0, 0, 255), -1)
        # cv2.putText(vertex_img, f"({x},{y})", (x + 10, y - 5),
        #             cv2.FONT_HERSHEY_SIMPLEX, 0.4, (255, 255, 0), 1)


    print(f"\nPattern {i + 1} - {len(vertices)} vertices:")
    print((vertices * scale).astype(int))
    patterns.append((vertices * scale).astype(int))

patterns = np.vstack(patterns)
density = 1
x_definitions = generate_grid_definitions(
    vertices=patterns,
    axis='x',
    var_prefix='GRID_x',
    density=density
)


y_definitions = generate_grid_definitions(
    vertices=patterns,
    axis='y',
    var_prefix='GRID_y',
    density=density
)

plt.figure(figsize=(15, 6))
plt.subplot(121), plt.imshow(cv2.cvtColor(contour_img, cv2.COLOR_BGR2RGB))
plt.title("Contours")
plt.subplot(122), plt.imshow(cv2.cvtColor(vertex_img, cv2.COLOR_BGR2RGB))
plt.title(f"Vertices (density={VERTEX_DENSITY})")
plt.tight_layout()
plt.show()


with open("vertices_output.txt", 'w') as f:
    for i, cnt in enumerate(filtered_contours):
        vertices = get_vertices(cnt)
        f.write(f"\n# Pattern {i + 1}\n")
        for (x, y) in (vertices * scale).astype(int):
            f.write(f"#IF $M3D vertex{{ z = {x} y = {y} }}\n")
print("Vertices saved to vertices_output.txt")

with open("grid_definitions.txt", 'w') as f:
    for line in x_definitions:
        f.write(line + "\n")  # Add newline after each definition
    for line in y_definitions:
        f.write(line + "\n")  # Add newline after each definition
print(f"Grid definitions successfully written to grid_definitions.txt")