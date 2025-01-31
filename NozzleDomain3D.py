import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mesh_tool import *  
from scipy.interpolate import interp1d

R_I = 1.0       # Inlet radius
R_star = 0.5    # Throat radius
r_star = 0.20
theta = np.radians(15)
L = 5.0       # Total nozzle length

d = (3 / 2) * R_star * np.tan(theta)
x_plus = L - d
r_plus = (5 * d**2) / (12 * R_star) + r_star
a = 0.2  # Cylindrical section length 

b = 2 * (x_plus + ((r_plus - R_I) / np.tan(theta)) - a)
e = (R_I - r_star) / np.tan(theta)
c = e - (b / 2) - (5 * d / 8)

x_values = np.linspace(0, L, 1000) 
r_values = np.zeros_like(x_values) 

# Calculate radius for each segment
for i, x in enumerate(x_values):
    if 0 <= x <= a:
        r_values[i] = R_I
    elif a <= x <= a + b:
        r_values[i] = R_I - (b * np.tan(theta) / 2) * (((x - a) / b) ** 3) * (2 - ((x - a) / b))
    elif a + b <= x <= a + b + c:
        r_values[i] = R_I + (a + (b / 2)) * np.tan(theta) - x * np.tan(theta)
    else:
        r_values[i] = ((L - x)**2 / (12*R_star))*(6 - ((L - x) / d) ** 2) + r_star

# Reflect the nozzle shape about the vertical line at x = L
x_mirror = 2 * L - x_values  # Reflect x_values about x = L
r_mirror = r_values  # r_values remain the same

# Stretch the mirrored x-values by a factor of 3
x_mirror_stretched = L + 3 * (x_mirror - L)

# Stretch the mirrored r-values by a factor of 3, and subtract the original max radius value to avoid shifting
r_mirror_stretched = 0.2 + 2 * (r_mirror - np.min(r_mirror))  # Subtract the min value to keep symmetry

N = 200
x_mirror_stretched = x_mirror_stretched[N:]
r_mirror_stretched = r_mirror_stretched[N:]

# Reverse the mirrored arrays
x_mirror_stretched_reversed = x_mirror_stretched[::-1]
r_mirror_stretched_reversed = r_mirror_stretched[::-1]

# Concatenate the original arrays with the reversed mirrored arrays
x_temp = np.concatenate((x_values, x_mirror_stretched_reversed))
r_temp = np.concatenate((r_values, r_mirror_stretched_reversed))

# Resolution based on x-resolution
nx = 200 # t needs to match c * width/nx with c = 1,2,3,..
# Remove duplicate x-values
unique_indices = np.unique(x_temp, return_index=True)[1]
x_temp_unique = x_temp[unique_indices]
r_temp_unique = r_temp[unique_indices]

M = 150
interp_func = interp1d(x_temp_unique, r_temp_unique, kind='cubic')
# Generate new x values with only the desired number of points
x_combined_unfiltered = np.linspace(x_temp.min(), x_temp.max(), nx + M)
r_combined_unfiltered = interp_func(x_combined_unfiltered)

x_combined = x_combined_unfiltered[:-M]
r_combined = r_combined_unfiltered[:-M]

# dimensions of box (m)
width = 30
height = 3
depth = 3
# inner box (m)
wi = x_combined[nx-1] - x_combined[0]
hi = r_combined[nx-1]
di = r_combined[nx-1]
print("x_combined",len(x_combined))
print("r_combined",len(r_combined))
print("hi",r_combined[nx-1])
print("wi",x_combined[nx-1] - x_combined[0])
print("di",di)
# # Resolution based on x-resolution
# nx = 40  # t needs to match c * width/nx with c = 1,2,3,..
ny = int((height / width) * nx)
nz = int((depth/width) * nx)
mesh = create_3d_mesh(nx, ny, nz, width, height, depth) 

# Interface region: solid with thickness t (in m)
t = 2 * width / nx  # needs to match the discretization
print("t",t,'\n')
step = 0
# Define regions
for e in mesh.elements:
    x, y, z = mesh.calc_barycenter(e)
    x = round(x, 2)
    y = round(y, 2)
    z = round(z, 2)
    # Get the corresponding radius value
    if step < len(r_combined):  # Ensure step is within bounds
        r_with_thickness = r_combined[step] + t
        r_with_thickness = round(r_with_thickness, 2)
        y_plus_t = r_with_thickness
        z_plus_t = r_with_thickness
        curve_in = np.sqrt(y**2 + z**2)

        cover_in = np.sqrt(y**2 + z**2)
        # curve_out = np.sqrt(y_plus_t**2 + z_plus_t**2)
        # if y >= r_combined[step] and y <= r_with_thickness:
        if curve_in >= r_combined[step] and curve_in <= r_with_thickness:
            e.region = 'solid'
        elif curve_in < r_combined[step]:
            e.region = 'void'
        elif cover_in >= r_combined[nx -1] + 1.5 and cover_in <= r_combined[nx -1] + t + 1.5 and x <= 14:    # FOR CASING
            e.region = 'solid'
    step = (step + 1) % len(r_combined)  # Reset step if it reaches the end

box_width = []    
box_height = []
box_depth = []
nozzle_curve = [] 
casing = []
# horizontal_side = []
# vertical_side = []
step = 0
for i, n in enumerate(mesh.nodes):
  x, y, z = n
  if x >= x_combined[nx-1] - t and x <= x_combined[nx-1] and y <= r_combined[nx-1] and z <= r_combined[nx-1]:
    box_width.append(i)
  if y >= r_combined[nx-1] - t and y <= r_combined[nx-1] and x <= x_combined[nx-1] and z <= r_combined[nx-1]:
    box_height.append(i)       
  if z >= r_combined[nx-1] - t and z <= r_combined[nx-1] and x <= x_combined[nx-1] and y <= r_combined[nx-1]:
    box_depth.append(i)
  # if y >= 0 and y <= t and z >= 0 and x >= 0:
  #   horizontal_side.append(i)
  # if z >= 0 and z <= t and y >= 0 and x >=0:
  #   vertical_side.append(i)

  if np.sqrt(y**2 + z**2) >= r_combined[nx -1] + 1.5 and np.sqrt(y**2 + z**2) <= r_combined[nx -1] + t + 1.5 and x <= 14:     # FOR CASING
     casing.append(i)

  if step == nx:
    step = 0

#   elif step < nx and y == r_combined[step] and x >= x_combined[0] and x <= x_combined[nx-1]:
  elif step < nx and np.sqrt(y**2 + z**2) >= r_combined[step] and np.sqrt(y**2 + z**2) <= r_combined[step] + t and x >= x_combined[0] and x <= x_combined[nx-1]:
    nozzle_curve.append(i)
  step += 1

mesh.bc.append(('box_width',box_width))
mesh.bc.append(('box_height',box_height))
mesh.bc.append(('box_depth',box_depth))
mesh.bc.append(('nozzle_curve',nozzle_curve))
mesh.bc.append(('casing',casing))                                                                                                 # FOR CASING
# mesh.bc.append(('horizontal_side',horizontal_side))
# mesh.bc.append(('vertical_side',vertical_side))
x_coords, y_coords, z_coords = zip(*mesh.nodes)
print(f"x range: {min(x_coords)} to {max(x_coords)}")
print(f"y range: {min(y_coords)} to {max(y_coords)}")
print(f"z range: {min(z_coords)} to {max(z_coords)}")   

print('size of box_width/box_height/box_depth/nozzle_curve',len(box_width),len(box_height),len(box_depth),len(nozzle_curve))
    
b = 'box3d-t_' + str(t) + '-nx_' + str(nx) + '-ny_' + str(ny) + '-nz_' + str(nz) + '.mesh'   

write_ansys_mesh(mesh, b)
print('created ', b,'\n')
# hint: cfs -m <mesh> <problem> -g will just write the mesh without simulation/optimization for check in ParaView


# Initialize a set to store unique nodes belonging to solid regions
# solid_nodes = set()
# Loop through the elements and collect nodes in solid regions
# for e in mesh.elements:
#     if e.region == 'solid':  # Check if the element is in the 'solid' region
#         solid_nodes.update(e.nodes)  # Add the nodes of this element to the set

# Extract the coordinates of solid nodes
# solid_coordinates = [(mesh.nodes[node][0], mesh.nodes[node][1], mesh.nodes[node][2]) for node in solid_nodes]
# Print the coordinates in the desired format
# print(f"<nodes name=\"nozzle_curve\">")
# for x, y, z in solid_coordinates:
#     print(f"    <coord x=\"{x}\" y=\"{y}\" z=\"{z}\" />")
# print(f"</nodes>")
# Optionally, write these coordinates to a file for further inspection
# with open("solid_region_coordinates.xml", "w") as f:
#     for x, y, z in solid_coordinates:
#         f.write(f"<nodes name=\"nozzle_curve\">\n")
#         f.write(f"    <coord x=\"{x}\" y=\"{y}\" z=\"{z}\" />\n")
#         f.write(f"</nodes>\n")

# print(f"Solid region coordinates have been saved to 'solid_region_coordinates.xml'.")

# Extract coordinates for the nozzle curve and casing curve regions
# nozzle_curve_coordinates = [mesh.nodes[i] for i in nozzle_curve]
# casing_coordinates = [mesh.nodes[i] for i in casing]                                                                        # FOR CASING

# Write nozzle curve coordinates to an XML file
# with open("nozzle_curve_coordinates.xml", "w") as f:
#     f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#     f.write("<nodeList>\n")
#     for x, y, z in nozzle_curve_coordinates:
#         f.write("<nodes name=\"nozzle_curve\">\n")
#         f.write(f"    <coord x=\"{x:.6f}\" y=\"{y:.6f}\" z=\"{z:.6f}\" />\n")
#         f.write("</nodes>\n")
#     f.write("</nodeList>")

# Write casing curve coordinates to a separate XML file
# with open("casing_curve_coordinates.xml", "w") as f:
#     f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#     f.write("<nodeList>\n")
#     for x, y, z in casing_coordinates:
#         f.write("<nodes name=\"casing_curve\">\n")
#         f.write(f"    <coord x=\"{x:.6f}\" y=\"{y:.6f}\" z=\"{z:.6f}\" />\n")
#         f.write("</nodes>\n")
#     f.write("</nodeList>")

# with open("horizontal_side_coordinates.xml", "w") as f:
#     f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#     f.write("<nodeList>\n")
#     for x, y, z in casing_coordinates:
#         f.write("<nodes name=\"horizontal_side\">\n")
#         f.write(f"    <coord x=\"{x:.6f}\" y=\"{y:.6f}\" z=\"{z:.6f}\" />\n")
#         f.write("</nodes>\n")
#     f.write("</nodeList>")

# with open("vertical_side_coordinates.xml", "w") as f:
#     f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#     f.write("<nodeList>\n")
#     for x, y, z in casing_coordinates:
#         f.write("<nodes name=\"vertical_side\">\n")
#         f.write(f"    <coord x=\"{x:.6f}\" y=\"{y:.6f}\" z=\"{z:.6f}\" />\n")
#         f.write("</nodes>\n")
#     f.write("</nodeList>")

# Bottom vertical side
# with open("vertical_side_bottom_coordinates.xml", "w") as f:
#     f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#     f.write("<nodeList>\n")
#     for x, y, z in casing_coordinates:
#         f.write("<nodes name=\"vertical_side_bottom\">\n")
#         f.write(f"    <coord x=\"{x:.6f}\" y=\"{-y:.6f}\" z=\"{z:.6f}\" />\n")
#         f.write("</nodes>\n")
#     f.write("</nodeList>")

print("Nozzle curve and casing curve coordinates have been written to their respective XML files.")





#-------------UNCOMMENT BELOW ALL TOGETHER-------------

# # Inserting nodes into problem .xml file 
# from cfs_utils import *
# import xml.etree.ElementTree as ET

# # Define the namespace for easier querying
# namespaces = {
#     'cfs': 'http://www.cfs++.org/simulation'  # Define the proper namespace
# }

# # Main XML file
# # main_file = open_xml("3D_distributed_load.xml")
# main_file = open_xml("F3D_distributed_load.xml")
# # External node files (replace these with your actual file paths)
# # node_files = ["casing_curve_coordinates.xml", "nozzle_curve_coordinates.xml", "horizontal_side_coordinates.xml", "vertical_side_coordinates.xml", "vertical_side_bottom_coordinates.xml"]
# node_files = ["nozzle_curve_coordinates.xml"]
# # node_files = ["casing_curve_coordinates.xml","nozzle_curve_coordinates.xml"]
# # Find the nodeList element
# node_list = main_file.find('.//cfs:nodeList', namespaces)
# if node_list is None:
#     raise ValueError("No <nodeList> element found in the XML file.")
# print("Found <nodeList> element in main XML.")

# # Clear existing nodes (optional, if you don't want duplicates)
# for child in list(node_list):
#     node_list.remove(child)
# print("Existing nodes cleared.")

# # Iterate over the external node files and add nodes to the main XML
# for node_file_path in node_files:
#     values = open_xml(node_file_path)
#     print(f"Processing file: {node_file_path}")

#     # Find all <nodes> elements in the external files
#     # Check if namespace is applied or not
#     new_nodes = values.findall('.//nodes')  # Without using the namespace, try this first

#     if not new_nodes:
#         # If no nodes found, try using the namespace
#         new_nodes = values.findall('.//cfs:nodes', namespaces)
#         if not new_nodes:
#             print(f"No <nodes> elements found in {node_file_path}")
#         else:
#             print(f"Found {len(new_nodes)} <nodes> elements in {node_file_path}.")
#     else:
#         print(f"Found {len(new_nodes)} <nodes> elements in {node_file_path}.")

#     # Append each <nodes> element to the <nodeList> of the main XML
#     for nodes_element in new_nodes:
#         node_list.append(nodes_element)
#         # print(f"Appended <nodes> with name: {nodes_element.get('name')}")

# output_file = "F3D_distributed_load.xml"
# main_file.write(output_file, encoding="UTF-8", xml_declaration=True)
# print(f"Merged XML saved as: {output_file}\n")

print("RUN: cfs -m ",b," F3D_distributed_load")
