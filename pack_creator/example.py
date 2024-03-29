from pack_creator import (
    Floe, State, circle_floe_shape, translate_floe_group, write_input_file,
    rectangle_floe_shape
)
import math

"""
Create basic pack with 5 regular polygon floes

Commands (from root directory) :
Creation : python3 pack_creator/example.py
Plot : python3 plot img -f io/inputs/in_example_pack.h5 --nocolor
"""

shapes = {
"F" : [(25.86, 0.0), (25.86, 47.57), (75.0, 47.57), (75.0, 59.93), (25.86, 59.93), (25.86, 92.36), (82.64, 92.36), (82.64, 104.71), (12.0, 104.71), (12.0, 0.0)],
"L" : [(76.14, 0.0), (76.14, 12.36), (24.57, 12.36), (24.57, 104.71), (10.71, 104.71), (10.71, 0.0)],
"O" : [(7.07, 37.21), (19.07, 12.5), (42.36, -1.79), (57.14, -1.79), (70.79, -1.79), (94.36, 11.14), (107.21, 36.29), (107.21, 52.21), (107.21, 67.93), (95.0, 92.71), (71.71, 106.57), (57.21, 106.57), (35.07, 106.57), (7.07, 77.07), (7.07, 51.0)],
"E" : [(89.71, 0.0), (89.71, 12.36), (25.43, 12.36), (25.43, 48.0), (83.36, 48.0), (83.36, 60.29), (25.43, 60.29), (25.43, 92.36), (87.29, 92.36), (87.29, 104.71), (11.57, 104.71), (11.57, 0.0)],
"D" : [(49.07, 0.0), (58.57, 0.0), (73.21, 3.57), (83.57, 10.5), (92.29, 22.57), (97.86, 41.0), (97.86, 52.93), (97.86, 66.93), (89.64, 88.71), (81.36, 95.71), (75.0, 101.14), (66.0, 103.21), (59.57, 104.71), (47.36, 104.71), (11.29, 104.71), (11.29, 0.0)],
"Y" : [(54.64, 0.0), (54.64, 44.36), (96.43, 104.71), (80.29, 104.71), (60.0, 73.93), (53.29, 63.64), (48.57, 55.43), (43.64, 64.29), (37.93, 73.14), (17.29, 104.71), (0.43, 104.71), (40.79, 44.36), (40.79, 0.0)],
"N" : [(24.43, 0.0), (24.43, 82.29), (79.43, 0.0), (93.64, 0.0), (93.64, 104.71), (80.36, 104.71), (80.36, 22.5), (25.36, 104.71), (11.14, 104.71), (11.14, 0.0)]
}


list_f = []

# Create rectangle obstacle contour

list_f.append(
    Floe([(-140, -200), (-140, 190), (840, 190), (840, -190), (-140, -190), (-140, -200), (850, -200), (850, 200), (-150, 200), (-150, -200)], State(pos=[0, 0]))
)

# Write FLOEDYN
pos_x = 0
for letter in "FLOEDYN":
    list_f += [
        Floe(shapes[letter], State(pos=[pos_x, 0])),
    ]
    pos_x += 100

# Add curcular shapes with increasing number of vertices
pos_x = 0
for nb in range(3, 11):
    list_f.append(
        Floe(circle_floe_shape(nb * 5, nb_vertices=nb), State(pos=[pos_x, -100]))
    )
    pos_x += nb * 15

# Translate all pack if needed
# translate_floe_group(list_f, 1000, 1000)

write_input_file(list_f, "in_example_pack")
