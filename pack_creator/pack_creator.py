import math
import h5py
from shapely.geometry import Polygon


class State:
    def __init__(self, pos=[0, 0], theta=0, speed=[0, 0], rot=0):
        self.pos = pos
        self.theta = theta
        self.speed = speed
        self.rot = rot

    def get_dataset_slice(self):
      return [
        self.pos[0],
        self.pos[1],
        self.theta,
        self.speed[0],
        self.speed[1],
        self.rot,
        0,
        self.pos[0],
        self.pos[1]
      ]


class Floe:
  def __init__(self, shape, state=None):
    """
    Shape must be described in counter clockwise order !
    """
    # re-center shape on mass center
    c = Polygon(shape).centroid
    # self.shape = shape # As point list [(x, y), ...]
    # self.state = state or State() # State object
    self.shape = [(x - c.x, y - c.y) for x, y in shape]
    if state:
        x, y = state.pos
        state.pos = [x + c.x, y + c.y]
        self.state = state
    else:
        self.state = State()


def get_window(floes_list):
    max_x = max(pt[0] + floe.state.pos[0] for floe in floes_list for pt in floe.shape)
    max_y = max(pt[1] + floe.state.pos[1] for floe in floes_list for pt in floe.shape)
    min_x = min(pt[0] + floe.state.pos[0] for floe in floes_list for pt in floe.shape)
    min_y = min(pt[1] + floe.state.pos[1] for floe in floes_list for pt in floe.shape)
    width = max(max_x - min_x, max_y - min_y) * 1.01
    x0 = (min_x + max_x) / 2
    y0 = (min_y + max_y) / 2
    return [x0 - width / 2 , x0 + width / 2, y0 - width / 2, y0 + width / 2]


def write_input_file(floes_list, filename):
    """Creates hdf5 input file for floe list"""
    with h5py.File(f"io/inputs/{filename}.h5", "w") as f:
        states_dset = f.create_dataset("floe_states", (1, len(floes_list), 9), dtype='float64') # 1 = t0, nb_floe, state_size
        grp = f.create_group("floe_shapes")
        for i, floe in enumerate(floes_list):
            shape_dset = grp.create_dataset(f"{i}", (len(floe.shape), 2), dtype='float64') # nb_points, 2 (space dim)
            shape_dset[...] = floe.shape
            # shape_dset.attrs['thickness'] = 1.2
            states_dset[0, i, ...] = floe.state.get_dataset_slice()
        win_dset = f.create_dataset("window", (4, ), dtype='float64')
        win_dset[...] = get_window(floes_list)


def circle_floe_shape(radius=1, nb_vertices=25):
  """Regular polygon floe, centered on (0, 0)"""
  f = []
  for k in range(nb_vertices):
    f.append((
      radius * math.cos((2 * k * math.pi) / nb_vertices), 
      radius * math.sin((2 * k * math.pi) / nb_vertices)
    ))
  return f

def rectangle_floe_shape(xlen, ylen):
   x = xlen / 2
   y = ylen / 2
   return [(x, -y), (x, y), (-x, y), (-x, -y)]

def translate_floe_group(floes, x, y):
    for floe in floes:
        floe.state.pos[0] += x
        floe.state.pos[1] += y

def complete_shape(shape, max_edge_length=10, close=True):
    """
    Complete the shape with points to have edges of length < max_edge_length
    """
    new_shape = []
    max_iter = len(shape) if close else len(shape) - 1
    for i in range(max_iter):
        new_shape.append(shape[i])
        next_i = (i + 1) % len(shape)
        edge = [shape[i], shape[next_i]]
        edge_length = math.sqrt((edge[1][0] - edge[0][0])**2 + (edge[1][1] - edge[0][1])**2)
        nb_points = int(edge_length / max_edge_length)
        for j in range(1, nb_points):
            new_shape.append((edge[0][0] + j * (edge[1][0] - edge[0][0]) / nb_points, edge[0][1] + j * (edge[1][1] - edge[0][1]) / nb_points))
    return new_shape

def resize_pack(list_floe, resize_coeff):
    for floe in list_floe:
        floe.shape = [(resize_coeff * x, resize_coeff * y) for x, y in floe.shape]
        floe.state.pos = [resize_coeff * x for x in floe.state.pos]

def translate_pack(list_floe, x, y):
    for floe in list_floe:
        floe.state.pos = [floe.state.pos[0] + x, floe.state.pos[1] + y]

def load_h5_input(filename, list_floe):
    with h5py.File(filename, "r") as f:
        grp = f["floe_shapes"]
        for i in range(len(grp)):
            shape = grp[f"{i}"][...]
            state = State(
                pos=f["floe_states"][0, i, :2],
                theta=f["floe_states"][0, i, 2],
                speed=f["floe_states"][0, i, 3:5],
                rot=f["floe_states"][0, i, 5]
            )
            list_floe.append(Floe(shape, state))
        win = f["window"][...]
        # translate_pack(list_floe, -win[0], -win[2])

def display_floes(list_floe):
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    for floe in list_floe:
        ax.plot([x + floe.state.pos[0] for x, y in floe.shape], [y + floe.state.pos[1] for x, y in floe.shape], 'b')
    ax.set_aspect('equal')
    plt.show()