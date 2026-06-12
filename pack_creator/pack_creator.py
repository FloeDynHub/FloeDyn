import math
import h5py
from shapely.geometry import Polygon, Point


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
    def __init__(self, shape, state=None, obstacle=False, thickness=None, cw=None):
        """
        Shape must be described in counter clockwise order !
        """
        self.shape = shape
        self.obstacle = obstacle
        self.thickness = thickness
        self.cw = cw
        self.state = state if state else State()
        self.recenter_shape()

    def recenter_shape(self):
        """
        Recenter shape on mass center
        """
        c = Polygon(self.shape).centroid
        self.shape = [(x - c.x, y - c.y) for x, y in self.shape]
        x, y = self.state.pos
        self.state.pos = [x + c.x, y + c.y]
    
    def reduce_shape(self, nb_points):
        """
        Reduce the shape to nb_points by removing points
        """
        self.shape = reduce_shape(self.shape, nb_points)
        self.recenter_shape()
    
    def max_diameter(self):
        """
        Returns the maximum diameter of the floe shape
        = largest possible distance between any two points on the polygon's boundary
        """
        max_diam = 0
        for i in range(len(self.shape)):
            for j in range(i + 1, len(self.shape)):
                p1 = Point(self.shape[i])
                p2 = Point(self.shape[j])
                dist = p1.distance(p2)
                if dist > max_diam:
                    max_diam = dist
        return max_diam
    
    def resize(self, resize_coeff):
        """
        Resize floe around its mass center by a given coefficient.
        """
        resize_floe(self, resize_coeff)
    

    def __repr__(self):
        return f"Floe(shape={self.shape}, state={self.state}, obstacle={self.obstacle}, thickness={self.thickness}, cw={self.cw})"

def get_window(floes_list):
    max_x = max(pt[0] + floe.state.pos[0] for floe in floes_list for pt in floe.shape)
    max_y = max(pt[1] + floe.state.pos[1] for floe in floes_list for pt in floe.shape)
    min_x = min(pt[0] + floe.state.pos[0] for floe in floes_list for pt in floe.shape)
    min_y = min(pt[1] + floe.state.pos[1] for floe in floes_list for pt in floe.shape)
    width = max(max_x - min_x, max_y - min_y) * 1.01
    x0 = (min_x + max_x) / 2
    y0 = (min_y + max_y) / 2
    return [x0 - width / 2 , x0 + width / 2, y0 - width / 2, y0 + width / 2]


def write_input_file(floes_list, filename, window=None):
    """Creates hdf5 input file for floe list"""
    with h5py.File(f"io/inputs/{filename}.h5", "w") as f:
        states_dset = f.create_dataset("floe_states", (1, len(floes_list), 9), dtype='float64') # 1 = t0, nb_floe, state_size
        grp = f.create_group("floe_shapes")
        for i, floe in enumerate(floes_list):
            shape_dset = grp.create_dataset(f"{i}", (len(floe.shape), 2), dtype='float64') # nb_points, 2 (space dim)
            shape_dset[...] = floe.shape
            if floe.obstacle:
                shape_dset.attrs['obstacle'] = 1
            if floe.thickness:
                shape_dset.attrs['thickness'] = floe.thickness
            if floe.cw:
                shape_dset.attrs['C_w'] = floe.cw
            states_dset[0, i, ...] = floe.state.get_dataset_slice()
        win_dset = f.create_dataset("window", (4, ), dtype='float64')
        win_dset[...] =  window or get_window(floes_list)


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

def reduce_shape(shape, nb_points):
    """
    Reduce the shape to nb_points by removing points
    """
    if len(shape) <= nb_points:
        return shape
    step = len(shape) // nb_points
    return [shape[i] for i in range(0, len(shape), step)]

from shapely.geometry import Point

def resize_pack(list_floe, resize_coeff):
    for floe in list_floe:
        floe.shape = [(resize_coeff * x, resize_coeff * y) for x, y in floe.shape]
        floe.state.pos = [resize_coeff * x for x in floe.state.pos]

def resize_floe(floe, resize_coeff):
    # Resize floe around its mass center
    c = Polygon(floe.shape).centroid
    floe.shape = [(c.x + resize_coeff * (x - c.x), c.y + resize_coeff * (y - c.y)) for x, y in floe.shape]

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
            obstacle = bool(grp[f"{i}"].attrs.get('obstacle', 0))
            thickness = grp[f"{i}"].attrs.get('thickness', None)
            cw = grp[f"{i}"].attrs.get('C_w', None)
            floe = Floe(shape.tolist(), state, obstacle, thickness, cw)
            list_floe.append(floe)
        win = f["window"][...]
        # translate_pack(list_floe, -win[0], -win[2])
    return list(win)

# ---- Obstacle-wall building helpers (OPTIMJAM scene design) -------------------------------------
# Walls/floors are made of SEPARATE obstacle blocks (per-block force readout, MPI-decomposable).
# Lessons from the in_dead_end INTER-storm forensics baked in here:
#  - block junctions must not present 90-degree corners on the contact face (a sliding circle's
#    polygon vertex catches them -> interpenetration -> dt collapse): contact-face corners are
#    CHAMFERED at 45 degrees, so each junction is a shallow V-groove ridden as a sloped contact;
#  - the inter-block gap only exists so the detector sees no contact at t=0 (obstacles never
#    move): keep it tiny;
#  - faces are subdivided so OptimizedFloe builds tight local disks (sharper proximity distances).

def unit(p, q):
    """Unit vector from p to q."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    n = math.hypot(dx, dy)
    return (dx / n, dy / n)


def chamfered_block(inner_a, inner_b, outer_a, outer_b, chamfer=0.5, subdiv=50):
    """Quadrilateral obstacle block with 45-degree chamfered corners on its INNER (contact) face
    and subdivided faces. inner_a->inner_b is the face the floes touch; outer_* are the matching
    back corners. Returns a CCW point list (Floe requirement)."""
    u = unit(inner_a, inner_b)    # along the contact face
    wa = unit(inner_a, outer_a)   # into the wall, a-end
    wb = unit(inner_b, outer_b)   # into the wall, b-end
    face_len = math.hypot(inner_b[0] - inner_a[0], inner_b[1] - inner_a[1])
    c = min(chamfer, face_len / 4)  # never let chamfers eat the whole face
    pa_side = (inner_a[0] + c * wa[0], inner_a[1] + c * wa[1])
    pa_face = (inner_a[0] + c * u[0], inner_a[1] + c * u[1])
    pb_face = (inner_b[0] - c * u[0], inner_b[1] - c * u[1])
    pb_side = (inner_b[0] + c * wb[0], inner_b[1] + c * wb[1])
    face = complete_shape([pa_face, pb_face], subdiv, close=False) + [pb_face]
    back = complete_shape([outer_b, outer_a], subdiv, close=False) + [outer_a]
    ring = [pa_side] + face + [pb_side] + back
    area2 = sum(ring[i][0] * ring[(i + 1) % len(ring)][1]
                - ring[(i + 1) % len(ring)][0] * ring[i][1] for i in range(len(ring)))
    if area2 < 0:
        ring.reverse()
    return ring


def blocks_along(polyline, thickness, offset='normal', gap=0.05, chamfer=0.5, subdiv=50):
    """One chamfered obstacle block per segment of `polyline` (the contact side).

    offset selects the direction in which the block body extends behind the contact face:
     - 'normal': each segment's CLOCKWISE normal — correct when the floes lie on the LEFT of the
       travel direction. Safe for straight parts and for curves bending AWAY from the offset side
       (bowl/cup/funnel: interior convex). Curves bending TOWARD the offset side would make
       adjacent block backs overlap — use a fixed direction there instead.
     - (dx, dy): fixed direction for every block — side edges stay parallel whatever the
       curvature (use for wavy/sinusoidal walls).
    Returns a list of obstacle Floe."""
    blocks = []
    for p, q in zip(polyline[:-1], polyline[1:]):
        u = unit(p, q)
        ia = (p[0] + (gap / 2) * u[0], p[1] + (gap / 2) * u[1])
        ib = (q[0] - (gap / 2) * u[0], q[1] - (gap / 2) * u[1])
        if offset == 'normal':
            d = (u[1], -u[0])  # clockwise normal of the travel direction
        else:
            n = math.hypot(offset[0], offset[1])
            d = (offset[0] / n, offset[1] / n)
        oa = (ia[0] + thickness * d[0], ia[1] + thickness * d[1])
        ob = (ib[0] + thickness * d[0], ib[1] + thickness * d[1])
        blocks.append(Floe(chamfered_block(ia, ib, oa, ob, chamfer, subdiv),
                           State(pos=[0, 0]), obstacle=True))
    return blocks


def check_no_initial_overlap(list_f):
    """Sanity check: no obstacle-obstacle nor floe-obstacle intersection at t=0 (the detector
    would never recover from a pre-existing interpenetration)."""
    polys = [Polygon([(x + f.state.pos[0], y + f.state.pos[1]) for x, y in f.shape])
             for f in list_f]
    obs_idx = [i for i, f in enumerate(list_f) if f.obstacle]
    flo_idx = [i for i, f in enumerate(list_f) if not f.obstacle]
    bad = [(i, j) for k, i in enumerate(obs_idx) for j in obs_idx[k + 1:]
           if polys[i].intersects(polys[j])]
    bad += [(i, j) for i in flo_idx for j in obs_idx if polys[i].intersects(polys[j])]
    if bad:
        raise SystemExit(f"ERROR: intersecting pairs at t=0: {bad[:10]}")
    print("Sanity check OK: no intersecting pair at t=0")


def save_preview(list_floe, path):
    """Headless PNG preview of the pack (same drawing as display_floes)."""
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 14))
    for floe in list_floe:
        xs = [floe.state.pos[0] + x for x, y in floe.shape] + [floe.state.pos[0] + floe.shape[0][0]]
        ys = [floe.state.pos[1] + y for x, y in floe.shape] + [floe.state.pos[1] + floe.shape[0][1]]
        ax.plot(xs, ys, 'r' if floe.obstacle else 'b', linewidth=0.6)
    ax.set_aspect('equal')
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)
    print(f"Preview saved: {path}")


def display_floes(list_floe, axes=None):
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    for floe in list_floe:
        cos_theta = math.cos(floe.state.theta)
        sin_theta = math.sin(floe.state.theta)
        shape_x = [
            floe.state.pos[0] + x * cos_theta - y * sin_theta
            for x, y in floe.shape
        ] + [floe.state.pos[0] + floe.shape[0][0] * cos_theta - floe.shape[0][1] * sin_theta]
        shape_y = [
            floe.state.pos[1] + x * sin_theta + y * cos_theta
            for x, y in floe.shape
        ] + [floe.state.pos[1] + floe.shape[0][0] * sin_theta + floe.shape[0][1] * cos_theta]
        ax.plot(shape_x, shape_y, 'b')
    ax.set_aspect('equal')
    if axes: # Explicit Axes limits
        ax.set_xlim(axes[0], axes[1])
        ax.set_ylim(axes[2], axes[3])
    plt.show()
