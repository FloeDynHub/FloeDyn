import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation
from math import *
import h5py

# def mkfloe(nb=15, d=1, t=(0,0)):
#     return [(d*cos(2*i*pi/nb) + t[0], d*sin(2*i*pi/nb) + t[1]) for i in range(nb+1)]



def init(data, ax):
    datasets = data[0].values()
    for i in range(len(datasets)):
        polygon = Polygon(datasets[i], True, color="white")
        ax.add_patch(polygon)
    return ax


def update(num, data, ax):
    datasets = data[num].values()
    for i in range(len(datasets)):
        polygon = Polygon(datasets[i], True, color="white")
        ax.patches[i].set_xy(polygon.get_xy())
    ax.set_title("t = %s" % data[num].name[1:])
    return ax,


def plot_floes(filename):
    if not filename:
        filename = "out"
    hdf5_file_name = "out/%s.h5" % filename
    file    = h5py.File(hdf5_file_name, 'r')
    fig, ax = plt.subplots()
    groups = file.values()
    ax = init(groups, ax)
    ax.axis('equal')
    ax.set_axis_bgcolor('#1B8EEF')
    anim = animation.FuncAnimation(fig, update, len(groups), fargs=(groups, ax),
        interval=50, blit=False)
    plt.show()

