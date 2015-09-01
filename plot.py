import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation
from math import *
import h5py


Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=3600)


def init(data, ax):
    for s in data.get("floe_outlines").values():
        polygon = Polygon(s[0], True, color="white")
        ax.add_patch(polygon)
    return ax


def update(num, data, ax):
    i=0
    for s in data.get("floe_outlines").values():
        polygon = Polygon(s[num], True, color="white")
        ax.patches[i].set_xy(polygon.get_xy())
        i+=1
    ax.set_title("t = %s" % data.get("time")[num])
    ax.axis('equal')
    ax.relim()
    ax.autoscale_view(True,True,True)
    return ax,


def plot_floes(filename, make_video=False):
    if not filename:
        filename = "out"
    hdf5_file_name = "io/%s.h5" % filename
    file    = h5py.File(hdf5_file_name, 'r')
    fig, ax = plt.subplots()
    groups = file.values()
    ax = init(file, ax)
    ax.axis('equal')
    ax.set_axis_bgcolor('#1B8EEF')
    anim = animation.FuncAnimation(fig, update, file.get("time").size, fargs=(file, ax),
        interval=50, blit=False)

    if make_video:
        anim.save('floes.mp4', writer=writer)
    else:
        plt.show()


