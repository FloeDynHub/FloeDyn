#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
""" Drawing floes from simulation output files """

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.animation as animation
from math import *
import h5py
import datetime
import os


Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=3600)


def _init(data, ax, begin=0):
    "initializes floes (Polygon creation) for matplotlib animation creation"
    for s in data.get("floe_outlines").values():
        polygon = Polygon(s[begin], True, color="white")
        ax.add_patch(polygon)
    return ax


def _update(num, data, ax, step=1, begin=0):
    "updates floes (Polygons) positions for matplotlib animation creation"
    indic = begin + step * num
    for i, s in enumerate(data.get("floe_outlines").itervalues()):
        ax.patches[i].set_xy(s[indic])
    ax.set_title("t = %s" % str(datetime.timedelta(seconds=int(data.get("time")[indic]))))
    ax.axis('equal')
    ax.relim()
    ax.autoscale_view(True,True,True)
    return ax,


def plot_floes(filename="out", step=1, make_video=False):
    "displays floes animation from hdf5 output file (simple plot or video creation)"
    hdf5_file_name = "io/%s.h5" % filename
    in_file    = h5py.File(hdf5_file_name, 'r')
    fig, ax = plt.subplots()
    ax = _init(in_file, ax)
    ax.axis('equal')
    ax.set_axis_bgcolor('#1B8EEF')
    anim = animation.FuncAnimation(fig, _update, in_file.get("time").size / step, fargs=(in_file, ax, step),
        interval=50, blit=False)

    if make_video:
        anim.save('%s.mp4' % filename, writer=writer)
    else:
        plt.show()


def plot_last_rec(filename):
    "displays floes at last time recorded in output file"
    if not filename:
        filename = "out"
    hdf5_file_name = "io/%s.h5" % filename
    file    = h5py.File(hdf5_file_name, 'r')
    fig, ax = plt.subplots()
    ax = _init(file, ax)
    ax.axis('equal')
    ax.set_axis_bgcolor('#1B8EEF')
    ax, = _update(file.get("time").size - 1, file, ax)
    plt.show();


# Parallel video creation :

def make_partial_floe_video(out_filename, data_chunk):
    "creates video directly from data"
    fig, ax = plt.subplots()
    ax = _init(data_chunk, ax)
    ax.axis('equal')
    ax.set_axis_bgcolor('#1B8EEF')
    anim = animation.FuncAnimation(fig, _update, len(data_chunk.get("time")), fargs=(data_chunk, ax), # TODO change chunk format to get size properly
        interval=50, blit=False)
    anim.save(out_filename, writer=writer)


def make_partial_floe_video_helper(t):
    return make_partial_floe_video(*t)


from multiprocessing import Pool, Process, cpu_count
from subprocess import call
def plot_floes_fast(filename="out", step=1):
    "Creates a video from hdf5 file output, using multiprocessing to go faster"
    hdf5_file_name = "io/%s.h5" % filename
    data_file    = h5py.File(hdf5_file_name, 'r')
    nb_step = len(data_file.get("time"))
    # Create multiple partial videos
    nb_process = cpu_count()
    trunk_size = nb_step / nb_process
    temp_dir = 'plot_tmp'
    call(['mkdir', temp_dir])
    partial_file_names = ["{}/{}.mpg".format(temp_dir, i) for i in range(nb_process)]
    L = [(
            partial_file_names[i],
            {
                "floe_outlines" : {k : dataset[i * trunk_size : (i+1) * trunk_size: step]
                    for k, dataset in data_file.get("floe_outlines").iteritems()},
                "time" : data_file.get("time")[i * trunk_size : (i+1) * trunk_size: step]
            }
         ) for i in range(nb_process)]
    p = Pool(nb_process)
    p.map(make_partial_floe_video_helper, L)
    p.close()
    p.join()
    # Concat all partial video
    out_filename = '{}.mp4'.format(filename)
    call(['ffmpeg',  '-i', 'concat:{}'.format("|".join(partial_file_names)), '-c', 'copy', out_filename])
    call(['rm', '-r', temp_dir])


