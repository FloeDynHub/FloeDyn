#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
""" Drawing floes from simulation output files """

import numpy as np
# import matplotlib
# matplotlib.use("GTKAgg")
from matplotlib import transforms, cm, animation, colors, gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection, PolyCollection

from math import *
import h5py
import datetime
import os

def filename_without_extension(path):
    return os.path.splitext(os.path.basename(path))[0]

def _init2(data, ax, begin=0):
    "initializes floes (Polygon creation) for matplotlib animation creation (v3 : shapes)"
    floes_collec = PolyCollection(data.get("floe_shapes"), linewidths=0.2,
                                 cmap=cm.YlOrRd, norm=colors.PowerNorm(gamma=0.3))
                                # cmap=cm.Paired)#MPI
    if getattr(OPTIONS, "color", True):
        floes_collec.set_array(data.get("impulses")[begin])
        floes_collec.set_clim([0, data["MAX_IMPULSE"]])
        plt.colorbar(floes_collec, ax = ax, fraction = 0.05) # display color bar
    else:
        floes_collec.set_array(np.zeros(len(data.get("floe_shapes"))))
    # ax.axes.get_yaxis().set_visible(False) # remove x axis informations
    ax.add_collection(floes_collec)

    if getattr(OPTIONS, "ghosts", False):
        ghosts_collec = PolyCollection(data.get("floe_shapes") * 8, linewidths=0.2, alpha=0.4,
                                     cmap=cm.YlOrRd, norm=colors.PowerNorm(gamma=0.3))
                                    # cmap=cm.Paired)#MPI
        ghosts_collec.set_clim(floes_collec.get_clim())
        ghosts_collec.set_array(np.zeros(len(data.get("floe_shapes")) * 8))
        ax.add_collection(ghosts_collec)

    # Display generator window
    w = data.get("WIN_WIDTH", 0)
    if w:
        w = data.get("WIN_WIDTH")
        ax.add_patch( Polygon(((w, w), (-w, w), (-w, -w), (w, -w)),
                      True, facecolor=None, alpha=0.2,
                      edgecolor="white", linewidth=0.2))
    # END Display generator window

    return ax


def _init2_dual(data, ax1, ax2, begin=0):
    "initializes floes (Polygon creation) for matplotlib animation creation (v3 : shapes)"
    ax1 = _init2(data, ax1)
    first_idx, _ = data["total_trunk_indices"]
    x = np.array(data.get("total_time")[0:first_idx])
    ys = [np.array([data.get("total_impulses")[i][idx] for i in range(len(x))]) for idx in data.get("special_indices", [])]
    for y in ys:
        ax2.plot(x,y)
    ax2.axis([data.get("total_time")[0], data.get("total_time")[-1], -100, data["MAX_IMPULSE"]])
    ax2.set_title("Norm of the contact impulses applied to the obstacle")
    # ax.set_xlabel("text") # bottom
    return ax1, ax2


def _init1(data, ax, begin=0):
    "initializes floes (Polygon creation) for matplotlib animation creation"
    for (i,s) in data.get("floe_outlines").iteritems():
        special = i=="1500"
        polygon = Polygon(s[begin],
                          True,
                          facecolor="#FFF7ED" if not special else "#FA9A50",
                          linewidth=0.2,
                          hatch="" if not special else "x")
        ax.add_patch(polygon)
    return ax

def _update2(num, data, ax, step=1, begin=0):
    """updates floes (Polygons) positions for matplotlib animation creation
    (version 3 : collection and compute outlines transforms"""
    indic = begin + step * num
    opt_ghosts, opt_color, opt_follow = (
        getattr(OPTIONS, opt, default) for opt, default in (("ghosts", False), ("color", True), ("follow", False))
    )

    verts = _transform_shapes(data.get("floe_shapes"), data.get("floe_states")[indic], opt_follow)
    ax.collections[0].set_verts(verts)

    if opt_ghosts:
        w = data.get("WIN_WIDTH", 0) * 2
        # ghosts_trans = [(w, 0), (w, w), (0, w), (-w, w), (-w, 0), (-w, -w), (0, -w), (w, -w)]
        ghosts_trans = [(i * w, j * w) for i in range(-1, 2) for j in range(-1, 2) if not i==j==0]
        # ghosts_trans = [(i * w, j * w) for i in range(-2, 3) for j in range(-2, 3) if not i==j==0] #MPI
        ghosts_verts = [v for trans in ghosts_trans for v in translate_group(verts, trans)]
        ax.collections[1].set_verts(ghosts_verts)

    if opt_color:
        ax.collections[0].set_array(data.get("impulses")[indic])
        if opt_ghosts:
            ax.collections[1].set_array(np.tile(data.get("impulses")[indic], 8))
    ax.set_title("t = {}".format(str(datetime.timedelta(seconds=int(data.get("time")[indic])))))
    if not STATIC_AXES:
        # ax.axis('equal')
        ax.relim()
        ax.autoscale_view(True,True,True)
    if opt_follow:
        ax.set_xlim(data.get("mass_center")[indic][0] + STATIC_AXES[0], data.get("mass_center")[indic][0] +  STATIC_AXES[1])
        ax.set_ylim(data.get("mass_center")[indic][1] + STATIC_AXES[2], data.get("mass_center")[indic][1] +  STATIC_AXES[3])
    # # MPI zones
    # ax.collections[0].set_clim([0, 9])
    # ax.collections[1].set_clim([0, 9])
    # W, H = ax.get_xlim()[1] - ax.get_xlim()[0], ax.get_ylim()[1] - ax.get_ylim()[0]
    # N = 6
    # WZ, HZ = W/N, H/N
    # ax.collections[0].set_array(np.array([(int((x[0] - ax.get_xlim()[0]) / WZ) + int((x[1] - ax.get_ylim()[0]) / HZ) * N)%9 for x in data.get("floe_states")[indic]]))
    # ghosts_trans = [(i * w, j * w) for i in range(-2, 3) for j in range(-2, 3) if not i==j==0]
    # gstates = [(x[0] + i, x[1] + j) for i,j in ghosts_trans for x in data.get("floe_states")[indic]]
    # ax.collections[1].set_array(np.array([(int((x[0] - ax.get_xlim()[0]) / WZ) + int((x[1] - ax.get_ylim()[0]) / HZ) * N)%9 for x in gstates]))
    return ax,

def _update2_dual(num, data, ax1, ax2, step=1, begin=0):
    first_idx, _ = data["total_trunk_indices"]
    indic = first_idx + num
    ax1, = _update2(num, data, ax1)
    for i, idx in enumerate(data.get("special_indices")):
        ax2.lines[i].set_xdata(np.append(ax2.lines[i].get_xdata(), data.get("total_time")[indic]))
        ax2.lines[i].set_ydata(np.append(ax2.lines[i].get_ydata(), data.get("total_impulses")[indic][idx]))
    return ax1, ax2


def _transform_shapes(floe_shapes, floe_states, follow=False):
    resp = floe_shapes
    def rotation_mat(theta):
        return np.array([[np.cos(theta), -np.sin(theta)], 
                         [np.sin(theta),  np.cos(theta)]])
    rots = [rotation_mat(x[2]) for x in floe_states]
    resp = [np.transpose(np.dot(rot, np.transpose(shape))) for rot, shape in zip(rots, floe_shapes)]
    pos_ids = (0,1) if follow else (7,8) # (7,8) contains translated position in fixed initial window
    resp = [np.add(shape, np.repeat([[x[pos_ids[0]], x[pos_ids[1]]]], len(shape), axis=0)) for x, shape in zip(floe_states, resp)]
    return resp

def translate_group(vertices, trans):
    return [np.add(shape, np.repeat([[trans[0], trans[1]]], len(shape), axis=0)) for shape in vertices]

def _update1(num, data, ax, step=1, begin=0):
    "updates floes (Polygons) positions for matplotlib animation creation"
    indic = begin + step * num
    for i, s in enumerate(data.get("floe_outlines").itervalues()):
        ax.patches[i].set_xy(s[indic])
    ax.set_title("t = {}".format(str(datetime.timedelta(seconds=int(data.get("time")[indic])))))
    if not STATIC_AXES:
        ax.axis('equal')
        ax.relim()
        ax.autoscale_view(True,True,True)
    return ax,


def plot_floes(filename="out", step=1, version=2, make_video=1, **kwargs):
    "displays floes animation from hdf5 output file (simple plot or video creation)"
    data_file    = h5py.File(filename, 'r')
    # Read Usefull data from file
    data_global = _get_useful_datas(data_file, step, version)
    fig, ax = plt.subplots()
    if getattr(OPTIONS, "hd", False):
        fig.set_size_inches(20, 15)
        fig.set_dpi(100)
    init, update = anim_fcts(version)
    ax = init(data_global, ax)
    if not STATIC_AXES:
        ax.axis('equal') # automatic scale
    else:
        ax.axis(STATIC_AXES)
    ax.set_axis_bgcolor('#162252')
    anim = animation.FuncAnimation(
        fig, update, len(data_global.get("time")), fargs=(data_global, ax),
        interval=1, blit=False)
    if make_video:
        anim.save('io/videos/{}.mp4'.format(filename_without_extension(filename)), writer=writer)
    else:
        plt.show()


def plot_one(filename, **kwargs):
    "displays floes at last time recorded in output file"
    init, update = anim_fcts(2)
    file    = h5py.File(filename, 'r')
    fig, ax = plt.subplots()
    idx = int(raw_input("Time index (-1 for last one) : "))
    num = idx if idx != -1 else data_global.get("time").size - 1 # TODO data_global not declared !!
    
    data_global = _get_useful_datas(file, 0, 2, num)
    ax = init(data_global, ax)
    ax.axis('equal')
    ax.set_axis_bgcolor('#162252') 
    ax, = update(0, data_global, ax)
    plt.show()


#############################
# Parallel video creation : #
#############################

def anim_fcts(version):
    d = {1: (_init1,_update1), 2: (_init2,_update2)}
    return d[version]

def make_partial_floe_video(out_filename, data_chunk, version):
    "creates video directly from data"
    print 1
    fig, ax = plt.subplots()
    print 2
    if getattr(OPTIONS, "hd", False):
        fig.set_size_inches(20, 15)
        fig.set_dpi(100)
    print 3
    init, update = anim_fcts(version)
    ax = init(data_chunk, ax)
    print 4
    if not STATIC_AXES:
        ax.axis('equal') # automatic scale
    else:
        ax.axis(STATIC_AXES)
    print 5
    ax.set_axis_bgcolor('#162252')
    anim = animation.FuncAnimation(
        fig, update, len(data_chunk.get("time")), fargs=(data_chunk, ax),
        interval=1, blit=False)
    print 6
    print out_filename
    anim.save(out_filename, writer=writer)

def make_partial_floe_video_helper(t):
    return make_partial_floe_video(*t)


def make_partial_floe_video_dual_plot(out_filename, data_chunk, version):
    "creates video directly from data, and adds obstacle impulse subplot under floes"
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    fig = plt.figure()
    if getattr(OPTIONS, "hd", False):
        fig.set_size_inches(16, 12)
        fig.set_dpi(100)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    init, update = _init2_dual, _update2_dual
    ax1, ax2 = init(data_chunk, ax1, ax2)
    if not STATIC_AXES:
        ax1.axis('equal') # automatic scale
    else:
        ax1.axis(STATIC_AXES)
    ax1.set_axis_bgcolor('#162252')
    anim = animation.FuncAnimation(
        fig, update, len(data_chunk.get("time")), fargs=(data_chunk, ax1, ax2),
        interval=50, blit=False)
    anim.save(out_filename, writer=writer)

def make_partial_floe_video_dual_plot_helper(t):
    return make_partial_floe_video_dual_plot(*t)


from multiprocessing import Pool, Process, cpu_count
from subprocess import call
def plot_floes_fast(filename="out", step=1, version=2, dual=False, **kwargs):
    "Creates a video from hdf5 file output, using multiprocessing to go faster"
    data_file    = h5py.File(filename, 'r')
    nb_step = int(ceil(len(data_file.get("time")) / step))
    # Create multiple partial videos
    trunks = _get_trunks(nb_step)
    nb_process = len(trunks)
    temp_dir = 'plot_tmp'
    call(['mkdir', temp_dir])
    partial_file_names = ["{}/{}.mpg".format(temp_dir, i) for i in range(nb_process)]
    # Read Usefull data from file
    data_global = _get_useful_datas(data_file, step, version)
    # Build trunks datas
    L = [( partial_file_names[i],
           _get_useful_trunk_datas(data_global, trunk, version),
            version ) for i,trunk in enumerate(trunks)]
    # Launch process pool 
    p = Pool(nb_process)
    partial_video_maker = make_partial_floe_video_helper if not dual else make_partial_floe_video_dual_plot_helper
    p.map(partial_video_maker, L)
    p.close()
    p.join()
    # Concat all partial video
    out_filename = 'io/videos/{}.mp4'.format(filename_without_extension(filename))
    call(['ffmpeg',  '-i', 'concat:{}'.format("|".join(partial_file_names)), '-c', 'copy', out_filename])
    call(['rm', '-r', temp_dir])


def _get_useful_datas(data_file, step, version, single_step="OFF"):
    d = {}
    # File datas
    file_time_dependant_keys =["time", "floe_states", "mass_center"]
    if single_step == "OFF":
        for key in file_time_dependant_keys:
            d[key] = data_file.get(key)[::step]
    else:
        for key in file_time_dependant_keys:
            d[key] = data_file.get(key)[single_step:single_step+1]
    d["total_time"] = d["time"]
    if version == 1:
        d["floe_outlines"] = {k : dataset[::step]
            for k, dataset in data_file.get("floe_outlines").iteritems()}
    elif version == 2:
        if data_file.get("floe_shapes") is not None:
            d["floe_shapes"] = [np.array(data_file.get("floe_shapes").get(k)) for k  in sorted(list(data_file.get("floe_shapes")), key=int)]
        else:
            d["floe_shapes"] = calc_shapes(data_file)
    # Other datas
    d["impulses"] = d["total_impulses"] = calc_impulses(d["floe_states"], 6) # Calc impulsions
    # d["impulses_1"] = calc_impulses(d["floe_states"], 18) # Calc impulsions
    # d["impulses"] = [np.zeros(len(d["floe_states"][0]))] * len(d["floe_states"]) # set impulses to 0
    # calc or set global max impulse for color range
    d["MAX_IMPULSE"] = max(np.amax(step_impulses) for step_impulses in d.get("impulses"))
    # d["MAX_IMPULSE"] = 4 * 1e7
    d["WIN_WIDTH"] = getattr(OPTIONS, "window", 0) / 2 or data_file.get("window")[1] # TODO more general (ok for square centered on 0)
    # Set static axes from window data
    w = data_file.get("window")
    w_width, w_length = w[1] - w[0], w[3] - w[2]
    global STATIC_AXES # TODO something better than this global
    STATIC_AXES = [w[0] - w_width/4, w[1] + w_width/4, w[2] - w_length/4, w[3] + w_length/4]
    # STATIC_AXES = [w[0] - w_width*2, w[1] + w_width*2, w[2] - w_length*2, w[3] + w_length*2] #MPI
    # d["window"] = data_file.get("window")
    d["special_indices"] = getattr(OPTIONS, "index", [])
    return d


def _get_useful_trunk_datas(data_global, trunk, version):
    """Slice global datas for a partial video"""
    d = {}
    trunked_keys = ["time", "floe_states", "impulses", "mass_center"]
    global_keys = ["total_time","total_impulses", "MAX_IMPULSE", "WIN_WIDTH", "special_indices"]
    for key in trunked_keys:
        d[key] = data_global.get(key)[trunk[0] : trunk[1]]
    for key in global_keys:
        d[key] = data_global.get(key)
    d["total_trunk_indices"] = (trunk[0], trunk[1])
    if version == 1:
        d["floe_outlines"] = {k : dataset[trunk[0] : trunk[1]]
            for k, dataset in data_global.get("floe_outlines").iteritems()}
    elif version == 2:
        d["floe_shapes"] = data_global.get("floe_shapes")
    return d


def calc_impulses(floe_states, n=1):
    "calcul floe received impulses between each step and step -n"
    imps = [np.array([state[6] for state in time_states]) for time_states in floe_states]
    imps = [np.subtract(imps[t], imps[max(t-n, 0)]) for t in range(len(imps))]
    # for imp in imps: # hide bad values
    #     np.putmask(imp, imp>1e6, 0)
    for i in range(len(imps)):
        if np.amax(imps[i]) > 1e8:
            imps[i] = imps[i-1]
    return imps


def calc_shapes(data):
    "calculate floe shapes (in relative frame) from outline and state"
    def rotation_mat(theta):
        return np.array([[np.cos(theta), -np.sin(theta)], 
                         [np.sin(theta),  np.cos(theta)]])
    resp = [np.array(data.get("floe_outlines").get(k)[0]) for k  in sorted(list(data.get("floe_outlines")), key=int)]
    resp = [np.add(shape, np.repeat([[-x[0], -x[1]]], len(shape), axis=0)) for x, shape in zip(data.get("floe_states")[0], resp)]
    rots = [rotation_mat(-x[2]) for x in data.get("floe_states")[0]]
    resp = [np.transpose(np.dot(rot, np.transpose(shape))) for rot, shape in zip(rots, resp)]
    return resp


def _get_trunks(nb_steps):
    nb_process = min(cpu_count(), max(nb_steps / 5, 1))
    resp = []
    trunk_size = nb_steps / nb_process
    trunk_rest = nb_steps % nb_process
    j = 0
    for i in range(nb_process):
        add = 1 if i < trunk_rest else 0
        resp.append((j, j + trunk_size + add))
        j += trunk_size + add
    return resp

####################
# Running module : #
####################

import time

bcolors = {
    "HEADER" : '\033[95m',
    "OKBLUE" : '\033[94m',
    "OKGREEN" : '\033[92m',
    "WARNING" : '\033[93m',
    "FAIL" : '\033[91m',
    "ENDC" : '\033[0m',
    "BOLD" : '\033[1m',
    "UNDERLINE" : '\033[4m'
}

def colored(s, key):
    return bcolors.get(key) + s + bcolors.get("ENDC")

def timeit(func):
    """decorator to measure function execution time"""
    def timed(*args, **kw):
        ts = time.time()
        result = func(*args, **kw)
        te = time.time()
        print colored('%r %2.2f s' % (func.__name__, te-ts), "OKGREEN")
        return result

    return timed

@timeit
def run():
    global OPTIONS
    import argparse
    parser = argparse.ArgumentParser(description='Drawing floes from simulation output files')
    parser.add_argument('function', metavar='Function', type=str, help='Function to call (anim, fast_vid, 1step)')
    parser.add_argument('-f', '--file', dest="filename", help='Output filename to consider')
    parser.add_argument('--step', type=int, default=1, help='Data reading step')
    parser.add_argument('-w', '--window', type=float, default=0, help='Display window')
    parser.add_argument('--dual', dest="dual", action="store_true", default=False, help="Plot some floe's impulse (use with -i)")
    parser.add_argument('-i', "--index", action="append", type=int, help='Indexes of special floes')
    parser.add_argument('--nocolor', dest='color', action="store_false", default=True, help='Color floes according to impulses')
    parser.add_argument('--follow', action="store_true", default=False, help="Automatically adjust axes to follow floes")
    parser.add_argument('--ghosts', action="store_true", default=False, help="Display ghost floes (PBC)")
    parser.add_argument('--hd', action="store_true", default=False, help="Make HD video")
    parser.add_argument('-c', '--codec', dest="codec", default=None, help='Video codec')
    OPTIONS = parser.parse_args()
    # print vars(OPTIONS)
    funcs = {
        "anim" : plot_floes,
        "mkvid" : plot_floes_fast,
        "img" : plot_one,
    }
    if OPTIONS.codec:
        writer.codec = "lib{}".format(OPTIONS.codec)
        # print writer.extra_args
        # writer.extra_args = ['-profile:v', 'high444', '-tune:v', 'animation', '-preset:v', 'slow', '-level', '4.0', '-x264opts', 'crf=18']
        writer.extra_args = ['-profile:v', 'high444', '-tune:v', 'animation', '-preset:v', 'slow', '-level', '4.0']

    funcs[OPTIONS.function](**vars(OPTIONS))


writer = animation.FFMpegWriter(fps=30, metadata=dict(artist='Me'), extra_args=None)#, codec="libx264", bitrate=32000

# STATIC_AXES = [-800, 200,-300, 300] # xmin, xmax, ymin, ymax # 4500 vs obs dual plot version
# STATIC_AXES = [-1000, 200,-400, 400] # xmin, xmax, ymin, ymax # 4500 vs obs
STATIC_AXES = [-1700, 1700, -1700, 1700] # xmin, xmax, ymin, ymax
# STATIC_AXES = False # if False, axes will automatically be adjusted to see all floes
OPTIONS = None


if __name__ == "__main__":
    run()

