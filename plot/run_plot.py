#!/usr/bin/env python
# -*- coding: utf-8 -*-
from utils import timeit
import argparse
from plotter import FloePlotter


@timeit
def run():
    parser = argparse.ArgumentParser(description='Drawing floes from simulation output files')
    parser.add_argument('function', metavar='Function', type=str, help='Function to call (anim, fast_vid, 1step)')
    parser.add_argument('-v,', '--version', type=int, default=2, help='init/update version') # TODO improve this design
    parser.add_argument('-f', '--file', dest="filename", help='Output filename to consider')
    parser.add_argument('--step', type=int, default=1, help='Data reading step')
    parser.add_argument('-w', '--window', action="store_true", dest="disp_window", default=False, help='Display window')
    parser.add_argument('--circles', action="store_true", dest="disp_circles", default=False, help='Display floes bounding circles')
    parser.add_argument('--nofloes', dest='disp_floes', action="store_false", default=True, help='Do not display floes')
    parser.add_argument('--dual', dest="dual", action="store_true", default=False, help="Plot some floe's impulse (use with -i)")
    parser.add_argument('-i', "--index", action="append", type=int, help='Indexes of special floes')
    parser.add_argument('--nocolor', dest='color', action="store_false", default=True, help='Color floes according to impulses')
    parser.add_argument('--follow', action="store_true", default=False, help="Automatically adjust axes to follow floes")
    parser.add_argument('--ghosts', action="store_true", default=False, help="Display ghost floes (PBC)")
    parser.add_argument('--crack', dest="fracture", action="store_true", default=False, help="Floes can fracture")
    parser.add_argument('--hd', action="store_true", default=False, help="Make HD video")
    parser.add_argument('-c', '--codec', dest="codec", default=None, help='Video codec')
    parser.add_argument('-a', '--axes', dest="static_axes", type=str, help='Initial Axes (xmin,xmax,ymin,ymax)', default=None)
    parser.add_argument('--mesh', '-m', action="store_true", dest="disp_mesh", default=False, help='Display mesh')
    # parser.add_argument('--colorby', '-C', dest="colorby", type=str, default='None', help='Displays mesh colored by specified data')
    OPTIONS = parser.parse_args()
    if OPTIONS.static_axes:
        try:
            OPTIONS.static_axes = [float(x) for x in OPTIONS.static_axes.split(",")]
            OPTIONS.static_axes[3]
        except ValueError as e:
            raise Exception("Axes : valeurs incorrectes ({})".format(e))
        except IndexError as e:
            raise Exception("Axes : valeurs manquantes ({})".format(e))

    plotter = FloePlotter(OPTIONS)
    plotter.do()


if __name__ == "__main__":
    run()