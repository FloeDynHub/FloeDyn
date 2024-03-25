#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import cm, animation, colors, gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PolyCollection

from multiprocessing import Pool, cpu_count
from subprocess import call

from utils import filename_without_extension, get_unused_path, check_path_existence, mkdir_path
import h5py
import datetime
import math
import os


class AxeManager(object):
    """
    Made to allow named patches and collections in matplotlib Axes objects
    """

    def __init__(self, ax):
        self.ax = ax
        self.collections = {}
        self.patch_lists = {}
        self.patches = {}

    def set_collection(self, key, obj):
        idx = len(self.ax.collections)
        self.ax.add_collection(obj)
        self.collections[key] = idx

    def set_patch(self, key, obj):
        idx = len(self.ax.collections)
        self.ax.add_patch(obj)
        self.patches[key] = idx

    def add_to_patch_list(self, key, obj):
        idx = len(self.ax.patches)
        self.ax.add_patch(obj)
        if key in self.patch_lists:
            self.patch_lists[key].append(idx)
        else:
            self.patch_lists[key] = [idx]

    def get_patch_list(self, key):
        return [self.ax.patches[i] for i in self.patch_lists[key]]

    def get_patch(self, key):
        return self.ax.patches[self.patches[key]]

    def get_collection(self, key):
        return self.ax.collections[self.collections[key]]

    def replace_collection(self, key, obj):
        idx = self.collections[key]
        del self.ax.collections[idx]
        self.set_collection(key, obj)


class FloePlotter(object):
    """ Drawing floes from simulation output files """

    def __init__(self, options):
        self.OPTIONS = options
        self.writer = animation.FFMpegWriter(fps=30, metadata=dict(artist='Me'), extra_args=None)#, codec="libx264", bitrate=32000
        self.colors = {
            "ocean" : '#162252',
            # "ocean" : 'white',
            "circles": 'white',
            "window": 'white'
        }

    def do(self):
        funcs = {
            "anim" : self.plot_floes,
            "mkvid" : self.plot_floes_fast,
            "img" : self.plot_one,
        }
        if self.OPTIONS.codec:
            self.writer.codec = "lib{}".format(self.OPTIONS.codec)
            # writer.extra_args = ['-profile:v', 'high444', '-tune:v', 'animation', '-preset:v', 'slow', '-level', '4.0', '-x264opts', 'crf=18']
            self.writer.extra_args = ['-profile:v', 'high444', '-tune:v', 'animation', '-preset:v', 'slow', '-level', '4.0']
        funcs[self.OPTIONS.function]()

    def get_final_video_path(self, video_ext="mp4"):
        filename = self.OPTIONS.outname or filename_without_extension(self.OPTIONS.filename)
        if not check_path_existence("io/videos"):
            print('Warning: the directory ''io/videos'' does not exist! It will be created')
            mkdir_path('io/videos')
        return get_unused_path('io/videos/{}.{}'.format(filename, video_ext))

    def _init1(self, data, ax):
        "initializes floes (Polygon creation) for matplotlib animation creation"
        for (i,s) in data.get("floe_outlines").items():
            special = i=="1500"
            polygon = Polygon(s[0],
                              True,
                              facecolor="#FFF7ED" if not special else "#FA9A50",
                              linewidth=0.2,
                              hatch="" if not special else "x")
            ax.add_patch(polygon)
        return ax

    def _update1(self, num, data, ax):
        "updates floes (Polygons) positions for matplotlib animation creation"
        begin, step = 0, 1
        indic = begin + step * num
        for i, s in enumerate(data.get("floe_outlines").values()):
            ax.patches[i].set_xy(s[indic])
        ax.set_title("t = {}".format(str(datetime.timedelta(seconds=int(data.get("time")[indic])))))
        if not data.get("static_axes"):
            ax.axis('equal')
            ax.relim()
            ax.autoscale_view(True,True,True)
        return ax,

    def _display_window(self, data, ax_mgr):
        # Display generator window
        w = data.get("window")
        if w:
            ax_mgr.set_patch("window",Polygon(
                ((w[0], w[2]), (w[0], w[3]), (w[1], w[3]), (w[1], w[2])),
                # True, facecolor=None, alpha=0.2, edgecolor="white", linewidth=0.2))
                True, facecolor="none", edgecolor=self.colors["window"], linewidth=1.5))
    
    def _display_mesh(self, data, ax_mgr):
        # Display mesh. Vous noterez au passage l'utilite de ce commentaire. 
        coord = data.get("floe_meshes_coord")
        connect = data.get("floe_meshes_connect")
        displayedData = []
        if (coord and connect):
            ax = ax_mgr.ax
            count = 0
            # building a list of np.array containing the node coordinates of each element 
            elements = []
            for iFloe in np.arange(len(coord)):
                for iElem in np.arange(len(connect[iFloe])):
                    elements.append(np.array([coord[iFloe][int(connect[iFloe][iElem][0])], coord[iFloe][int(connect[iFloe][iElem][1])], coord[iFloe][int(connect[iFloe][iElem][2])]]))
                    displayedData.append(0)
            meshes_collec = PolyCollection(
                elements,
                linewidths=0.5,
                edgecolors="black",
                cmap=cm.YlOrRd, 
                # norm=colors.PowerNorm(gamma=0.3)
                )
            meshes_collec.set_array(displayedData)
            meshes_collec.set_clim([0, max(displayedData)])
            plt.colorbar(meshes_collec, ax = ax, fraction = 0.05)
            ax_mgr.set_collection("meshes", meshes_collec)


    def _init_floes(self, data, ax_mgr):
        "initializes floes (Polygon creation) for matplotlib animation creation (v3 : shapes)"
        ax = ax_mgr.ax
        ax.set_facecolor(self.colors["ocean"])
        floes_collec = PolyCollection(
            data.get("floe_shapes"),
            # linewidths=0.2,
            linewidths=0.05,
            # edgecolors="red", facecolors="none") # no facecolor mode
            cmap=cm.YlOrRd, norm=colors.PowerNorm(gamma=0.3))
            # cmap=cm.Paired)#MPI
        if getattr(self.OPTIONS, "color", True):
        # if getattr(self.OPTIONS, "color", True) and getattr(self.OPTIONS, "colorby", 'None'):
            floes_collec.set_array(data.get("impulses")[0])
            floes_collec.set_clim([0, data["MAX_IMPULSE"]])
            plt.colorbar(floes_collec, ax = ax, fraction = 0.05) # display color bar
        else:
            floes_collec.set_array(np.zeros(len(data.get("floe_shapes")))) # comment for no facecolor mode
        # ax.axes.get_yaxis().set_visible(False) # remove x axis informations

        # if --disp_mesh is asked, we first check that a mesh is actually given in the output file         
        opt_mesh = (getattr(self.OPTIONS, "disp_mesh"))
        if data.get("floe_meshes_connect") is None:
            if opt_mesh:
                print('No mesh found in the output file, unable to plot it.')
            opt_mesh = False 
        # the floe shapes are plotted only if no mesh is shown. 
        if not opt_mesh:
            ax_mgr.set_collection("floes", floes_collec)

        if getattr(self.OPTIONS, "ghosts", False):
            ghosts_collec = PolyCollection(data.get("floe_shapes") * 8, linewidths=0.2, alpha=0.4,
                                         cmap=cm.YlOrRd, norm=colors.PowerNorm(gamma=0.3))
                                        # cmap=cm.Paired)#MPI
            ghosts_collec.set_clim(floes_collec.get_clim())
            ghosts_collec.set_array(np.zeros(len(data.get("floe_shapes")) * 8))
            ax_mgr.set_collection("floe_ghosts", ghosts_collec)

    def _init_circles(self, data, ax_mgr, begin=0):
        # floes_collec = CircleCollection(np.zeros(shape=(len(data.get(floe_shapes)), 2),
        #                                 linewidths=0.2, facecolors=None, edgecolors="white"))
        # ax.add_collection(floes_collec)
        ax = ax_mgr.ax
        for shape in data.get("floe_shapes"):
            radius = max([math.sqrt(x**2+y**2) for x,y in shape])
            # circle = Circle((0,0), radius, facecolor='none', edgecolor="white", linewidth=2)
            circle = Circle((0,0), radius,
                facecolor='none', edgecolor=self.colors["circles"],
                linewidth=1
            )
            ax_mgr.add_to_patch_list("floe_bounding_circles", circle)
        # return ax

    def _init2(self, data, ax_mgr):
        if getattr(self.OPTIONS, "disp_floes"):
            self._init_floes(data, ax_mgr)
        if getattr(self.OPTIONS, "disp_circles"):
            self._init_circles(data, ax_mgr)
        if getattr(self.OPTIONS, "disp_window"):
            self._display_window(data, ax_mgr)
        if getattr(self.OPTIONS, "disp_mesh"):
            self._display_mesh(data, ax_mgr)
        # return ax

    def _update_floes(self, num, data, ax_mgr):
        """updates floes (Polygons) positions for matplotlib animation creation
        (version 3 : collection and compute outlines transforms"""
        ax = ax_mgr.ax
        begin, step = 0, 1
        indic = begin + step * num
        opt_ghosts, opt_color, opt_follow, opt_fracture, opt_mesh = (
            getattr(self.OPTIONS, opt) for opt in ["ghosts", "color", "follow", "fracture", "disp_mesh"]
        )
        if data.get("floe_meshes_connect") is None:
            # if opt_mesh:
            #     print('No mesh found in the output file, unable to plot it.') # should already have been displayed in _init_floes  
            opt_mesh = False 
        # if opt_colorby != 'None': 
        #     avail_data_list = ['Ux', 'Uy', 'U', 'index']
        #     avail_data_list = ['index', 'Ux', 'Uy', 'U']
        #     # this list should match the one in void HDF5ManagerWithMesh<TFloeGroup, TDynamicsMgr>::save_step(real_type time, const dynamics_mgr_type& dynamics_manager)
        #     if (opt_colorby not in avail_data_list):
        #         if indic == 1:
        #             print("No {} color data available".format(opt_colorby))
        #             print("please use one of {}".format(avail_data_list))
        #         opt_colorby = 'None'
        #     else:
        #         opt_mesh = True
        #         color_data_id = avail_data_list.index(opt_colorby)

        nb_floes = len(data.get("floe_shapes"))
        verts = self._transform_shapes(data.get("floe_shapes"), data.get("floe_states")[indic], opt_follow) # [0:min(num,nb_floes)]

        # display a warning if there are disabled floes in the database but no -c / --crack option is provided (not critical, but a bunch of weird floes may appear as well)
        if (not opt_fracture) and (indic == 1) and (abs(np.sum(data.get("floe_states")[indic,:,9]) - len(data.get("floe_states")[indic,:,9])) > 0.1):
            print('There are {} disabled floes at the beginning, did you forget the --crack option?'.format(int(len(data.get("floe_states")[indic,:,9]) - np.sum(data.get("floe_states")[indic,:,9]))))

        if not opt_mesh:
            # displaying floe shapes
            if opt_fracture:
                # Filter floe_shapes according to floe_state "alive" because nb of floes may have change
                verts = [v for i, v in enumerate(verts) if data.get("floe_states")[indic][i][9] == 1]
            ax_mgr.get_collection("floes").set_verts(verts)
            # coloring floe shapes
            if opt_color:
                ax_mgr.get_collection("floes").set_array(data.get("impulses")[indic])
                if opt_ghosts:
                    ax_mgr.get_collection("floe_ghosts").set_verts(verts)
        else:
            # displaying floe meshes and colors when appropriate
            connect = data.get("floe_meshes_connect")
            coord = self._transform_meshes(data.get("floe_meshes_coord"), data.get("floe_states")[indic], data.get("mass_center")[indic], opt_follow)
            if opt_fracture:
                # filtering disabled floes 
                connect = [v for i, v in enumerate(connect) if data.get("floe_states")[indic][i][9] == 1]
                coord = [v for i, v in enumerate(coord) if data.get("floe_states")[indic][i][9] == 1]            
            elements = []

            displayedData = []
            if data.get("floe_elem_data") is not None: 
                color_data = data.get("floe_elem_data")
                color_data_indic = indic 
                cdTemp = np.zeros(np.array(color_data).shape)
                k = 0
                for iFloe in range(0, nb_floes):
                    if data.get("floe_states")[indic][iFloe][9] == 1:
                        cdTemp[:,k,:] = color_data[:,iFloe,:]
                        k = k+1
                color_data = cdTemp
            else:
                nElem = max([len(c) for i,c in enumerate(connect)])
                color_data = np.zeros([1, len(coord), 1, nElem])
                # color_data_id=0
                color_data_indic=0


            for iFloe in np.arange(len(coord)):
                for iElem in np.arange(len(connect[iFloe])):
                    elements.append(np.array([coord[iFloe][int(connect[iFloe][iElem][0])], coord[iFloe][int(connect[iFloe][iElem][1])], coord[iFloe][int(connect[iFloe][iElem][2])]]))
                    displayedData.append(color_data[color_data_indic, iFloe, iElem])
                    # displayedData.append(color_data[color_data_indic, iFloe, iElem, color_data_id])
                    # displayedData.append(color_data[color_data_indic, iFloe, color_data_id, iElem])

            globalMin = np.min(np.min(np.min(color_data[:, :, :])))
            globalMax = np.max(np.max(np.max(color_data[:, :, :])))
            
            ax_mgr.get_collection("meshes").set_verts(elements)
            ax_mgr.get_collection("meshes").set_clim([globalMin, globalMax])
            ax_mgr.get_collection("meshes").set_array(displayedData)

        if opt_ghosts:
            # attempt to fix bug: Matthias
            w = data.get("window")
            w_width, w_length = w[1] - w[0], w[3] - w[2]
            ghosts_trans = [(i * w_width, j * w_length) for i in range(-1, 2) for j in range(-1, 2) if not i==j==0]
            ghosts_verts = [v for trans in ghosts_trans for v in self.translate_group(verts, trans)]
            # original line: ghosts_verts = [v for trans in d["ghosts_trans"] for v in translate_group(verts, trans)] 
            # End of attempt
            ax_mgr.get_collection("floe_ghosts").set_verts(ghosts_verts)

        if opt_color:
            impulses = [v for i, v in enumerate(data.get("impulses")[indic]) if data.get("floe_states")[indic][i][9] == 1]
            ax_mgr.get_collection("floes").set_array(impulses)
            if opt_ghosts:
                ax_mgr.get_collection("floe_ghosts").set_array(np.tile(impulses, 8))
        ax.set_title("t = {}".format(str(datetime.timedelta(seconds=int(data.get("time")[indic])))))
        if not data.get("static_axes"):
            # ax.axis('equal')
            ax.relim()
            ax.autoscale_view(True,True,True)
        if opt_follow:
            axes = data.get("static_axes")
            ax.set_xlim(data.get("mass_center")[indic][0] + axes[0], data.get("mass_center")[indic][0] +  axes[1])
            ax.set_ylim(data.get("mass_center")[indic][1] + axes[2], data.get("mass_center")[indic][1] +  axes[3])
        # # begin MPI zones
        # ax.collections[0].set_clim([0, 9])
        # ax.collections[1].set_clim([0, 9])
        # W, H = ax.get_xlim()[1] - ax.get_xlim()[0], ax.get_ylim()[1] - ax.get_ylim()[0]
        # N = 6
        # WZ, HZ = W/N, H/N
        # ax.collections[0].set_array(np.array([(int((x[0] - ax.get_xlim()[0]) / WZ) + int((x[1] - ax.get_ylim()[0]) / HZ) * N)%9 for x in data.get("floe_states")[indic]]))
        # ghosts_trans = [(i * w, j * w) for i in range(-2, 3) for j in range(-2, 3) if not i==j==0]
        # gstates = [(x[0] + i, x[1] + j) for i,j in ghosts_trans for x in data.get("floe_states")[indic]]
        # ax.collections[1].set_array(np.array([(int((x[0] - ax.get_xlim()[0]) / WZ) + int((x[1] - ax.get_ylim()[0]) / HZ) * N)%9 for x in gstates]))
        # # end MPI zones
        # return ax,
    

    def _update_circles(self, num, data, ax_mgr):
        # ax = ax_mgr.ax
        begin, step = 0, 1
        indic = begin + step * num
        for state, circle in zip(data.get("floe_states")[indic], ax_mgr.get_patch_list("floe_bounding_circles")):
            circle.center = state[0], state[1]
        # return ax,

    def _update2(self, num, data, ax_mgr):
        """updates floes (Polygons) positions for matplotlib animation creation
        (version 3 : collection and compute outlines transforms"""
        if getattr(self.OPTIONS, "disp_floes"):
            self._update_floes(num, data, ax_mgr)
        if getattr(self.OPTIONS, "disp_circles"):
            self._update_circles(num, data, ax_mgr)
        # if getattr(self.OPTIONS, "disp_mesh"):
        #     self._update_meshes(num, data, ax_mgr)
        

    def _init2_dual(self, data, ax1, ax2):
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

    def _update2_dual(self, num, data, ax1, ax2):
        first_idx, _ = data["total_trunk_indices"]
        indic = first_idx + num
        ax1, = self._update2(num, data, ax1)
        for i, idx in enumerate(data.get("special_indices")):
            ax2.lines[i].set_xdata(np.append(ax2.lines[i].get_xdata(), data.get("total_time")[indic]))
            ax2.lines[i].set_ydata(np.append(ax2.lines[i].get_ydata(), data.get("total_impulses")[indic][idx]))
        return ax1, ax2

    def _transform_shapes(self, floe_shapes, floe_states, follow=False):
        resp = floe_shapes
        def rotation_mat(theta):
            return np.array([[np.cos(theta), -np.sin(theta)], 
                             [np.sin(theta),  np.cos(theta)]])
        rots = [rotation_mat(x[2]) for x in floe_states]
        resp = [np.transpose(np.dot(rot, np.transpose(shape))) for rot, shape in zip(rots, floe_shapes)]
        pos_ids = (0,1) if follow else (7,8) # (7,8) contains translated position in fixed initial window
        resp = [np.add(shape, np.repeat([[x[pos_ids[0]], x[pos_ids[1]]]], len(shape), axis=0)) for x, shape in zip(floe_states, resp)]
        return resp
    
    def _transform_meshes(self, floe_coords, floe_states, floe_mass_centers, follow=False):
        pos_ids = (0,1) if follow else (7,8) # (7,8) contains translated position in fixed initial window
        coord = floe_coords
        def rotation_mat(theta):
            return np.array([[np.cos(theta), -np.sin(theta)], 
                             [np.sin(theta),  np.cos(theta)]])
        rots = [rotation_mat(x[2]) for x in floe_states]
        coord = [np.transpose(np.dot(rot, np.transpose(mesh))) for rot, mesh in zip(rots, floe_coords)]
        coord = [np.add(mesh, np.repeat([[x[pos_ids[0]], x[pos_ids[1]]]], len(mesh), axis=0)) for x, mesh in zip(floe_states, coord)]
        return coord

    def translate_group(self, vertices, trans):
        return [np.add(shape, np.repeat([[trans[0], trans[1]]], len(shape), axis=0)) for shape in vertices]

    def get_anim_fcts(self, version=None):
        d = {
            1: (self._init1,self._update1),
            2: (self._init2,self._update2),
            3: (self._init_circles, self._update_circles) # obsolete with -v 2 --circles --nofloes
        }
        return d[version or self.OPTIONS.version]

    def plot_floes(self, make_video=1):
        "displays floes animation from hdf5 output file (simple plot or video creation)"
        data_file    = h5py.File(self.OPTIONS.filename, 'r')
        # Read Usefull data from file
        data_global = self._get_useful_data(data_file)
        fig, ax = plt.subplots()
        ax_mgr = AxeManager(ax)
        if getattr(self.OPTIONS, "hd", False):
            fig.set_size_inches(20, 15)
            fig.set_dpi(100)
        init, update = self.get_anim_fcts()
        init(data_global, ax_mgr)
        # print ax_mgr.patches
        if not data_global.get("static_axes"):
            ax.axis('equal') # automatic scale
        else:
            ax.axis(data_global.get("static_axes"))
        anim = animation.FuncAnimation(
            fig, update, len(data_global.get("time")), fargs=(data_global, ax_mgr),
            interval=1, blit=False)
        if make_video:
            anim.save(self.get_final_video_path(), writer=self.writer)
        else:
            plt.show()

    def plot_one(self):
        "displays floes at last time recorded in output file"
        init, update = self.get_anim_fcts()
        file    = h5py.File(self.OPTIONS.filename, 'r')
        fig, ax = plt.subplots()
        ax_mgr = AxeManager(ax)
        if file.get("time") and len(file.get("time")):
            idx = int(input("Time index (-1 for last one) : "))
        else: # ex : input files
            idx = 0
        num = idx

        data_global = self._get_useful_data(file, num, img=True)
        init(data_global, ax_mgr)
        if not data_global.get("static_axes"):
            ax.axis('equal') # automatic scale
        else:
            ax.axis(data_global.get("static_axes"))
        update(0, data_global, ax_mgr)
        if getattr(self.OPTIONS, "hd", False):
            fig.set_size_inches(60, 45)
            fig.set_dpi(100)
            plt.savefig('{}.png'.format(os.path.splitext(self.OPTIONS.filename)[0]))
            # plt.savefig('{}.eps'.format(self.OPTIONS.filename), format='eps')#, dpi=1000)
        else:
            plt.show()


    #############################
    # Parallel video creation : #
    #############################

    def make_partial_floe_video(self, out_filename, data_chunk, version):
        "creates video directly from data"
        # print(1)
        fig, ax = plt.subplots()
        # print(2)
        if getattr(self.OPTIONS, "hd", False):
            fig.set_size_inches(20, 15)
            fig.set_dpi(100)
        # print(3)
        init, update = self.get_anim_fcts(version)
        ax_mgr = AxeManager(ax)
        init(data_chunk, ax_mgr)
        # print(4)
        if not data_chunk.get("static_axes"):
            ax.axis('equal') # automatic scale
        else:
            ax.axis(data_chunk.get("static_axes"))
        # print(5)
        # ax.set_axis_bgcolor('#162252')
        anim = animation.FuncAnimation(
            fig, update, len(data_chunk.get("time")), fargs=(data_chunk, ax_mgr),
            interval=1, blit=False)
        # print(6)
        print(out_filename)
        anim.save(out_filename, writer=self.writer)

    def make_partial_floe_video_helper(self, t):
        return self.make_partial_floe_video(*t)


    def make_partial_floe_video_dual_plot(self, out_filename, data_chunk, version):
        "creates video directly from data, and adds obstacle impulse subplot under floes"
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        fig = plt.figure()
        if getattr(self.OPTIONS, "hd", False):
            fig.set_size_inches(16, 12)
            fig.set_dpi(100)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        init, update = _init2_dual, _update2_dual
        ax1, ax2 = init(data_chunk, ax1, ax2)
        if not data_chunk.get("static_axes"):
            ax1.axis('equal') # automatic scale
        else:
            ax1.axis(data_chunk.get("static_axes"))
        # ax1.set_axis_bgcolor('#162252')
        anim = animation.FuncAnimation(
            fig, update, len(data_chunk.get("time")), fargs=(data_chunk, ax1, ax2),
            interval=50, blit=False)
        anim.save(out_filename, writer=self.writer)

    def make_partial_floe_video_dual_plot_helper(self, t):
        return self.make_partial_floe_video_dual_plot(*t)

    def plot_floes_fast(self, dual=False, **kwargs):
        "Creates a video from hdf5 file output, using multiprocessing to go faster"
        data_file    = h5py.File(self.OPTIONS.filename, 'r')
        nb_step = int(math.ceil(len(data_file.get("time")) // self.OPTIONS.step))
        # Create multiple partial videos
        trunks = self._get_trunks(nb_step)
        nb_process = len(trunks)
        temp_dir = 'plot_tmp'
        call(['mkdir', temp_dir])
        partial_file_names = ["{}/{}.mpg".format(temp_dir, i) for i in range(nb_process)]
        # Read Usefull data from file
        data_global = self._get_useful_data(data_file)
        # Build trunks datas
        L = [( partial_file_names[i],
               self._get_useful_trunk_data(data_global, trunk),
                self.OPTIONS.version ) for i,trunk in enumerate(trunks)]
        # Launch process pool 
        p = Pool(nb_process)
        partial_video_maker = self.make_partial_floe_video_helper if not dual else make_partial_floe_video_dual_plot_helper
        p.map(partial_video_maker, L)
        p.close()
        p.join()
        # Concat all partial video
        out_filename = self.get_final_video_path()
        call(['ffmpeg',  '-i', 'concat:{}'.format("|".join(partial_file_names)), '-c', 'copy', out_filename])
        call(['rm', '-r', temp_dir])


    def _get_useful_data(self, data_file, single_step="OFF", img=False):
        d = {}
        # File datas
        if not img:
            file_time_dependant_keys =["time", "floe_states", "mass_center", "floe_elem_data"]
        else:
            file_time_dependant_keys =["time", "floe_states"] # allow input files to be plotted
        if single_step == "OFF":
            for key in file_time_dependant_keys:
                if data_file.get(key) is not None:
                    d[key] = data_file.get(key)[::self.OPTIONS.step]
        else:
            for key in file_time_dependant_keys:
                if img and key == "time" and not "time" in data_file: # plot input files
                    d[key] = [0]
                else:
                    d[key] = [data_file.get(key)[single_step]]
        if not img:
            d["total_time"] = d["time"]
        if self.OPTIONS.version == 1:
            d["floe_outlines"] = {k : dataset[::self.OPTIONS.step]
                for k, dataset in data_file.get("floe_outlines").items()}
        elif self.OPTIONS.version >= 2:
            if data_file.get("floe_shapes") is not None:
                d["floe_shapes"] = [np.array(data_file.get("floe_shapes").get(k)) for k in sorted(list(data_file.get("floe_shapes")), key=int)]
                if data_file.get("floe_meshes_coord") is not None:
                    d["floe_meshes_coord"] = [np.array(data_file.get("floe_meshes_coord").get(k)) for k in sorted(list(data_file.get("floe_meshes_coord")), key=int)]
                if data_file.get("floe_meshes_connect") is not None:
                    d["floe_meshes_connect"] = [np.array(data_file.get("floe_meshes_connect").get(k)) for k in sorted(list(data_file.get("floe_meshes_connect")), key=int)]
            else:
                d["floe_shapes"] = self.calc_shapes(data_file)
        # Other datas
        d["impulses"] = d["total_impulses"] = self.calc_impulses(d["floe_states"], 6) # Calc impulsions
        # d["impulses"] = [np.zeros(len(d["floe_states"][0]))] * len(d["floe_states"]) # set impulses to 0
        # calc or set global max impulse for color range
        d["MAX_IMPULSE"] = max(np.amax(step_impulses) for step_impulses in d.get("impulses"))
        # Set static axes from window data
        w = data_file.get("window")
        w_width, w_length = w[1] - w[0], w[3] - w[2]
        d["static_axes"] = self.OPTIONS.static_axes
        if not d["static_axes"]:
            d["static_axes"] = [w[0] - w_width/2, w[1] + w_width/2, w[2] - w_length/2, w[3] + w_length/2]
            # d["static_axes"] = [w[0] - w_width*2, w[1] + w_width*2, w[2] - w_length*2, w[3] + w_length*2] #MPI bigger domain
        d["window"] = list(data_file.get("window", None))
        if getattr(self.OPTIONS, "ghosts"):
            d["ghosts_trans"] = [(i * w_width, j * w_length) for i in range(-1, 2) for j in range(-1, 2) if not i==j==0]
            # d["ghosts_trans"] = [(i * w, j * w) for i in range(-2, 3) for j in range(-2, 3) if not i==j==0] #MPI, more ghosts
        d["special_indices"] = getattr(self.OPTIONS, "index", [])
        return d


    def _get_useful_trunk_data(self, data_global, trunk):
        """Slice global datas for a partial video"""
        d = {}
        trunked_keys = ["time", "floe_states", "impulses", "mass_center"]
        global_keys = [
            "total_time","total_impulses", "MAX_IMPULSE", "window",
            "special_indices", "static_axes", "ghosts_trans"
        ]
        print(trunk)
        for key in trunked_keys:
            d[key] = data_global.get(key)[trunk[0] : trunk[1]]
        for key in global_keys:
            d[key] = data_global.get(key)
        d["total_trunk_indices"] = (trunk[0], trunk[1])
        if self.OPTIONS.version == 1:
            d["floe_outlines"] = {k : dataset[trunk[0] : trunk[1]]
                for k, dataset in data_global.get("floe_outlines").items()}
        elif self.OPTIONS.version >= 2:
            d["floe_shapes"] = data_global.get("floe_shapes")
        return d


    def calc_impulses(self, floe_states, n=1):
        "calcul floe received impulses between each step and step -n"
        imps = [np.array([state[6] for state in time_states]) for time_states in floe_states]
        imps = [np.subtract(imps[t], imps[max(t-n, 0)]) for t in range(len(imps))]
        # for imp in imps: # hide bad values
        #     np.putmask(imp, imp>1e6, 0)
        # for i in range(len(imps)):
        #     if np.amax(imps[i]) > 1e8:
        #         imps[i] = imps[i-1]
        return imps


    def calc_shapes(self, data):
        "calculate floe shapes (in relative frame) from outline and state"
        def rotation_mat(self, theta):
            return np.array([[np.cos(theta), -np.sin(theta)], 
                             [np.sin(theta),  np.cos(theta)]])
        resp = [np.array(data.get("floe_outlines").get(k)[0]) for k  in sorted(list(data.get("floe_outlines")), key=int)]
        resp = [np.add(shape, np.repeat([[-x[0], -x[1]]], len(shape), axis=0)) for x, shape in zip(data.get("floe_states")[0], resp)]
        rots = [rotation_mat(-x[2]) for x in data.get("floe_states")[0]]
        resp = [np.transpose(np.dot(rot, np.transpose(shape))) for rot, shape in zip(rots, resp)]
        return resp


    def _get_trunks(self, nb_steps):
        nb_process = min(cpu_count(), max(nb_steps // 5, 1))
        resp = []
        trunk_size = nb_steps // nb_process
        trunk_rest = nb_steps % nb_process
        j = 0
        for i in range(nb_process):
            add = 1 if i < trunk_rest else 0
            resp.append((j, j + trunk_size + add))
            j += trunk_size + add
        return resp
