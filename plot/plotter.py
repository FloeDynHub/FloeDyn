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

    def _init_floes(self, data, ax_mgr):
        "initializes floes (Polygon creation) for matplotlib animation creation (v3 : shapes)"
        ax = ax_mgr.ax
        ax.set_facecolor(self.colors["ocean"])
        floes_collec = PolyCollection(
            data.get("floe_shapes"),
            linewidths=0.05,
            cmap=self.get_cmap(data),
            norm=colors.PowerNorm(gamma=0.3))
            # cmap=cm.Paired) # MPI visu
        if getattr(self.OPTIONS, "color", True):
            floes_collec.set_array(data.get(data["color_key"])[0])
            floes_collec.set_clim([0, data.get(data["max_color_key"])])
            plt.colorbar(floes_collec, ax = ax, fraction = 0.05) # display color bar
        else:
            floes_collec.set_array(np.zeros(len(data.get("floe_shapes")))) # comment for no facecolor mode
        # ax.axes.get_yaxis().set_visible(False) # remove x axis informations
        ax_mgr.set_collection("floes", floes_collec)

        if getattr(self.OPTIONS, "ghosts", False):
            ghosts_collec = PolyCollection(data.get("floe_shapes") * 8, linewidths=0.2, alpha=0.4,
                                         cmap=self.get_cmap(data), norm=colors.PowerNorm(gamma=0.3))
                                        # cmap=cm.Paired)#MPI
            ghosts_collec.set_clim(floes_collec.get_clim())
            ghosts_collec.set_array(np.zeros(len(data.get("floe_shapes")) * 8))
            ax_mgr.set_collection("floe_ghosts", ghosts_collec)
    
    def get_cmap(self, data):
        return cm.YlOrRd if data["color_key"] == "impulses" else cm.GnBu

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

        # return ax

    def _update_floes(self, num, data, ax_mgr):
        """updates floes (Polygons) positions for matplotlib animation creation
        (version 3 : collection and compute outlines transforms"""
        ax = ax_mgr.ax
        begin, step = 0, 1
        indic = begin + step * num
        opt_ghosts, opt_color, opt_follow, opt_fracture = (
            getattr(self.OPTIONS, opt) for opt in ["ghosts", "color", "follow", "fracture"]
        )

        nb_floes = len(data.get("floe_shapes"))
        verts = self._transform_shapes(data.get("floe_shapes"), data.get("floe_states")[indic], opt_follow) # [0:min(num,nb_floes)]
        if opt_fracture:
            # Filter floe_shapes according to floe_state "alive" because nb of floes may have change
            verts = [v for i, v in enumerate(verts) if data.get("floe_states")[indic][i][9] == 1]
        ax_mgr.get_collection("floes").set_verts(verts)

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
            impulses = [v for i, v in enumerate(data.get(data["color_key"])[indic]) if len(data.get("floe_states")[indic][i]) < 10 or data.get("floe_states")[indic][i][9] == 1]
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
        if num > 0 and num % 10 == 0:
            progress = num / len(data.get("time"))
            print("progress::{:.2f}%".format(progress * 100))
        if getattr(self.OPTIONS, "disp_floes"):
            self._update_floes(num, data, ax_mgr)
        if getattr(self.OPTIONS, "disp_circles"):
            self._update_circles(num, data, ax_mgr)

    def _init2_dual(self, data, ax_mgr1, ax_mgr2):
        "initializes floes (Polygon creation) for matplotlib animation creation (v3 : shapes)"
        self._init2(data, ax_mgr1)
        self.init_chart(data, ax_mgr2) if not self.OPTIONS.bar else self.init_bar(data, ax_mgr2)
        return ax_mgr1.ax, ax_mgr2.ax
    
    def init_chart(self, data, ax_mgr2):
        first_idx, _ = data["total_trunk_indices"]
        x = np.array(data.get("total_time")[0:first_idx])
        ys = [np.array([data.get("total_impulses")[i][idx] for i in range(len(x))]) for idx in data.get("special_indices", [])]
        ax2 = ax_mgr2.ax
        for y in ys:
            ax2.plot(x,y)
        ax2.axis([data.get("total_time")[0], data.get("total_time")[-1], -100, data["MAX_IMPULSE"]])
        ax2.set_title("Norm of the contact impulses applied to the obstacle")

    def init_bar(self, data, ax_mgr2):
        first_idx, _ = data["total_trunk_indices"]
        special_indices = data.get("special_indices")
        impulse_key = "total_impulses" if not self.OPTIONS.cumul else "cumul_impulses"
        max_key = "MAX_IMPULSE" if not self.OPTIONS.dual else "max_received_impulse_special_indices"
        initial_heights = [data.get(impulse_key)[first_idx][i] for i in special_indices]
        ax2 = ax_mgr2.ax
        self.bars = ax2.bar(range(len(special_indices)), initial_heights)
        ax2.set_xticks(range(len(special_indices)))
        ax2.set_xlim(-1, len(special_indices))
        ax2.set_ylim(-100, data[max_key])
        return self.bars

    def _update2_dual(self, num, data, ax_mgr1, ax_mgr2):
        first_idx, _ = data["total_trunk_indices"]
        indic = first_idx + num
        self._update2(num, data, ax_mgr1)
        self.update_chart(num, data, ax_mgr2) if not self.OPTIONS.bar else self.update_bar(num, data, ax_mgr2)
            
    def update_chart(self, num, data, ax_mgr2):
        first_idx, _ = data["total_trunk_indices"]
        indic = first_idx + num
        ax2 = ax_mgr2.ax
        for i, idx in enumerate(data.get("special_indices")):
            ax2.lines[i].set_xdata(np.append(ax2.lines[i].get_xdata(), data.get("total_time")[indic]))
            ax2.lines[i].set_ydata(np.append(ax2.lines[i].get_ydata(), data.get("total_impulses")[indic][idx]))

    def update_bar(self, num, data, ax_mgr2):
        first_idx, _ = data["total_trunk_indices"]
        indic = first_idx + num
        special_indices = data.get("special_indices")
        impulse_key = "total_impulses" if not self.OPTIONS.cumul else "cumul_impulses"
        for i, idx in enumerate(special_indices):
            new_height = data.get(impulse_key)[indic][idx]
            self.bars[i].set_height(new_height)
        return self.bars

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
        with h5py.File(self.OPTIONS.filename, 'r') as data_file:
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
            plt.savefig('{}.png'.format(os.path.splitext(self.OPTIONS.filename)[0]), bbox_inches='tight')
            # plt.savefig('{}.eps'.format(self.OPTIONS.filename), format='eps')#, dpi=1000)
        else:
            plt.show()


    #############################
    # Parallel video creation : #
    #############################

    def make_partial_floe_video(self, out_filename, data_chunk, version):
        "creates video directly from data"
        fig, ax = plt.subplots()
        if getattr(self.OPTIONS, "hd", False):
            fig.set_size_inches(20, 15)
            fig.set_dpi(100)
        init, update = self.get_anim_fcts(version)
        ax_mgr = AxeManager(ax)
        init(data_chunk, ax_mgr)
        if not data_chunk.get("static_axes"):
            ax.axis('equal') # automatic scale
        else:
            ax.axis(data_chunk.get("static_axes"))
        anim = animation.FuncAnimation(
            fig, update, len(data_chunk.get("time")), fargs=(data_chunk, ax_mgr),
            interval=1, blit=False)
        print(out_filename)
        anim.save(out_filename, writer=self.writer)

    def make_partial_floe_video_helper(self, t):
        return self.make_partial_floe_video(*t)


    def make_partial_floe_video_dual_plot(self, out_filename, data_chunk, version):
        "creates video directly from data, and adds obstacle impulse subplot under floes"
        try:
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            fig, ax = plt.subplots()
            if getattr(self.OPTIONS, "hd", False):
                fig.set_size_inches(16, 12)
                fig.set_dpi(100)
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1])
            ax_mgr1 = AxeManager(ax1)
            ax_mgr2 = AxeManager(ax2)
            init, update = self._init2_dual, self._update2_dual
            init(data_chunk, ax_mgr1, ax_mgr2)
            if not data_chunk.get("static_axes"):
                ax1.axis('equal') # automatic scale
            else:
                ax1.axis(data_chunk.get("static_axes"))
            # ax1.set_axis_bgcolor('#162252')
            anim = animation.FuncAnimation(
                fig, update, len(data_chunk.get("time")), fargs=(data_chunk, ax_mgr1, ax_mgr2),
                interval=50, blit=False)
            anim.save(out_filename, writer=self.writer)
        except Exception as e:
            print(f"Erreur dans le processus enfant : {e}")

    def make_partial_floe_video_dual_plot_helper(self, t):
        return self.make_partial_floe_video_dual_plot(*t)

    def plot_floes_fast(self, **kwargs):
        "Creates a video from hdf5 file output, using multiprocessing to go faster"
        with h5py.File(self.OPTIONS.filename, 'r') as data_file:
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
        partial_video_maker = self.make_partial_floe_video_helper if not self.OPTIONS.dual else self.make_partial_floe_video_dual_plot_helper
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
            file_time_dependant_keys =["time", "floe_states", "mass_center"]
        else:
            file_time_dependant_keys =["time", "floe_states"] # allow input files to be plotted
        if single_step == "OFF":
            for key in file_time_dependant_keys:
                d[key] = data_file.get(key)[:10000:self.OPTIONS.step]
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
                d["floe_shapes"] = [np.array(data_file.get("floe_shapes").get(k)) for k  in sorted(list(data_file.get("floe_shapes")), key=int)]
            else:
                d["floe_shapes"] = self.calc_shapes(data_file)
        # Other datas
        d["cumul_impulses"] = [np.array([state[6] for state in time_states]) for time_states in d["floe_states"]]
        d["impulses"] = d["total_impulses"] = self.calc_impulses(d["cumul_impulses"], 6) # Calc impulsions
        # Calculate speed norm for each floe at each time step
        d["speeds"] = [np.linalg.norm([np.array(state)[:,3], np.array(state)[:,4]], axis=0) for state in d.get("floe_states")]
        # calc or set global max impulse for color range
        d["MAX_IMPULSE"] = max(np.amax(step_impulses) for step_impulses in d.get("impulses"))
        d["MAX_SPEED"] = max(np.amax(step_speeds) for step_speeds in d.get("speeds"))
        d["color_key"] = "impulses" if not getattr(self.OPTIONS, "speed_color", True) else "speeds"
        d["max_color_key"] = "MAX_IMPULSE" if not getattr(self.OPTIONS, "speed_color", True) else "MAX_SPEED"
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
        if self.OPTIONS.dual:
            d["special_indices"] = getattr(self.OPTIONS, "index") or []
            d["max_received_impulse_special_indices"] = max(np.amax(
                [step_impulses[i] for i in d["special_indices"]])
                for step_impulses in d.get("cumul_impulses"))
        return d


    def _get_useful_trunk_data(self, data_global, trunk):
        """Slice global datas for a partial video"""
        d = {}
        trunked_keys = ["time", "floe_states", "impulses", "mass_center"]
        global_keys = [
            "total_time","total_impulses", "MAX_IMPULSE", "max_received_impulse_special_indices", "window",
            "special_indices", "static_axes", "ghosts_trans", "cumul_impulses"
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


    def calc_impulses(self, cumul_impulses, n=1):
        "calcul floe received impulses between each step and step -n"
        imps = [np.subtract(cumul_impulses[t], cumul_impulses[max(t-n, 0)]) for t in range(len(cumul_impulses))]
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
