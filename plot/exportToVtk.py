#!/usr/bin/env python

def export_floes(data_file, num, data):
    export_mesh(data_file, 0,data)
    return True 

def export_mesh(data_file, num, data):
    # Exports mesh. 
    begin, step = 0, 1
    indic = begin + step * num

    coord = data.get("floe_meshes_coord")
    connect = data.get("floe_meshes_connect")
    nb_floes = len(data.get("floe_shapes"))
    if (coord and connect):
        count = 0
        elements = []
        for iFloe in np.arange(len(coord)):
            for iElem in np.arange(len(connect[iFloe])):
                elements.append(np.array([coord[iFloe][int(connect[iFloe][iElem][0])], coord[iFloe][int(connect[iFloe][iElem][1])], coord[iFloe][int(connect[iFloe][iElem][2])]]))
        
        coord = transform_meshes(data_file, data.get("floe_meshes_coord"), data.get("floe_states")[indic], data.get("mass_center")[indic], False)  
        # connect = [v for i, v in enumerate(connect) if data.get("floe_states")[indic][i][9] == 1]
        # coord = [v for i, v in enumerate(coord) if data.get("floe_states")[indic][i][9] == 1]            


        color_data = data.get("floe_elem_data")
        color_data_indic = indic 
        cdTemp = np.zeros(np.array(color_data).shape)
        k = 0
        for iFloe in range(0, nb_floes):
            # print(data.get("floe_states").shape)
            if data.get("floe_states")[indic,iFloe,9] == 1:
            # if data.get("floe_states")[indic][iFloe][9] == 1:
                cdTemp[:,k,:] = color_data[:,iFloe,:]
                k = k+1
        color_data = cdTemp
        return True 
    else:
        print('could not load coord and connect')
        return False

def get_useful_data(data_file, single_step="OFF", img=False):
    d = {}
    # File datas
    file_time_dependant_keys =["time", "floe_states", "mass_center", "floe_elem_data"]
    
    for key in file_time_dependant_keys:
        d[key] = data_file.get(key)[1]
    d["total_time"] = d["time"]
    if data_file.get("floe_shapes") is not None:
        d["floe_shapes"] = [np.array(data_file.get("floe_shapes").get(k)) for k in sorted(list(data_file.get("floe_shapes")), key=int)]
        if data_file.get("floe_meshes_coord") is not None:
            d["floe_meshes_coord"] = [np.array(data_file.get("floe_meshes_coord").get(k)) for k in sorted(list(data_file.get("floe_meshes_coord")), key=int)]
        if data_file.get("floe_meshes_connect") is not None:
            d["floe_meshes_connect"] = [np.array(data_file.get("floe_meshes_connect").get(k)) for k in sorted(list(data_file.get("floe_meshes_connect")), key=int)]
        
    # d["impulses"] = d["total_impulses"] = self.calc_impulses(d["floe_states"], 6) # Calc impulsions
    # d["MAX_IMPULSE"] = max(np.amax(step_impulses) for step_impulses in d.get("impulses"))
    w = data_file.get("window")
    w_width, w_length = w[1] - w[0], w[3] - w[2]
    print('hello')
    print(data_file.get("floe_states").shape)
    print(len(d["floe_shapes"]))
    return d

def transform_meshes(data_file, floe_coords, floe_states, floe_mass_centers, follow=False):
    pos_ids = (0,1) if follow else (7,8) # (7,8) contains translated position in fixed initial window
    coord = floe_coords
    def rotation_mat(theta):
        return np.array([[np.cos(theta), -np.sin(theta)], 
                         [np.sin(theta),  np.cos(theta)]])
    # for x in data_file.get("floe_states"):
    #     print('et voilà x : {}'.format(x))
    rots = [rotation_mat(x[2]) for x in data_file.get("floe_states")]
    coord = [np.transpose(np.dot(rot, np.transpose(mesh))) for rot, mesh in zip(rots, floe_coords)]
    coord = [np.add(mesh, np.repeat([[x[pos_ids[0]], x[pos_ids[1]]]], len(mesh), axis=0)) for x, mesh in zip(floe_states, coord)]
    return coord


def read_data(data_file):
    d = {}
    file_time_dependant_keys =["time", "floe_states", "mass_center", "floe_elem_data", "floe_node_data"]
    for key in file_time_dependant_keys:
        d[key] = data_file.get(key)

    d["floe_shapes"] = [np.array(data_file.get("floe_shapes").get(k)) for k in sorted(list(data_file.get("floe_shapes")), key=int)]
    d["floe_meshes_coord"] = [np.array(data_file.get("floe_meshes_coord").get(k)) for k in sorted(list(data_file.get("floe_meshes_coord")), key=int)]
    d["floe_meshes_connect"] = [np.array(data_file.get("floe_meshes_connect").get(k)) for k in sorted(list(data_file.get("floe_meshes_connect")), key=int)]
    
    return d 
    
        
    # d["impulses"] = d["total_impulses"] = self.calc_impulses(d["floe_states"], 6) # Calc impulsions
    # d["MAX_IMPULSE"] = max(np.amax(step_impulses) for step_impulses in d.get("impulses"))
    # w = data_file.get("window")
    # w_width, w_length = w[1] - w[0], w[3] - w[2]
    # print('hello')
    # print(data_file.get("floe_states").shape)
    # print(len(d["floe_shapes"]))
    # return d

def export_vtk(filename, coord, connect, data, dataname):
    nNodes = coord.shape[0]
    nElem = connect.shape[0]
    if connect.shape[1] != 3:
        print('I did not expect the second dimension og connect to be {}'.format(connect.shape[1]))
    
    # print('trying to export a mesh containing {} nodes and {} elements'.format(nNodes, nElem))
    # coord = np.array([[0,0],[0,1],[1,0],[1,1]])
    # connect = np.array([[0,1,2],[1,3,2]])
    # elemType = np.array([5,5])
    elemType = 5*np.ones([nElem, 1]) # vtk element type 5 corresponds to linear triangles. To be improved 
    nTot = 4*nElem # again, I expect only linear triangles 

    f = open(filename, "w")
    f.write("# vtk DataFile Version 2.0\n")
    f.close()
    f = open(filename, "a") 
    f.write("VTK from FloeDyn\n")
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    f.write("POINTS {} float\n".format(nNodes))
    np.savetxt(f, np.concatenate((coord, np.zeros((nNodes, 1))), axis=1), fmt='%4.6f', delimiter=' ',)
    f.write("\nCELLS {} {}\n".format(nElem, nTot))
    np.savetxt(f, np.concatenate((3*np.ones((nElem, 1)),connect), axis=1), fmt='%d', delimiter=' ',)
    f.write("\nCELL_TYPES {}\n".format(nElem))
    np.savetxt(f, elemType, fmt='%d', delimiter=' ',)
    f.close()
    return True 


def add_vtk_elem_data(filename, coord, connect, data, dataname):
    # nNodes = data.shape[0]
    nElem = data.shape[0]
    # if connect.shape[1] != 3:
    #     print('I did not expect the second dimension of connect to be {}'.format(connect.shape[1]))
    # if data.shape[0] != nElem:
    #     print('Dimension mismatch, I got {} instead of {}'.format(data.shape[0], nElem))
    f = open(filename, "a")
    f.write("\nCELL_DATA {}\n".format(nElem))
    f.write("SCALARS {} float 1\nLOOKUP_TABLE default\n".format(dataname)) 
    np.savetxt(f, data, fmt='%f', delimiter=' ',)
    # # f.write("SCALARS {} float 1\nLOOKUP_TABLE default\n".format(dataname2))
    # # np.savetxt(f, data, fmt='%d', delimiter=' ',)
    f.close()
    return True 

def add_vtk_node_data(filename, coord, connect, data, dataname):
    # nNodes = data.shape[0]
    nNodes = coord.shape[0]
    # print('est ce que la taille ce serait pas {} par hasard ?'.format(coord.shape[0]))
    # nElem = connect.shape[0]
    # if connect.shape[1] != 3:
    #     print('I did not expect the second dimension of connect to be {}'.format(connect.shape[1]))
    if data.shape[0] == nNodes:
        f = open(filename, "a")
        f.write("\nPOINT_DATA {}\n".format(nNodes))
        f.write("SCALARS {} float 1\nLOOKUP_TABLE default\n".format(dataname))
        np.savetxt(f, data, fmt='%f', delimiter=' ',)
        # # f.write("SCALARS {} float 1\nLOOKUP_TABLE default\n".format(dataname2))
        # # np.savetxt(f, data, fmt='%d', delimiter=' ',)
        # print("Writing scalar data, shape is {}".format(data.shape))
        f.close()
    elif data.shape[0] == nNodes*2:
        f = open(filename, "a")
        f.write("\nPOINT_DATA {}\n".format(nNodes))
        # f.write("SCALARS {} float 2\nLOOKUP_TABLE default\n".format(dataname)) # VECTORS would take three components, but we're in 2D 
        f.write("VECTORS {} float\n".format(dataname)) # VECTORS take three components, but we're in 2D 
        # np.savetxt(f, data.reshape([nNodes, 2]), fmt='%f %f', delimiter=' ',)
        np.savetxt(f, np.concatenate((data.reshape([nNodes, 2]), np.zeros([nNodes,1])), axis=1), fmt='%f %f %f', delimiter=' ',)
        f.close()
    else:
        print('Dimension mismatch, I got {} instead of {}'.format(data.shape[0], nNodes))
    return True 

def data_at_t(data, t, displacement=False):
    
    def rotation_mat(theta):
        return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta),  np.cos(theta)]])
    
    pos_ids = (7,8) # (7,8) contains translated position in fixed initial window

    nFloes = data.get("floe_elem_data")[0].shape[0]
    connect = np.zeros((0,3))
    coord = np.zeros((0,2))
    elem_field = np.zeros((0,1))
    node_field = np.zeros((0,1))
    shift = 0
    for iFloe in range(0, nFloes):
        if data.get("floe_states")[t,iFloe,9] == 1:
            nNodes = data.get("floe_meshes_coord")[iFloe].shape[0]
            nElem = data.get("floe_meshes_connect")[iFloe].shape[0]
            r = rotation_mat(data.get("floe_states")[t, iFloe, 2])
            # concaténation des coordonnées et rotation / translation 
            coordShifted = data.get("floe_meshes_coord")[iFloe]
            coordShifted = np.transpose(np.dot(r, np.transpose(coordShifted)))
            coordShifted = np.add(coordShifted, [[data.get("floe_states")[t, iFloe, pos_ids[0]], data.get("floe_states")[t, iFloe, pos_ids[1]]]])
            coord = np.concatenate((coord, coordShifted), axis=0)
            # concatenation des connectivités 
            connect = np.concatenate((connect, data.get("floe_meshes_connect")[iFloe]+shift))
            shift += nNodes
            dataTemp = np.zeros((nElem, 1))
            # print(f"nElem = {nElem}")
            temp = data.get("floe_elem_data")[t, iFloe, :nElem]
            if len(temp) != nElem:
                dataTemp[:nElem,0] = np.zeros([nElem,])
            else:
                dataTemp[:nElem,0] = temp
            elem_field = np.concatenate((elem_field, dataTemp))
            dataTempNode = np.zeros((nNodes*2, 1))
            temp = data.get("floe_node_data")[t, iFloe, :nNodes*2]
            if len(temp) != nNodes:
                dataTempNode[:nNodes*2,0] = np.zeros([nNodes*2,])
            else:
                dataTempNode[:nNodes*2,0] = data.get("floe_node_data")[t, iFloe, :nNodes*2]
            # the solution is computed in the reference frame of the floe, so it has to be rotated as well 
            dataTempNode = dataTempNode.reshape([nNodes, 2])
            dataTempNode = np.transpose(np.dot(r, np.transpose(dataTempNode)))
            dataTempNode = dataTempNode.reshape([2*nNodes, 1])
            
            node_field = np.concatenate((node_field, dataTempNode))
            
                
    d = {}
    d["coord"] = coord
    d["connect"] = connect
    d["elem_field"] = elem_field
    d["node_field"] = node_field
    return d

import numpy as np
from subprocess import call
from utils import filename_without_extension, get_unused_path, check_path_existence, mkdir_path
import h5py
import datetime
import math
import sys 
if len(sys.argv) == 1:
    filename = "../io/outputs/0_test.h5"
    print(f"you did not specify output file name. Default file is {filename}")
else:
    filename = (sys.argv[1])
    print("Converting {}.".format(filename))
data_file = h5py.File(filename, 'r')
rawData = read_data(data_file)

nNodes = rawData.get("floe_meshes_coord")[0].shape[0]
nElem = rawData.get("floe_meshes_connect")[0].shape[0]
nTime = len(rawData.get("time"))
time = np.array(rawData.get("time"))
data = rawData 
# print(data.get("floe_states").shape)
# print(data.get("floe_meshes_coord")[0].shape)
# print(data.get("floe_meshes_connect")[0].shape)
# print(data.get("floe_elem_data")[0].shape)
# export_vtk(filename, data.get("floe_meshes_coord")[0], data.get("floe_meshes_connect")[0], np.zeros((nNodes, 1)), dataname)
print('the mesh contains {} elements and {} nodes.'.format(nElem, nNodes))
print('size de node data : {}'.format(data["floe_node_data"].shape))
print('size de elem data : {}'.format(data["floe_elem_data"].shape))

# loop over the time steps, faut translate/rotate puis extract only useful floes at give time step 
for iTime in range(0, nTime):
    data = data_at_t(rawData, iTime, iTime > 0)
    timing = time[iTime]
    # print('exporting time {} '.format(timing))
    filename = '../io/outputs/vtk/test_{}.vtk'.format(iTime)
    dataname = 'trucTruc'
    nNodes = data["coord"].shape[0]
    nElem = data["connect"].shape[0]
    export_vtk(filename, data["coord"], data["connect"], np.zeros((nNodes, 1)), dataname)
    # add_vtk_node_data(filename, data["coord"], data["connect"], np.arange((nNodes, 1)), 'trucTrucNoeuds')
    # add_vtk_node_data(filename, data["coord"], data["connect"], data["node_field"], 'displacement_x')
    add_vtk_node_data(filename, data["coord"], data["connect"], data["node_field"], 'displacement_x')
    # print(data["node_field"])
    add_vtk_elem_data(filename, data["coord"], data["connect"], data["elem_field"], 'sigma_von_mises')
    # add_vtk_elem_data(filename, data["coord"], data["connect"], np.zeros((nElem, 1)), 'trucTrucElem')

# print(data.get("time") )

# int(data.get("time")[indic])
# ax.set_title("t = {}".format(str(datetime.timedelta(seconds=int(data.get("time")[indic])))))
