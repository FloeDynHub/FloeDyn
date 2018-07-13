#! /usr/bin/env python
# encoding: utf-8

print('→ loading hdf5')

from waflib.Configure import conf

def options(opt):
    opt.add_option('--hdf5dir', action='store', default='', dest='hdf5')

@conf
def read_hdf5(ctx):
    ctx.start_msg('Checking for the variable HDF5')
    if ctx.options.hdf5:
        ctx.env.HDF5 = ctx.options.hdf5
        ctx.end_msg(ctx.env.HDF5)
    else:
        ctx.env.HDF5 = '/usr'
        #ctx.end_msg('HDF5 is not set')

def configure(conf): 
    conf.read_hdf5()
    hdf5_libs = ['hdf5', 'hdf5_cpp']
    conf.env.LIBPATH_HDF5   = [conf.env.HDF5 + "lib"]
    conf.env.INCLUDES_HDF5   = [conf.env.HDF5 + "include"]
    for libname in hdf5_libs:
        conf.check_cxx(lib = libname, use = 'HDF5')

