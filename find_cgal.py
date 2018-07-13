#! /usr/bin/env python
# encoding: utf-8

import os
print('→ loading cgal')

from waflib.Configure import conf

def options(opt):
    opt.add_option('--cgaldir', action='store', default='', dest='cgal')

@conf
def read_cgal(ctx):
    ctx.start_msg('Checking for the variable CGAL')
    if ctx.options.cgal:
        ctx.env.CGAL = ctx.options.cgal
        ctx.end_msg(ctx.env.CGAL)
    else:
        ctx.env.CGAL = '/usr'
        #ctx.end_msg('CGAL is not set')

def configure(conf): 
    conf.read_cgal()
    cgal_libs = ['CGAL']
    conf.env.LIBPATH_CGAL   = [os.path.join(conf.env.CGAL, 'lib')]
    conf.env.INCLUDES_CGAL   = [os.path.join(conf.env.CGAL, 'include')]
    for libname in cgal_libs:
        conf.check_cxx(lib = libname, use = 'CGAL')

