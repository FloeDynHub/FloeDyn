#! /usr/bin/env python
# encoding: utf-8

import os
print('→ loading gmp')

from waflib.Configure import conf

def options(opt):
    opt.add_option('--gmpdir', action='store', default='', dest='gmp')

@conf
def read_gmp(ctx):
    ctx.start_msg('Checking for the variable GMP')
    if ctx.options.gmp:
        ctx.env.GMP = ctx.options.gmp
        ctx.end_msg(ctx.env.GMP)
    else:
        ctx.env.GMP = '/usr'
        #ctx.end_msg('GMP is not set')

def configure(conf): 
    conf.read_gmp()
    gmp_libs = ['gmp']
    conf.env.LIBPATH_GMP   = [os.path.join(conf.env.GMP, 'lib')]
    conf.env.INCLUDES_GMP   = [os.path.join(conf.env.GMP, 'include')]
    for libname in gmp_libs:
        conf.check_cxx(lib = libname, use = 'GMP')

