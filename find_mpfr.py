#! /usr/bin/env python
# encoding: utf-8

import os
print('→ loading mpfr')

from waflib.Configure import conf

def options(opt):
    opt.add_option('--mpfrdir', action='store', default='', dest='mpfr')

@conf
def read_mpfr(ctx):
    ctx.start_msg('Checking for the variable MPFR')
    if ctx.options.mpfr:
        ctx.env.MPFR = ctx.options.mpfr
        ctx.end_msg(ctx.env.MPFR)
    else:
        ctx.env.MPFR = '/usr'
        #ctx.end_msg('MPFR is not set')

def configure(conf): 
    conf.read_mpfr()
    mpfr_libs = ['mpfr']
    conf.env.LIBPATH_MPFR   = [os.path.join(conf.env.MPFR, 'lib')]
    conf.env.INCLUDES_MPFR   = [os.path.join(conf.env.MPFR, 'include')]
    for libname in mpfr_libs:
        conf.check_cxx(lib = libname, use = 'MPFR')

