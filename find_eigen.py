#! /usr/bin/env python
# encoding: utf-8

print('→ loading eigen')

from waflib.Configure import conf

def options(opt):
    opt.add_option('--eigendir', action='store', default='', dest='eigen')

@conf
def read_eigen(ctx):
    ctx.start_msg('Checking for the variable EIGEN')
    if ctx.options.eigen:
        ctx.env.EIGEN = ctx.options.eigen
        ctx.end_msg(ctx.env.EIGEN)
    else:
        ctx.env.EIGEN = '/usr'
        #ctx.end_msg('EIGEN is not set')

def configure(ctx): 
        ctx.read_eigen()
