#! /usr/bin/env python
# encoding: utf-8

print('→ loading myboost')

from waflib.Configure import conf

def options(opt):
    opt.add_option('--boostdir', action='store', default='', dest='boost')

@conf
def read_boost(ctx):
    ctx.start_msg('Checking for the variable BOOST')
    if ctx.options.boost:
        ctx.env.BOOST = ctx.options.boost
        ctx.end_msg(ctx.env.BOOST)
    else:
        ctx.env.BOOST = '/usr'
        #ctx.end_msg('BOOST is not set')

def configure(ctx): 
        ctx.read_boost()
