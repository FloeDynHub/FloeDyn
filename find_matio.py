#! /usr/bin/env python
# encoding: utf-8

print('→ loading matio')

from waflib.Configure import conf

def options(opt):
    opt.add_option('--matiodir', action='store', default='', dest='matio')

@conf
def read_matio(ctx):
    ctx.start_msg('Checking for the variable MATIO')
    if ctx.options.matio:
        ctx.env.MATIO = ctx.options.matio
        ctx.end_msg(ctx.env.MATIO)
    else:
        ctx.end_msg('MATIO is not set')

def configure(ctx): 
        ctx.read_matio()
