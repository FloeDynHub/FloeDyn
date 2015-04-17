#! /usr/bin/env python
# encoding: utf-8

# Written by Quentin Jouet

"""
Documentation : https://waf.io/book/

require extra waf tool boost.py (TODO : auto detect ?)

USAGE
Build and run unit tests : waf test
(--name <str> : restrict to test filenames containing str)

Build a single test ("^STEST_" files) : waf TEST --name <filename>
"""

import os
import fnmatch
from subprocess import call

top = "."
out = "build"

test_target = 'catchtest'
TEST_target = "STEST"


def options(opt):
    opt.load('compiler_cxx boost')
    # opt.add_option('--run', action='store_true', default=False, dest='run')
    opt.add_option('--name', action='store', default="", dest='name')
    opt.add_option('--target', action='store', default="", dest='target')


def configure(conf):
    conf.check_waf_version(mini='1.8.8')
    conf.load('compiler_cxx boost')
    conf.check_cfg(atleast_pkgconfig_version='0.0.0')
    conf.check_boost(lib='system filesystem')
    conf.check_cfg(
        package='matio', args=['matio >= 1.5.2', '--cflags', '--libs'],
        msg="Checking for 'matio 1.5.2'", mandatory=False)


def list_test_cpp(pattern=""):
    """return all .cpp filename in tests/ with name begining by 'test'"""
    cpp_list = ["tests/main.cpp"]
    for root, dirnames, filenames in os.walk('tests'):
        for filename in fnmatch.filter(filenames, 'test*.cpp'):
            if pattern in filename:
                cpp_list.append(os.path.join(root, filename))
    return cpp_list


def find_STEST_cpp(pattern=""):
    """
    return a .cpp filname if only 1 is found in tests/containing
    the given pattern
    """
    cpp_list = []
    for root, dirnames, filenames in os.walk('tests'):
        for filename in fnmatch.filter(filenames, 'STEST*.cpp'):
            if pattern in filename:
                cpp_list.append(os.path.join(root, filename))
    if len(cpp_list) > 1:
        print(len(cpp_list), " corresponding tests found :")
        print(" ".join(cpp_list))
    elif len(cpp_list) == 0:
        print("No corresponding test found")
    else:
        return cpp_list


def clean_tests(ctx):
    print('cleaning tests...')
    ctx.exec_command('waf clean --target unittests')


def run_tests(ctx):
    print('running tests...')
    # ctx.exec_command("%s/%s" % (out, test_target))
    call("%s/%s" % (out, test_target))  # for terminal coloration !


OPTION_DICT = {
    "includes": ['./src', '/usr/local/include'],
    "lib": ['boost_timer', 'boost_chrono', 'boost_system', 'matio'],
    "libpath": ['/usr/lib', '/usr/local/lib'],
    "linkflags": ['-g'],
    "cxxflags": ['-std=c++11', "-Wall"],
}


def build(bld):
    opts = OPTION_DICT
    if bld.options.target == "unittests":  # using Catch framework
        opts["source"] = list_test_cpp(bld.options.name)
        opts["target"] = test_target
        bld.program(**opts)
        # running tests
        bld.add_post_fun(run_tests)
    elif bld.options.target == "TEST":
        opts["source"] = find_STEST_cpp(bld.options.name) or bld.fatal('STOP')
        opts["target"] = TEST_target
        bld.program(**opts)
    else:
        print("nothing to build yet")


def test(ctx):
    """Build, run and clean tests using Catch Framework"""
    name_opt = "--name %s" % ctx.options.name if ctx.options.name else ""
    ctx.exec_command('waf build --target unittests %s' % name_opt)
    # ctx.exec_command('waf run_tests')
    ctx.exec_command('waf clean_tests')


def TEST(ctx):
    """Build a single test ("^STEST_" files) (to run manually)"""
    ctx.exec_command('waf build --target TEST --name %s' % ctx.options.name)
    print("to run the test : ./build/%s <args>" % TEST_target)
