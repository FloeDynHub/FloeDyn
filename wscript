#! /usr/bin/env python
# encoding: utf-8

# Written by Quentin Jouet

"""
Documentation : https://waf.io/book/

USAGE
Build main product : ./waf --target <targ>
    whith <targ> = FLOE or FLOE_PBC (PBC for Periodic Boundary Conditions)

Build and run unit tests : ./waf test
    options:
        --name <str> : restrict to test filenames containing str)

Build a single test ("^STEST_" files) : ./waf TEST --name <filename>

Compilation options (available for tests too) :
    --debug to get debug product
    --omp to get omp parallelisation (when compiler allows it)

"""

import os
import fnmatch
from subprocess import call, Popen, PIPE, STDOUT
import time

top = "."
out = "build"

test_target = 'catchtest'
TEST_target = "STEST"

def timeit(func):
    """decorator to measure function execution time"""
    def timed(*args, **kw):
        ts = time.time()
        result = func(*args, **kw)
        te = time.time()
        print '%r %2.2f s' % (func.__name__, te-ts)
        return result

    return timed


def options(opt):
    opt.load('compiler_cxx gxx boost')
    # opt.add_option('--run', action='store_true', default=False, dest='run')
    opt.add_option('--name', action='store', default="", dest='name')
    opt.add_option('--target', action='store', default="", dest='target')
    opt.add_option(
        '--debug', action='store_true', default=False, dest='debug')
    opt.add_option('--gcc', action='store_true', default=False, dest='gcc')
    opt.add_option('--icc', action='store_true', default=False, dest='icc')
    opt.add_option('--omp', action='store_true', default=False, dest='omp')


def configure(conf):
    conf.check_waf_version(mini='1.8.8')
    if conf.options.gcc:
        conf.load('g++')
    elif conf.options.icc:
        conf.load('icpc')
    else:
        conf.load('compiler_cxx')
    conf.check_cfg(atleast_pkgconfig_version='0.0.0')
    conf.check_boost(lib='system filesystem', mandatory=False)
    conf.check_cfg(
        package='matio', args=['matio >= 1.5.2', '--cflags', '--libs'],
        msg="Checking for 'matio 1.5.2'", mandatory=False)


def list_test_cpp(pattern=""):
    """return all .cpp filename in tests/ with name begining by 'test'"""
    # use glob instead ?
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
    ctx.exec_command('./waf clean --target unittests')


@timeit
def run_tests(ctx):
    print('running tests...')
    call("%s/%s" % (out, test_target))


def get_option_dict(debug=True):
    OPTION_DICT = {
        "includes": [ '../src',
                      '/usr/local/include',
                      '/usr/include',
                      '/usr/local/include/eigen3',
                    ] + [path for path in os.environ["PATH"].split(":") if not "bin" in path],
        "lib": ['boost_system',
                'matio',
                "hdf5",
                "hdf5_cpp",
                "CGAL", "gmp", "mpfr", "boost_thread"
                ],
        "libpath": ["/usr/local/lib", "/usr/lib"] + os.environ.get("LD_LIBRARY_PATH", "/").split(":")
    }
    if debug:
        OPTION_DICT.update({
            "linkflags": ['-g'],
            "cxxflags": [
                '-std=c++11',
                 '-O0',
                 "-Wall", #"-Wextra",
            ],
            "defines": []
        })
    else:
        OPTION_DICT.update({
            "linkflags": [],
            "cxxflags": [
                '-std=c++11',
                 "-O3",
                 # "-march=native", # g++ fails with this
                 "-mtune=native",
                 "-Wall",# "-Wextra",
             ],
            "defines": ["NDEBUG"]
        })
    return OPTION_DICT


def build(bld):
    opts = get_option_dict(bld.options.debug)
    if bld.options.omp:
        opts["linkflags"].append("-fopenmp")
        opts["cxxflags"].append("-fopenmp")
    if bld.options.target == "unittests":  # using Catch framework
        opts["source"] = list_test_cpp(bld.options.name)
        opts["target"] = test_target
        bld.program(**opts)
    elif bld.options.target == "TEST":
        opts["source"] = find_STEST_cpp(bld.options.name) or bld.fatal('STOP')
        opts["target"] = TEST_target
        bld.program(**opts)
    elif bld.options.target in ["FLOE", "FLOE_PBC"]:
        opts["source"] = ["product/FLOE.cpp"]
        opts["target"] = bld.options.target
        if bld.options.target == "FLOE_PBC":
            opts["defines"].append('PBC')
        bld.program(**opts)
    else:
        print("nothing to build yet")


def test(ctx):
    """Build, run and clean tests using Catch Framework"""
    # TODO find a simple way to forward options
    omp_opt = " --omp" if ctx.options.omp else ""
    name_opt = " --name %s" % ctx.options.name if ctx.options.name else ""
    debug_opt = " --debug" if ctx.options.debug else ""
    err = ctx.exec_command('./waf build --target unittests%s%s%s' % (
        name_opt, omp_opt, debug_opt))
    if not err:
        run_tests(ctx)
        # clean_tests(ctx)


def TEST(ctx):
    """Build a single test ("^STEST_" files) (to run manually)"""
    omp_opt = " --omp" if ctx.options.omp else ""
    name_opt = " --name %s" % ctx.options.name if ctx.options.name else ""
    debug_opt = " --debug" if ctx.options.debug else ""
    ctx.exec_command('./waf build --target TEST%s%s%s' % (
        name_opt, omp_opt, debug_opt))
    print("to run the test : ./build/%s <args>" % TEST_target)

