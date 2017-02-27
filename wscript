#! /usr/bin/env python
# encoding: utf-8

# Written by Quentin Jouet

"""
Documentation : https://waf.io/book/

USAGE
Build main product : ./waf --target <targ>
    with <targ> = FLOE or FLOE_PBC (PBC for Periodic Boundary Conditions)

Build and run unit tests : ./waf TESTS
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
    # Build opts :
    # opt.add_option('--run', action='store_true', default=False, dest='run')
    opt.add_option('--name', action='store', default="", dest='name')
    opt.add_option('-t','--target', action='store', default="", dest='target')
    opt.add_option(
        '--debug', action='store_true', default=False, dest='debug')
    opt.add_option('--omp', action='store_true', default=False, dest='omp')
    # configure opts :
    opt.add_option('--gcc', action='store_true', default=False, dest='gcc')
    opt.add_option('--icc', action='store_true', default=False, dest='icc')
    # opt.add_option('--mpi', action='store_true', default=False, dest='mpi')


def configure(conf):
    conf.check_waf_version(mini='1.8.8')
    if conf.options.gcc:
        conf.load('g++')
        # if conf.options.mpi:
        #     conf.env['CXX'] = 'mpicxx'
    elif conf.options.icc:
        conf.load('icpc')  
    else:
        conf.load('compiler_cxx')
    conf.check_cfg(atleast_pkgconfig_version='0.0.0')
    conf.check_boost(lib='system filesystem', mandatory=False)
    conf.check_cfg(
        package='matio', args=['matio >= 1.5.2', '--cflags', '--libs'],
        msg="Checking for 'matio 1.5.2'", mandatory=False)


def recursive_file_finder(folder="", pattern=""):
    file_list = []
    for root, dirnames, filenames in os.walk(folder):
        for filename in fnmatch.filter(filenames, pattern):
            file_list.append(os.path.join(root, filename))
    return file_list


def list_test_cpp(pattern=""):
    """return all .cpp filename in tests/ with name begining by 'test'"""
    cpp_list = ["tests/main.cpp"] + recursive_file_finder('tests', 'test*{}*.cpp'.format(pattern) )
    return cpp_list


def find_STEST_cpp(pattern=""):
    """
    return a .cpp filname if only 1 is found in tests/containing
    the given pattern
    """
    cpp_list = recursive_file_finder('tests', 'STEST*{}*.cpp'.format(pattern) )
    if len(cpp_list) > 1:
        print(len(cpp_list), " corresponding tests found :")
        print(" ".join(cpp_list))
    elif len(cpp_list) == 0:
        print("No corresponding test found")
    else:
        return cpp_list


def clean_tests(ctx):
    print('cleaning tests...')
    ctx.exec_command('./waf clean --target TESTS')


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
                      # '/usr/local/include/siconos',
                    ], #+ [path for path in os.environ["PATH"].split(":") if not "bin" in path],
        "lib": ['boost_system',
                'boost_program_options',
                'matio',
                "hdf5",
                "hdf5_cpp",
                "CGAL", "gmp", "mpfr", "boost_thread",
                # "siconos_numerics"
                ],
        "libpath": ["/usr/local/lib", "/usr/lib"], # + os.environ.get("LD_LIBRARY_PATH", "/").split(":"),
        "framework": ["Accelerate"],
        "frameworkpath" : ["/System/Library/Frameworks"]
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
                 "-Wall", "-Wextra", #"-Wshadow",
                 "-Wno-unused-parameter",
                 "-Wno-gnu-anonymous-struct", "-Wno-nested-anon-types", # floe/geometry/geometries/point.hpp
                 "-pedantic"
             ],
            "defines": ["NDEBUG"]
        })
    OPTION_DICT["cxxflags"].extend(os.environ.get("CFLAGS", "").split(" "))
    OPTION_DICT["linkflags"].extend(os.environ.get("LDFLAGS", "").split(" "))
    return OPTION_DICT

import subprocess

def build(bld):
    opts = get_option_dict(bld.options.debug)
    if bld.options.omp:
        opts["linkflags"].append("-fopenmp")
        opts["cxxflags"].append("-fopenmp")
    if "MPI" in bld.options.target:
        opts["linkflags"].extend(["-lmpi_cxx", "-lmpi"])
        opts["defines"].append('MPIRUN')
        opts["cxxflags"].extend(subprocess.check_output(["mpicc", "--showme:compile"]).strip().split(" "))
        opts["linkflags"].extend(subprocess.check_output(["mpicc", "--showme:link"]).strip().split(" "))
        # print opts["linkflags"]
    if bld.options.target in ["FLOE", "FLOE_PBC", "FLOE_MPI"]:
        opts["source"] = ["product/FLOE.cpp"] + recursive_file_finder("src/floe", "*.cpp")
        opts["target"] = bld.options.target
        if bld.options.target == "FLOE_PBC":
            opts["defines"].append('PBC')
        bld.program(**opts)
    elif "FLOE" in bld.options.target:
        print "TARGET", bld.options.target
        opts["source"] = ["product/{}.cpp".format(bld.options.target)] + recursive_file_finder("src/floe", "*.cpp")
        opts["target"] = bld.options.target
        bld.program(**opts)
    elif bld.options.target == "TESTS":  # using Catch framework
        opts["source"] = list_test_cpp(bld.options.name) + recursive_file_finder("src/floe", "*.cpp")
        opts["target"] = test_target
        bld.program(**opts)
    elif bld.options.target == "TEST":
        opts["source"] = find_STEST_cpp(bld.options.name) or bld.fatal('STOP')
        opts["target"] = TEST_target
        bld.program(**opts)
    else:
        print("Nothing to build.")


def forward_options(opt_list, options):
    def my_str(var):
        return " " + var if isinstance(var, basestring) else ""
    def format_opt(opt):
        opt_val = getattr(options, opt, "")
        return "--{}{}".format(opt, my_str(opt_val)) if opt_val else ""
    return " ".join([format_opt(opt) for opt in opt_list ])


def TESTS(ctx):
    """Build, run and clean tests using Catch Framework"""
    err = ctx.exec_command('./waf build --target TESTS {}'.format(
        forward_options(["omp", "name", "debug"], ctx.options)))
    if not err:
        run_tests(ctx)
        # clean_tests(ctx)


def TEST(ctx):
    """Build a single test ("^STEST_" files) (to run manually)"""
    ctx.exec_command('./waf build --target TEST {}'.format(
        forward_options(["omp", "name", "debug"], ctx.options)))
    print("to run the test : ./build/%s <args>" % TEST_target)

