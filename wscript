#! /usr/bin/env python
# encoding: utf-8

# Written by Quentin Jouet

"""
Documentation : https://waf.io/book/

USAGE
configuration with: ./waf configure
eventually if one or more dependencies and/or libraries are not found,
one can use:        ./waf configure --depdir=path_to_dep_include
example with eigen: ./waf configure --eigendir=/usr/local/
or:                 ./waf configure --eigendir=/usr/local/Cellar/eigen/3.3.5/

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
from waflib.Configure import conf

out = "build"

test_target = 'catchtest'
TEST_target = "STEST"

def timeit(func):
    """decorator to measure function execution time"""
    def timed(*args, **kw):
        ts = time.time()
        result = func(*args, **kw)
        te = time.time()
        print( '%r %2.2f s' % (func.__name__, te-ts) )
        return result

    return timed


floedyn_deps = {
    'gmp' : ['gmp'],
    'boost' : ['boost_system', 'boost_program_options'],
    'eigen' : [], # header only
    'matio' : ['matio'],
    'hdf5'  : ['hdf5_cpp'],
    'cgal'  : [''],
    'mpfr'  : ['mpfr'],
    'cereal' : [''], # header only
    }

floedyn_includes = {
    'gmp' : [],
    'boost' : [],
    'eigen' : ['eigen3'],
    'matio' : [],
    'hdf5'  : ['hdf5/serial/'],
    'cgal'  : [],
    'mpfr'  : [],
    'cereal' : [],
    }


def options(opt):
    opt.load('compiler_cxx gxx')
    # Build opts :
    # opt.add_option('--run', action='store_true', default=False, dest='run')
    opt.add_option('--name', action='store', default="", dest='name')
    opt.add_option('-t','--target', action='store', default="", dest='target')
    opt.add_option('--debug', action='store_true', default=False, dest='debug', help='to active assertions')
    opt.add_option('-D','--defmacro', action='store', default="", dest='defmacro',
     help="MULTIOUTPUT: for using multiple output files. The size of a floe selection is required.\n"

        "LCPSTATS: for saving LCP statistics. The LCP max number saved is required (still below to 30000 to avoid memory overflow!)" )
    opt.add_option('--omp', action='store_true', default=False, dest='omp')
    # configure opts :
    opt.add_option('--gcc', action='store_true', default=False, dest='gcc')
    opt.add_option('--icc', action='store_true', default=False, dest='icc')
    opt.add_option('--default-search-path', action='store', default='/usr', dest='default_search_path')
    opt.add_option('--with-nix', action='store', default=False, dest='nix_build_inputs')
    opt.add_option('--enable-mpi', action='store_true', default=False, dest='mpi_on')
    # #for dep in floedyn_deps:
    # #    opt.load('find_' + dep, tooldir='.')
    # # opt.add_option('--mpi', action='store_true', default=False, dest='mpi')

    for dep in floedyn_deps:
        optname = '--' + dep + 'dir'
        opt.add_option(optname, action='store', default=None, dest=dep)

@conf
def configure_package(conf, name, required_libs=None, includes_suffix=None):
    print(f'---> Start conf for package {name} ...')
    searchpath = getattr(conf.env, name, conf.env.default_search_path)
    libpath_name = 'LIBPATH_' + name.upper()
    includespath_name = 'INCLUDES_' + name.upper()
    if isinstance(searchpath, list):
        lib_paths = []
        includes_list = []
        for p in searchpath:
            lib_paths.append(os.path.join(p, 'lib'))
            includes_list.append(os.path.join(p, 'include'))
            for suffix in includes_suffix:
                includes_list.append(os.path.join(os.path.join(p, 'include'), suffix))

    else:
        lib_paths = [os.path.join(searchpath, 'lib')]
        incpath = os.path.join(searchpath, 'include')
        includes_list = [incpath]
        for suffix in includes_suffix:
            includes_list.append(os.path.join(incpath, suffix))
    setattr(conf.env, libpath_name, lib_paths)
    setattr(conf.env, includespath_name, includes_list)
    if required_libs is None:
        required_libs = [name]

    print(f'  * Required libs for {name} : {required_libs}')
    print(f'  * Search for package {name} in {searchpath}')
    # conf.env.LIBPATH_BOOST   = [os.path.join(conf.env.BOOST, 'lib')]
    # conf.env.INCLUDES_BOOST   = [os.path.join(conf.env.BOOST,'include')]
    #boost_required_libs = ['boost_system', 'boost_program_options']
    #boost_optional_libs = ['boost_thread'] # Fix this later ...
    #for libname in required_libs:
    #    conf.check_cxx(lib = libname, use = name.upper())




    
def configure(conf):
    # Check waf version
    conf.check_waf_version(mini='1.8.8')
    # Check compiler
    if conf.options.gcc:
        conf.load('g++')
        if conf.options.mpi_on:
            conf.env['CXX'] = 'mpicxx'
    elif conf.options.icc:
        conf.load('icpc')
        if conf.options.mpi_on:
            conf.env['CXX'] = 'mpicpc'
    else:
        conf.load('compiler_cxx')
        if conf.options.mpi_on:
             conf.env['CXX'] = 'mpicxx'


    conf.env.default_search_path = conf.options.default_search_path.split()
    
    for dep in floedyn_deps:
        value = getattr(conf.options, dep,
                        conf.options.default_search_path)
        if value is None:
            if conf.options.nix_build_inputs:
                depdirs = conf.options.nix_build_inputs.split()
                for d in depdirs:
                    if d.find(dep) > -1:
                        value = d
                if value is None:
                    value = conf.env.default_search_path
            else:
                value = conf.env.default_search_path
        setattr(conf.env, dep, value)
    print("Default search path for libraries and headers of dependencies : ", conf.env.default_search_path)

    for dep in floedyn_deps:
        configure_package(conf, dep, floedyn_deps[dep], floedyn_includes[dep])
        for libname in floedyn_deps[dep]:
            conf.check_cxx(lib = libname, use = dep.upper())

    print(f'- CXX compiler is {conf.env.CXX}')
    #'boost', ['boost_system', 'boost_program_options'])
    #configure_package(conf, 'matio',['matio'])
    #configure_package(conf, 'matio',['matio'])
    
    # Boost setup
    #conf.load('find_boost', tooldir='.')
    # conf.env.BOOST = conf.env.default_search_path
    # print("BOOST PATH == ", conf.env.BOOST)
    # conf.env.LIBPATH_BOOST   = [os.path.join(conf.env.BOOST, 'lib')]
    # conf.env.INCLUDES_BOOST   = [os.path.join(conf.env.BOOST,'include')]
    # boost_required_libs = ['boost_system', 'boost_program_options']
    # #boost_optional_libs = ['boost_thread'] # Fix this later ...
    # for libname in boost_required_libs:
    #     conf.check_cxx(lib = libname, use = 'BOOST')


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
        "includes": ['../src'
                    ], #+ [path for path in os.environ["PATH"].split(":") if not "bin" in path],
        "lib": ['boost_system',
                'boost_program_options',
                'matio',
                "hdf5",
                "hdf5_cpp",
                "gmp", "mpfr", "boost_thread"
                ],
        "framework": ["Accelerate"],
        "frameworkpath" : ["/System/Library/Frameworks"]
    }
    if debug:
        OPTION_DICT.update({
            "linkflags": ['-g'],
            "cxxflags": [
                '-g',
                '-std=c++14',
                 '-O0',
                 "-Wall", #"-Wextra",
            ],
            "defines": []
        })
    else:
        OPTION_DICT.update({
            "linkflags": [],
            "cxxflags": [
                '-std=c++14',
                 "-O3",
                 # "-march=native", # g++ fails with this
                 "-mtune=native",
                 "-Wall", "-Wextra", #"-Wshadow",
                 "-Wno-unused-parameter", "-Wno-unused-local-typedef",
                 "-Wno-gnu-anonymous-struct", "-Wno-nested-anon-types", # floe/geometry/geometries/point.hpp
                 "-Wno-redeclared-class-member",
                 #"-isystem /usr/local/include/boost/",
                 # "-pedantic"
             ],
            "defines": ["NDEBUG"]
        })
    OPTION_DICT["cxxflags"].extend(os.environ.get("CFLAGS", "").split(" "))
    OPTION_DICT["linkflags"].extend(os.environ.get("LDFLAGS", "").split(" "))
    return OPTION_DICT


import subprocess

def build(bld):
    opts = get_option_dict(bld.options.debug)
    opts['use']= []
    bld.options.install_path = '${PREFIX}'
    for dep in floedyn_deps:
        opts['use'].append(dep.upper())
    if bld.options.omp:
        opts["linkflags"].append("-fopenmp")
        opts["cxxflags"].append("-fopenmp")
    if "MULTIOUTPUT" in bld.options.defmacro:
        print("compilation with multiple output files.")
        opts["defines"].append('MULTIOUTPUT')
    if "LCPSTATS" in bld.options.defmacro:
        print("compilation with LCP statistics storage.")
        opts["defines"].append('LCPSTATS')
    if "MPI" in bld.options.target:
        opts["linkflags"].extend(["-lmpi"])
        opts["defines"].append('MPIRUN')
        opts["cxxflags"].extend(subprocess.check_output(["mpicc", "--showme:compile"]).strip().split(b" "))
        opts["linkflags"].extend(subprocess.check_output(["mpicc", "--showme:link"]).strip().split(b" "))
    if bld.options.target in ["FLOE", "FLOE_PBC", "FLOE_MPI"]:
        opts["source"] = ["product/FLOE.cpp"] + recursive_file_finder("src/floe", "*.cpp")
        opts["target"] = bld.options.target
        if bld.options.target == "FLOE_PBC":
            opts["defines"].append('PBC')
        bld.program(**opts)
    elif "FLOE" in bld.options.target:
        print( "TARGET", bld.options.target )
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

