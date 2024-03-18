

Todo: a few words about floedyn ...

For more information, please visit the [FloeDyn documentation](https://github.com/FloeDynHub/FloeDyn-documentation).

# Download and install

## Prerequisites

* python3
* C++ compiler
* gmp
* boost, boost_system, boost_program_options
 (works for 1.72, fails for 1.75, to be investigate ...)
* eigen
* matio
* hdf5 (c++)
* cgal
* mpfr
* cereal

## Get sources

```
> git clone git@github.com:FloeDynHub/FloeDyn.git
```


## Configure

Assuming that <source_dir> is  the path to the Floe_Cpp repository
Run: 

```
> cd <source_dir>
> python3 ./waf configure --prefix=<where-you-want-to-install> 
```

Results:

* Creates a build directory
* Inspects your system to check if all dependencies are satisfied
* Generates (in build) all files required to build and install the project


By default, waf looks for dependencies in /usr/. To specify another location (e.g. on clusters using guix or nix), try

```
> python3 ./waf configure --prefix=<where-you-want-to-install> --default-search-path=<other-path>
# eg --default-search-path=$GUIX_PROFILE on Luke or Dahu
```


More about waf: https://waf.io/apidocs/tutorial.html

## Build

```
> python3 ./waf --target FLOE -v
# -v is verbose mode, optional
# use --target FLOE_MPI for the parallel version
```


## Execute

```
mpirun -np 2 <path-to>/build/FLOE_MPI -i <input.h5> -t <nb_time_steps>
```

Some h5 files are available in Floe_Cpp/io/inputs/.





# Build and install on Gricad clusters

## Install all dependencies using guix inside a 'floedyn' profile:

```
source /applis/site/guix-start.sh
guix package  -m ./Floe_Cpp/manifest_floedyn.scm -p $GUIX_USER_PROFILE_DIR/floedyn
```

## Or work inside a guix environment

```
source /applis/site/guix-start.sh
guix environment  -m ./Floe_Cpp/manifest_floedyn.scm --pure
```

## Configure and build floedyn

In both case, use the same process as described above (waf configure and so on)
with the --default-search-path=$GUIX_PROFILE option.


## Execute

An example of a oar script is given in [run_simu.sh](./run_simu.sh).

Try 
```
oarsub -S "./run_simu.sh $HOME/Floe_Cpp/io/inputs/in_single_small_floe.h5 20"
```
