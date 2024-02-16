# Nix usage for floecpp package

At the beginning of each terminal session:


```
source /applis/site/nix.sh
```

## First case : development

i.e. prepare environment and run configure/build explicitely

```
sh ./start_env_devel4nix.sh
```

This will start a shell environment with all required dependencies for floecpp installed.

Then run waf configure, build, ... as usual.

You can use the configure option --with-nix to automatically add all searched directories for dependencies:

```
./waf configure --with-nix="$buildInputs"

# check:
echo $buildInputs

```

## Build and install soft system-wide

Check list of local nix packages :

```
nix-env -f . -qaPs
```

Install:

```
nix-env -f . -iA floecpp-dev

which FLOE
/home/perignon/.nix-profile/bin/FLOE
```



