{ pkgs ? import <nixpkgs> {} }:
with pkgs;


let


boost-dev-meta = (pkgs.boost163.meta // {outputsToInstall =["out" "dev"]; });
mpfr-dev-meta = (pkgs.mpfr.meta // {outputsToInstall =["out" "dev"]; });
gmp-dev-meta = (pkgs.gmp.meta // {outputsToInstall =["out" "dev"]; });

in


{

inherit pkgs;


python-env = pkgs.python36.withPackages (ps: with ps; [pip  ipython ]);
boost-dev = pkgs.boost163 // {meta = boost-dev-meta;};
mpfr-dev = pkgs.mpfr // {meta = mpfr-dev-meta;};
gmp-dev = pkgs.gmp // {meta = gmp-dev-meta;};

floecpp = callPackage ./floecpp.nix {};

# Run :
# source /applis/site/nix.sh
# nix-env -f . -iA python-env boost-dev mpfr-dev
# nix-env -f "<nixpkgs>" -i cgal hdf5-cpp gmp-with-cxx matio eigen

}