{ pkgs ? import <nixpkgs> {} }:
with pkgs;


let


boost-dev-meta = (pkgs.boost.meta // {outputsToInstall =["out" "dev"]; });
mpfr-dev-meta = (pkgs.mpfr.meta // {outputsToInstall =["out" "dev"]; });

in


{

inherit pkgs;


python-env = pkgs.python36.withPackages (ps: with ps; [pip  ipython ]);
boost-dev = pkgs.boost // {meta = boost-dev-meta;};
mpfr-dev = pkgs.mpfr // {meta = mpfr-dev-meta;};

# Run :
# nix-env -f . -iA python-env boost-dev mpfr-dev
# nix-env -f "<nixpkgs>" -i cgal hdf5-cpp gmp-with-cxx matio eigen

}
