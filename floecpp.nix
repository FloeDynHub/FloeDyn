{ stdenv, fetchurl, cmake, pkgconfig,
  pkgs ? import <nixpkgs> {},
  gcc ? pkgs.gcc ,
  pythonX ? pkgs.python36,
  boost ? pkgs.boost,
  gmp ? pkgs.gmp,
  version ? "1.0.0",
  }:
  
with pkgs;
with stdenv.lib;

let
	pythonenv = pythonX.withPackages (ps: with ps; [pip numpy ipython h5py matplotlib lxml scipy pytest]);
        boost-dev-meta = (pkgs.boost.meta // {outputsToInstall =["out" "dev"]; });
        boost-dev = pkgs.boost // {meta = boost-dev-meta;};
        mpfr-dev-meta = (pkgs.mpfr.meta // {outputsToInstall =["out" "dev"]; });
        mpfr-dev = pkgs.mpfr // {meta = mpfr-dev-meta;};
in

stdenv.mkDerivation rec {
 name = "floecpp-${version}";
 
 enableParallelBuilding = true;	

 nativeBuildInputs = [
    pkgconfig
    pythonX
    gcc
    mpfr-dev
    pythonX.pkgs.wrapPython];
 
 buildInputs = [
    gmp
    boost-dev
    cgal
    eigen
    gmpxx
    hdf5-cpp
    matio
    ];
  
 hardeningDisable = [ "format" ];
 src = ./.;

  configurePhase = ''
  ls $cgal
  ${python.interpreter} waf configure --prefix=$out --default-search-path=$buildInputs
  '';

  buildPhase = "${python.interpreter} waf --target FLOE";

  installPhase = ''
    #${python.interpreter} waf install
  '';

  meta = with stdenv.lib; {
    homepage = https://gricad-gitlab.univ-grenoble-alpes.fr/Mo_MIZ/Floe_Cpp;
    description = "Ice granular model ";
    license = licenses.asl20;
    platforms = platforms.all;
  };


}
