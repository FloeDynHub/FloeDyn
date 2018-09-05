{ stdenv, fetchurl, cmake, pkgconfig,
  pkgs ? import <nixpkgs> {},
  gcc ? pkgs.gcc ,
  pythonX ? pkgs.python36,
  boost ? pkgs.boost,
  gmp ? pkgs.gmp,
  version ? "1.0.0",
  devel_mode ? true,
  }:
  
with pkgs;
with stdenv.lib;

let
	pythonenv = pythonX.withPackages (ps: with ps; [pip numpy ipython h5py matplotlib lxml scipy pytest]);
        boost-dev-meta = (pkgs.boost.meta // {outputsToInstall =["out" "dev"]; });
        boost-dev = pkgs.boost // {meta = boost-dev-meta;};
        mpfr-dev-meta = (pkgs.mpfr.meta // {outputsToInstall =["out" "dev"]; });
        mpfr-dev = pkgs.mpfr // {meta = mpfr-dev-meta;};
	cereal = callPackage ./cereal.nix {};
in

stdenv.mkDerivation rec {
 name = "floecpp-${version}";
 
 enableParallelBuilding = true;	

 nativeBuildInputs = [
    pkgconfig
    pythonX
    gcc
    pythonenv
    pythonX.pkgs.wrapPython];
 
 buildInputs = [
    mpfr-dev
    gmp
    boost-dev
    cgal
    eigen
    gmpxx
    hdf5-cpp
    matio
    cereal
    ];
  
 hardeningDisable = [ "format" ];

 src = if devel_mode then ./.
 else ./.; # ssh keys issue
    #fetchgitPrivate "$HOME/.ssh" {
    ##url = "git@gricad-gitlab.univ-grenoble-alpes.fr:Mo_MIZ/Floe_Cpp.git";
    #rev = "a5189d1ba98e09612fd17a032131d6345308352c";
    #sha256 = "1wjghd5vbk42mb7jyjz4b0gmhw21d7zf4ihxwqh7nw9im3isa56s";
    #};

  configurePhase = ''
  
  ${pythonX.interpreter} $src/waf configure --prefix=$out --with-nix="$buildInputs"
  '';

  buildPhase = "${pythonX.interpreter} $src/waf --target FLOE";

  installPhase = ''
    mkdir $out
    ${pythonX.interpreter} $src/waf install --target=FLOE
  '';

  meta = with stdenv.lib; {
    homepage = https://gricad-gitlab.univ-grenoble-alpes.fr/Mo_MIZ/Floe_Cpp;
    description = "Ice granular model ";
    license = licenses.asl20;
    platforms = platforms.all;
  };


}
