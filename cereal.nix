{ stdenv, fetchFromGitHub,
  pkgs ? import <nixpkgs> {},
}:
  
with pkgs;
with stdenv.lib;

stdenv.mkDerivation rec {

 version = "v1.2.2";
 name = "cereal-${version}";
 
 src = fetchFromGitHub {
   owner = "USCiLab";
   repo = "cereal";
   rev = version;
   sha256 = "1ckr8r03ggg5pyzg8yw40d5ssq40h5najvyqlnxc85fxxp8rnrx4";
   };

 dontBuild = true;

 installPhase = ''
    mkdir -p $out
    cp -R $src/include $out/
  '';


}
