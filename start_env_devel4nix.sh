#! /bin/sh

nix-shell -E 'with import <nixpkgs> {}; callPackage ./floecpp.nix {}' --pure
