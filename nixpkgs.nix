# The nixpkgs of both shells: the last commit of nixos-24.11. The Makefile
# and the scripts also pass it as -I nixpkgs=..., so nix-shell takes its
# own bash from it and needs no channel.
import (fetchTarball {
  url = "https://github.com/NixOS/nixpkgs/archive/50ab793786d9de88ee30ec4e4c24fb4236fc2674.tar.gz";
  sha256 = "1s2gr5rcyqvpr58vxdcb095mdhblij9bfzaximrva2243aal3dgx";
})
