# Cortex-M4 shell for the pqm4 builds and board runs (NUCLEO-L4R5ZI). The
# pqm4 targets of the Makefile enter it by themselves. mkShell adds the
# host gcc and gcov, which scripts/ct_cover.py uses next to the Arm tools.
#
# nixpkgs 24.11 (nixpkgs.nix) gives binutils 2.42; 2.44 rejects the TWFalcon assembly
# (functions without '.type sym, %function').
let
  pkgs = import ./nixpkgs.nix { };
in
pkgs.mkShell {
  packages = with pkgs; [
    git
    gnumake
    gcc-arm-embedded-13
    stlink
    openocd
    rsync
    (python3.withPackages (ps: [ ps.pyserial ps.tqdm ]))
  ];
}
