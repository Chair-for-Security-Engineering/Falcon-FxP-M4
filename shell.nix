# Native shell: compilers, MPFR/GMP (test_fpr, gen-const, the prec backend),
# the Boehm GC (prec), valgrind (make ctgrind) and rsync (the constant-time
# scripts). The Makefile enters it by itself when a library is missing.
# Cortex-M4 shell: shell-cm4.nix.
let
  pkgs = import ./nixpkgs.nix { };
in
pkgs.mkShell {
  packages = with pkgs; [
    git
    gnumake
    clang
    gcc
    rsync
  ];
  buildInputs = with pkgs; [
    mpfr
    gmp
    boehmgc
    valgrind
  ];
}
