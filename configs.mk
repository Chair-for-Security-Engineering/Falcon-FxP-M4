# The paper's configurations, for make CONFIG=<name>: the native backend,
# the pqm4 implementation, and the knobs that differ from the defaults
# in knobs.h. Works for native builds, 'make precision' and pqm4.

CONFIG_ieee-base                  = ieee     ieee-m4f
CONFIG_ieee-base-fftnew           = ieee     ieee-m4f    FFT_NEW=1
CONFIG_ieee-base-origshake        = ieee     ieee-m4f    SHAKE_HUANG=0
CONFIG_ieee-split                 = ieee     ieee-m4f    SPLIT=1
CONFIG_ieee-rec                   = ieee     ieee-m4f    SPLIT=1 SYMP_ROOT=1
CONFIG_ieee-rec-fftnew            = ieee     ieee-m4f    SPLIT=1 SYMP_ROOT=1 FFT_NEW=1
CONFIG_ieee-c                     = ieee     ieee-c
CONFIG_tw-m4f                     = tw       tw-m4f
CONFIG_tw-c                       = tw       tw-c
CONFIG_fxp-split                  = fxp32.64 fxp3264-cm4 TIES_UP=0
CONFIG_fxp-unfused                = fxp32.64 fxp3264-cm4 RECIP_FUSED=0 TIES_UP=0
CONFIG_fxp-rec                    = fxp32.64 fxp3264-cm4
CONFIG_fxp-rec-fftnew             = fxp32.64 fxp3264-cm4 FFT_NEW=1
CONFIG_fxp-nonsplit               = fxp32.64 fxp3264-cm4 SPLIT=0 TIES_UP=0
CONFIG_fxp-nonsplit-tiesup        = fxp32.64 fxp3264-cm4 SPLIT=0
CONFIG_fxp-nonsplit-tiesup-fftnew = fxp32.64 fxp3264-cm4 SPLIT=0 FFT_NEW=1
CONFIG_fxp-c                      = fxp32.64 fxp3264-c   L10_64=0 TIES_UP=0 TWIDDLE_U64=0
