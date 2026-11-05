# SPDX-License-Identifier: Apache-2.0 or CC0-1.0
ifndef _CONFIG
_CONFIG :=

ifndef SRCDIR
SRCDIR := $(CURDIR)
endif

##############################
# Include retained variables #
##############################

RETAINED_VARS :=

CONFIG := obj/.config.mk

-include $(CONFIG)

$(CONFIG):
	@echo "  GEN     $@"
	$(Q)[ -d $(@D) ] || mkdir -p $(@D)
	@echo "# These variables are retained and can't be changed without a clean" > $@
	@$(foreach var,$(RETAINED_VARS),echo "$(var) := $($(var))" >> $@; echo "LAST_$(var) := $($(var))" >> $@;)

###############
# Some Macros #
###############
objs = $(addprefix obj/,$(addsuffix .o,$(1)))
hashprofobjs = $(addprefix obj/hashprof/,$(addsuffix .o,$(1)))
hostobjs = $(addprefix obj-host/,$(addsuffix .o,$(1)))

Q ?= @

PLATFORM ?=

ifeq (,$(PLATFORM))
$(error No PLATFORM specified (see README.md for a list of supported platforms)!)
endif

ifeq (,$(wildcard $(SRCDIR)/mk/$(PLATFORM).mk))
$(error Unknown platform!)
endif

# The platform file should set all necessary CC/CXX/AR/LD variables
include mk/$(PLATFORM).mk

RETAINED_VARS += PLATFORM

Q ?= @

HOST_CC := cc
HOST_AR := ar
HOST_LD := $(HOST_CC)

################
# Common Flags #
################

SYSROOT ?= $(shell $(CC) --print-sysroot)

CFLAGS += \
	-std=gnu99 \
	--sysroot=$(SYSROOT) \
	-I$(SRCDIR)/common \
	-I$(SRCDIR)/mupq/common

DEBUG ?=
OPT_SIZE ?=
LTO ?=
AIO ?= 1
MUPQ_ITERATIONS ?= 1
# [Falcon_FxP_M4] MUPQ_EXTRA_CPPFLAGS: extra preprocessor defines
# forwarded from the repository root Makefile (the SPLIT/SYMP_ROOT/
# SYMP_TREE/GAMMA/GAMMA_NTT variant knobs and EXTRA_CFLAGS). Retained
# so a change forces a rebuild like the other build variables.
MUPQ_EXTRA_CPPFLAGS ?=
# [Falcon_FxP_M4] MUPQ_DEFAULT_DEF: an implementation config.mk's
# overridable baked default -- emits -D$(1)=$(2) unless the macro is
# already set through MUPQ_EXTRA_CPPFLAGS (an explicit variant knob on
# the root make command line), which then takes precedence.
MUPQ_DEFAULT_DEF = $(if $(findstring $(1)=,$(MUPQ_EXTRA_CPPFLAGS)),,-D$(1)=$(2))

RETAINED_VARS += DEBUG OPT_SIZE LTO AIO MUPQ_ITERATIONS MUPQ_EXTRA_CPPFLAGS

ifeq ($(DEBUG),1)
CFLAGS += -O0 -g3
else ifeq ($(OPT_SIZE),1)
CFLAGS += -Os -g3
else
CFLAGS += -O3 -g3
endif

ifeq ($(LTO),1)
CFLAGS += -flto
LDFLAGS += -flto
endif

CPPFLAGS += -DMUPQ -DMUPQ_NAMESPACE=$(MUPQ_NAMESPACE) -DMUPQ_ITERATIONS=$(MUPQ_ITERATIONS)
# [Falcon_FxP_M4] append the forwarded defines
CPPFLAGS += $(MUPQ_EXTRA_CPPFLAGS)

CFLAGS += \
	-Wall -Wextra -Wshadow \
	-MMD \
	-fno-common \
	$(CPPFLAGS)

LDFLAGS += \
	-Lobj \
	-Wl,--gc-sections

LDLIBS += -lm

HOST_CPPFLAGS += -DMUPQ_NAMESPACE=$(MUPQ_NAMESPACE)

HOST_CFLAGS += \
	-O2 -g3 \
	-I$(SRCDIR)/mupq/common \
	-Wall -Wextra -Wshadow \
	-MMD \
	-fno-common \
	$(HOST_CPPFLAGS)

HOST_LDFLAGS += \
	-Lobj-host

# Check if the retained variables have been changed
define VAR_CHECK
ifneq ($$(origin LAST_$(1)),undefined)
ifneq "$$($(1))" "$$(LAST_$(1))"
$$(info Variable $(1) changed, forcing rebuild!)
.PHONY: $(CONFIG)
endif
endif
endef

$(foreach VAR,$(RETAINED_VARS),$(eval $(call VAR_CHECK,$(VAR))))
endif
