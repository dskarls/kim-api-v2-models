# This Source Code Form is subject to the terms of the Mozilla Public License,
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.

# load all basic KIM make configuration
KIM_API_BUILD_CONFIG = kim-api-v2-build-config
ifeq ($(shell $(KIM_API_BUILD_CONFIG) --version 2> /dev/null),)
  $(error $(KIM_API_BUILD_CONFIG) utility is not available.  Something is wrong with your KIM API package setup)
  endif
  include $(shell $(KIM_API_BUILD_CONFIG) --master-config)

# set model specific details
MODEL_NAME := Pair_Johnson_Fe__MO_857282754307_003
MODEL_CREATE_FUNCTION_NAME := create
MODEL_CREATE_FUNCTION_LANG := fortran

LOCALOBJ = Pair_Johnson_Fe.o

LOCALCLEAN =

# APPEND to compiler option flag lists
#FFLAGS   += -g
#CFLAGS   +=
#CXXFLAGS +=
#LDFLAGS  +=

# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.StandAloneModel
