#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2013--2014, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#    Ellad B. Tadmor
#    Valeriu Smirichinski
#

#
# Portions Copyright (c) 2014, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Mingjian Wen
#

#
#  Release: This file is part of the the paper:
#           Interpolation effects in tabulated potentials,
#           ubmitted to Modelling Simulation Mater. Sci. Eng.
#


# load all basic KIM make configuration
KIM_API_BUILD_CONFIG = kim-api-v2-build-config
ifeq ($(shell $(KIM_API_BUILD_CONFIG) --version 2> /dev/null),)
  $(error $(KIM_API_BUILD_CONFIG) utility is not available.  Something is wrong with your KIM API package setup)
endif
include $(shell $(KIM_API_BUILD_CONFIG) --master-config)

# set model specific details
MODEL_NAME := Pair_Morse_Modified_MacDonaldMacDonald_Cu__MO_034823476734_001
MODEL_CREATE_FUNCTION_NAME := model_create
MODEL_CREATE_FUNCTION_LANG := c

LOCALOBJ = Pair_Morse_Modified_MacDonaldMacDonald_Cu.o

LOCALCLEAN =

CFLAGS   += -std=c99 # just to allow c++ style comments (//)


# load remaining KIM make configuration
include $(KIM_DIR)/$(builddir)/Makefile.StandAloneModel
