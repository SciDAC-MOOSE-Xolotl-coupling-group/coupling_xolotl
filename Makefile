###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif
# Try to use PETSc submodule if PETSC_DIR is not set
PETSC_DIR          ?=$(MOOSE_DIR)/petsc
PETSC_ARCH         ?=arch-moose

# Kokkos
KOKKOS_DIR         ?= $(CURDIR)/kokkos

# Xolotl
XOLOTL_DIR         ?= $(CURDIR)/xolotl

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework

ADDITIONAL_SRC_DEPS := $(XOLOTL_DIR)/install/include/interface.h

include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

# Darwin
ifneq (,$(findstring darwin,$(libmesh_HOST)))
	lib_suffix := dylib
else
	lib_suffix := so
endif

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).

ALL_MODULES         := no

CHEMICAL_REACTIONS  := no
CONTACT             := no
FLUID_PROPERTIES    := no
HEAT_CONDUCTION     := no
MISC                := no
NAVIER_STOKES       := no
PHASE_FIELD         := no
RDG                 := no
RICHARDS            := no
SOLID_MECHANICS     := no
STOCHASTIC_TOOLS    := no
TENSOR_MECHANICS    := no
XFEM                := no
POROUS_FLOW         := no

include $(MOOSE_DIR)/modules/modules.mk

# List XOLOTL as a dependency
# Use ADDITIONAL flags to link XOLOTL
XOLOTL_DEPEND_LIBS     := $(XOLOTL_DIR)/install/lib/libxolotlInterface.$(lib_suffix)
# -Wl,-rpath trikcy is used for load XOLOTL properly from executable
ADDITIONAL_LIBS        += -L$(XOLOTL_DIR)/install/lib -Wl,-rpath,$(XOLOTL_DIR)/install/lib -lxolotlInterface
ADDITIONAL_INCLUDES    += -I$(XOLOTL_DIR)/install/include

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := coupling_xolotl
BUILD_EXEC         := yes
GEN_REVISION       := no
# DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here

$(ADDITIONAL_SRC_DEPS): $(XOLOTL_DEPEND_LIBS)

# TODO: should list all source files as a dependency
# Then if source codes change, make will try to call "cmake"
# "cmake" should build lib with updated source codes
$(XOLOTL_DEPEND_LIBS): $(XOLOTL_DIR)/xolotl/solver/src/Solver.cpp
	cd kokkos; \
	mkdir build; \
	cd build; \
	cmake -DCMAKE_INSTALL_PREFIX=$(KOKKOS_DIR)/install -DKokkos_ENABLE_SERIAL=ON \
	-DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DKokkos_ENABLE_TESTS=OFF -DKokkos_CXX_STANDARD=17 .. ; \
	make install; \
	cd ../../xolotl; \
	mkdir build; \
	cd build; \
	cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_DIR=$(KOKKOS_DIR)/install -DPETSC_DIR=$(PETSC_DIR) \
	-DPETSC_ARCH=$(PETSC_ARCH) -DCMAKE_CXX_COMPILER=mpicxx \
	-DHDF5_ROOT=$(PETSC_DIR)/$(PETSC_ARCH) -DBOOST_ROOT=$(PETSC_DIR)/$(PETSC_ARCH) \
	-DBUILD_SHARED_LIBS=yes -DCMAKE_CXX_FLAGS_RELEASE="-O3 -fPIC" -DBUILD_TESTING=OFF \
	-DCMAKE_INSTALL_PREFIX=$(XOLOTL_DIR)/install ..; \
	make; make install \
