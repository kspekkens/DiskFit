# Makefile for DiskFit
# If you're not using gfortran, please modify the file to set the following variables
# according to your compiler/setup: FC FFLAGS link_flags link_libs

# Variables for implicit Fortran rule
# FC: Fortran Compiler
# FFLAGS: Flags for compiling Fortran files
FC     := gfortran

# The `-mfpmath=sse -msse2` flags cause the resulting executable
# to have more consistent behaviour across machines and CPU architectures.
# See https://gcc.gnu.org/wiki/FloatingPointMath
# -O1 enables optimizations, however upon research, no optimizations enabled will cause any
# inconsistencies in floating point operations, and O1 results in about 3-5x speedup.
# See https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
FFLAGS := -O1 -mfpmath=sse -msse2 -march=x86-64 -mtune=generic
FFLAGS += -ICOMMON
FFLAGS += -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter

# Variables for final linking stage
link_flags := 
link_libs  := 
# If not specified, do a static link
static_link ?= 1
ifeq ($(static_link),1)
	link_flags += -static
endif

# Target executables (must have <file>.f in the current directory)
targets := DiskFit bootlace
obj_targets := $(patsubst %,%.o,$(targets))

# We have to put the module files separately so that they are built first
src_modules  := $(wildcard MODULES/*.f)
obj_modules  := $(patsubst %.f,%.o,$(src_modules))
src_common   := $(wildcard COMMON/*.f)
obj_common   := $(patsubst %.f,%.o,$(src_common))
src_3rdparty := $(wildcard LINPACK/*.f NRCODE/*.f)
obj_3rdparty := $(patsubst %.f,%.o,$(src_3rdparty))

# We put libcfitsio here as a dependency so that it will be built first,
# before all the normal object files. It's the most complicated step,
# so we should get it out of the way first.
all: CFITSIO/libcfitsio.a $(targets)

# Pattern rule to build the targets by linking the .o file with all the
# other object files, as well as libcfitsio.a
%: %.o $(obj_modules) $(obj_common) $(obj_3rdparty) CFITSIO/libcfitsio.a
	$(FC) $(FFLAGS) $(link_flags) -o $@ $^ $(link_libs)

# We use the implicit rule for compiling Fortran to object files.

# Indicate that all other DiskFit objects depend on the modules.
# This forces the modules to be built first.
# Additionally, recompile all objects if the commons.h file changes
$(obj_targets) $(obj_common): $(obj_modules) COMMON/commons.h

# Build licfitsio.a using the makefile in the subdirectory.
CFITSIO/libcfitsio.a:
	$(MAKE) -C CFITSIO libcfitsio.a

# New target for generating assembly source from all Fortran files.
# Used only for debugging.
all-asm: $(patsubst %.o,%.S,$(obj_targets) $(obj_modules) $(obj_common) $(obj_3rdparty))
%.S: %.f
	$(FC) $(FFLAGS) -S -o $@ $^

clean:
	rm -f $(obj_targets) $(obj_modules) $(obj_common) $(obj_3rdparty) *.mod
	$(MAKE) -C CFITSIO clean
clean-all: clean
	rm -f $(targets)
	$(MAKE) -C CFITSIO clean-all

# Phony targets: targets that don't represent actual files, but instead
# they are symbolic names.
.PHONY: all all-asm clean clean-all
# Adjust suffixes to enable only the implicit rules for .f files.
.SUFFIXES:    # Delete the default suffixes
.SUFFIXES: .f # Define our suffix list (just .f)
# Prevent compilation of .f into executables directly, forcing our
# custom compilation rule to be used.
%: %.f
# Mark all intermediate files as secondary, meaning that they won't be
# deleted after compilation.
.SECONDARY:
