# Makefile for DiskFit
# If you're not using gfortran, please modify the file to set the following variables
# according to your compiler/setup: FC FFLAGS link_flags link_libs

# Variables for implicit Fortran rule
# FC: Fortran Compiler
# FFLAGS: Flags for compiling Fortran files
FC     := gfortran

# The `-frounding-math -fsignaling-nans -mfpmath=sse -msse2` flags cause the resulting executable
# to have full IEEE compliance and more consistent behaviour across machines and CPU architectures.
# See https://gcc.gnu.org/wiki/FloatingPointMath
# -O1 enables optimizations, however upon research, no optimizations enabled will cause any
# inconsistencies in floating point operations, and O1 results in about 3-5x speedup.
# See https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
FFLAGS := -O1 -frounding-math -fsignaling-nans -mfpmath=sse -msse2 -march=x86-64 -mtune=generic
FFLAGS += -ICOMMON

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
src_modules := $(wildcard MODULES/*.f)
obj_modules := $(patsubst %.f,%.o,$(src_modules))
src_common := $(wildcard COMMON/*.f LINPACK/*.f NRCODE/*.f)
obj_common := $(patsubst %.f,%.o,$(src_common))

all: $(targets)

%: CFITSIO/libcfitsio.a %.o $(obj_modules) $(obj_common)
	$(FC) $(FFLAGS) $(link_flags) -o $@ $^ $(link_libs)

CFITSIO/libcfitsio.a:
	$(MAKE) -C CFITSIO libcfitsio.a

# Indicate that all other objects depend on the modules.
# This forces the modules to be built first.
$(obj_targets) $(obj_common): $(obj_modules)

# We use the implicit rule for compiling Fortran to object files.

.PHONY: all clean clean-all
clean:
	rm -f $(obj_targets) $(obj_modules) $(obj_common) *.mod
	$(MAKE) -C CFITSIO clean
clean-all: clean
	rm -f $(targets)
	$(MAKE) -C CFITSIO clean-all
