# Makefile for DiskFit
# If you're not using gfortran, please modify the file to set the following variables
# according to your compiler/setup: FC FFLAGS link_flags link_libs

# Variables for implicit Fortran rule
# FC: Fortran Compiler
# FFLAGS: Flags for compiling Fortran files
FC := gfortran
CC := gcc

# The `-mfpmath=sse -msse2` flags cause the resulting executable
# to have more consistent behaviour across machines and CPU architectures.
# See https://gcc.gnu.org/wiki/FloatingPointMath
# -O1 enables optimizations, however upon research, no optimizations enabled will cause any
# inconsistencies in floating point operations, and O1 results in about 3-5x speedup.
# See https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
FFLAGS := -O1 -mfpmath=sse -msse2 -march=x86-64 -mtune=generic
FFLAGS += -ICOMMON
FFLAGS += -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter

CFLAGS := -O1 -mfpmath=sse -msse2 -march=x86-64 -mtune=generic
CFLAGS += -ICOMMON
CFLAGS += -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter

# Variables for final linking stage
linker := $(FC)
link_flags := 
link_libs  := 
# If not specified, do a static link
static_link ?= 1

platform := unknown
ifeq ($(OS),Windows_NT)
	# It's unknown how to treat this - let the user figure it out if they want
	# to compile for Windows
	platform := win
else
	# Run `uname` to check what OS this is
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		platform := linux
		# We don't want to static link with libquadmath, so instead of
		# using gfortran to link, we use gcc and explicitly link with
		# libgcc, libfortran, and libm.
		ifeq ($(static_link),1)
			linker := $(FC)
			link_flags += -nodefaultlibs -static
			link_libs += -lgfortran -lm -lc -lgcc -lgcc_eh
		endif
	endif
	ifeq ($(UNAME_S),Darwin)
		platform := macos
		# We can't static link on MacOS just by doing `-static` since
		# Apple has made it impossible to statically link against their
		# libSystem. We still do a dynamic link, but we statically link
		# against all other libraries that we need.
		# We need to manually locate the libraries to link statically
		# and specify the archive or object files on the command line,
		# instead of the `-l` option, with the exception of -lSystem, which we
		# dynamically link with.
		ifeq ($(static_link),1)
			linker := $(FC)
			link_flags += -nodefaultlibs -static-libgcc
			link_libs += $(shell sh TOOLS/DarwinFindLibraryArchive.sh) -lSystem
		endif
	endif
endif

# Target executables (must have <file>.f in the current directory)
targets := DiskFit bootlace makemap set_mask
obj_targets := $(patsubst %,%.o,$(targets))

# List of required Numerical Recipes code files - for printing error msg
nrcode_files_required_msg := brent.f f1dim.f hpsort.f linmin.f mnbrak.f
nrcode_files_required_msg += moment.f powell.f ran1.f select.f 
src_nrcode_required_msg := $(patsubst %,NRCODE/%,$(nrcode_files_required_msg))
nrcode_patch_required_msg := brent.f f1dim.f linmin.f mnbrak.f
nrcode_patch_required_msg += moment.f powell.f ran1.f
nrcode_patch_required_msg_with_suffix := $(patsubst %,%.dblpatch,$(nrcode_patch_required_msg))
src_nrcode_patch_required_msg := $(patsubst %,NRCODE/%,$(nrcode_patch_required_msg_with_suffix))

# Numerical Recipes code files used for compiling
nrcode_files_dbl    := brent_dbl.f f1dim_dbl.f linmin_dbl.f mnbrak_dbl.f
nrcode_files_dbl    += moment_dbl.f powell_dbl.f ran1_dbl.f
nrcode_files_no_dbl := hpsort.f select.f
src_nrcode_dbl      := $(patsubst %,NRCODE/%,$(nrcode_files_dbl))
src_nrcode_no_dbl   := $(patsubst %,NRCODE/%,$(nrcode_files_no_dbl))
src_nrcode          := $(src_nrcode_dbl) $(src_nrcode_no_dbl)

# We have to put the module files separately so that they are built first
src_modules    := $(wildcard MODULES/*.f)
obj_modules    := $(patsubst %.f,%.o,$(src_modules))
src_common     := $(wildcard COMMON/*.f)
obj_common     := $(patsubst %.f,%.o,$(src_common))
src_common_c   := $(wildcard COMMON/*.c)
obj_common_c   := $(patsubst %.c,%.o,$(src_common_c))
src_3rdparty   := $(wildcard LINPACK/*.f Burkardt/*.f fftpak/*.f)
src_3rdparty   += $(src_nrcode)
obj_3rdparty   := $(patsubst %.f,%.o,$(src_3rdparty))
obj_nontargets := $(obj_modules) $(obj_common) $(obj_common_c) $(obj_3rdparty)

# We put libcfitsio here as a dependency so that it will be built first,
# before all the normal object files. It's the most complicated step,
# so we should get it out of the way first.
all: $(targets)

# Pattern rule to build the targets by linking the .o file with all the
# other object files, as well as libcfitsio.a
$(targets): %: %.o $(obj_nontargets) CFITSIO/libcfitsio.a
	$(linker) $(link_flags) -o $@ $^ $(link_libs)

# We use the implicit rule for compiling Fortran to object files.

# Indicate that all other DiskFit objects depend on the modules.
# This forces the modules to be built first.
# Additionally, recompile all objects if the commons.h file changes
$(obj_targets) $(obj_common): $(obj_modules) COMMON/commons.h

# Indicate that NRCODE _dbl files can be built using a patch file.
$(src_nrcode_dbl): %_dbl.f: %.f | %.f.dblpatch
	patch -o "$@" "$<" "$<.dblpatch" || (rm "$@" && false)

$(src_nrcode_required_msg):
	@echo "ERROR: Please ensure the following Numerical Recipes code files"
	@echo "       are provided in the CODE/NRCODE subdirectory:"
	@for x in $(nrcode_files_required_msg) ; do echo "           $$x" ; done
	@false

$(src_nrcode_patch_required_msg):
	@echo "ERROR: Please ensure the following Numerical Recipes code patches"
	@echo "       are provided in the CODE/NRCODE subdirectory:"
	@for x in $(nrcode_patch_required_msg_with_suffix) ; do echo "           $$x" ; done
	@false


# Build libcfitsio.a using the makefile in the subdirectory.
CFITSIO/libcfitsio.a:
	$(MAKE) -C CFITSIO libcfitsio.a

# New target for generating assembly source from all Fortran files.
# Used only for debugging.
all-asm: $(patsubst %.o,%.S,$(obj_targets) $(obj_nontargets))
%.S: %.f
	$(FC) $(FFLAGS) -S -o $@ $^

clean:
	rm -f $(obj_targets) $(obj_nontargets) *.mod
clean-all: clean
	rm -f $(targets)
	$(MAKE) -C CFITSIO clean-all
clean-src:
	$(MAKE) -C CFITSIO clean-src
clean-asm:
	rm -f $(patsubst %.o,%.S,$(obj_targets) $(obj_nontargets))

# Phony targets: targets that don't represent actual files, but instead
# they are symbolic names.
.PHONY: all all-asm clean clean-all clean-src clean-asm
# Adjust suffixes to enable only the implicit rules for .f files.
.SUFFIXES:          # Delete the default suffixes
.SUFFIXES: .o .f .c # Define our suffix list (build .f/.c into .o)
# Prevent compilation of .f or .o into executables directly, forcing our
# custom compilation rule to be used.
%: %.f
%: %.o
# Mark all intermediate files as secondary, meaning that they won't be
# deleted after compilation.
.SECONDARY:
