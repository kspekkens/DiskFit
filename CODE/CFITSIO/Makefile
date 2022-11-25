cfitsio_src := https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3340.tar.gz
cfitsio_tarball := cfitsio.tar.gz
cfitsio_dir := cfitsio_src
target := libcfitsio.a

# We copy the flags from the main Makefile, but we change O1 to O2 since that's what's used
# already by cfitsio.
cflags := -O2 -mfpmath=sse -msse2 -march=x86-64 -mtune=generic

ifeq ($(platform),macos)
	cflags += -Dunix=1
endif

all: $(target)

$(target): $(cfitsio_dir)/libcfitsio.a
	cp $^ $@

$(cfitsio_dir)/libcfitsio.a: $(cfitsio_dir)/Makefile
	$(MAKE) -C $(cfitsio_dir) libcfitsio.a

$(cfitsio_dir)/Makefile: $(cfitsio_dir)/configure-patched
	cd $(cfitsio_dir) && CFLAGS="$(cflags)" ./configure-patched --enable-sse2

$(cfitsio_dir)/configure-patched: $(cfitsio_dir) clean configure.patch
	cp $(cfitsio_dir)/configure $(cfitsio_dir)/configure-patched
	patch $(cfitsio_dir)/configure-patched configure.patch

$(cfitsio_dir): $(cfitsio_tarball)
	mkdir -p "$@"
	tar -xz -C "$@" --strip-components=1 -f $^

$(cfitsio_tarball):
	wget -O "$@" "$(cfitsio_src)"

.PHONY: $(cfitsio_dir) all clean clean-all
clean:
	[ ! -e "$(cfitsio_dir)/Makefile" ] || $(MAKE) -C $(cfitsio_dir) clean
	rm -f $(cfitsio_dir)/Makefile $(cfitsio_dir)/configure-patched
	rm -f $(target) $(cfitsio_dir)/libcfitsio.a
clean-all: clean
clean-src:
	rm -rf $(cfitsio_dir)
	rm -f $(cfitsio_tarball)
