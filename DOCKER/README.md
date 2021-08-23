# Docker build system for Linux
This folder contains the Docker setup for building DiskFit inside of an
Alpine Linux Docker container. This is done since Alpine uses musl libc,
which is released under a permissive MIT license allowing us to
statically link against it.

In order to build DiskFit with Docker, follow these steps:
1. Install Docker on your build machine.
1. Download the DiskFit source code and extract it on the build machine.
1. Add the required Numerical Recipes code files under CODE/NRCODE.
1. Change into this DOCKER folder.
1. Run the following command: `bash build-linux.sh`

From there, Docker will automatically fetch the Alpine Linux image,
install the required packages (including gfortran), and build DiskFit.
Once it is done, all the DiskFit executables will be copied into the
DOCKER directory.
