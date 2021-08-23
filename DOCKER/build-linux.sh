#!/usr/bin/env bash

function panic {
    echo "ERROR: $1" >&2
    exit 1
}

# Exit on errors
set -e

# Get the DiskFit base directory (containing CODE, DOCKER, etc.)
DF_BASE=$(dirname $0)/..
# If the directory doesn't exist for some reason, panic
[ -d "$DF_BASE" ] || panic "DiskFit base directory not valid."

# Read the list of targets from the DiskFit Makefile so we know what
# files to copy from the build container
TARGETS=$(sed -n -e "/targets := /{" -e "s/targets := //" -e "p" -e "q" -e "}" < "$DF_BASE/CODE/Makefile")
# If there are no targets, panic
[ -z "$TARGETS" ] && panic "Unable to enumerate targets from Makefile."

# Run the Dockerfile to actually build DiskFit
docker build -t diskfit_build -f "$DF_BASE/DOCKER/Dockerfile" "$DF_BASE"

# Spin up a new Docker container so we can copy the executable files out
DOCKER_CONTAINER_ID=$(docker create diskfit_build)

# Don't exit on errors after this to avoid leaving around a useless
# Docker container
set +e

# For each target, copy it from the container or give a warning
for x in $TARGETS ; do
    docker cp "$DOCKER_CONTAINER_ID":/build/CODE/"$x" "$DF_BASE/DOCKER/$x" \
        || echo "WARNING: Failed to copy target $x"
done

# Once all files are copied, delete the container
docker rm "$DOCKER_CONTAINER_ID" >/dev/null
