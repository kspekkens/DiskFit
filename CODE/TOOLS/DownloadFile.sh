#!/bin/sh
panic () {
    echo "ERROR: $1" >&2
    exit 1
}
URL="$1"
OUTPUT="$2"
if which curl > /dev/null ; then
    echo "Downloading using curl..."
    curl -o "$OUTPUT" "$URL"
elif which wget > /dev/null ; then
    echo "Downloading using wget..."
    wget -O "$OUTPUT" "$URL"
else
    panic "Neither curl nor wget installed. Please install one of curl or wget."
fi
