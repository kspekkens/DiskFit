#!/usr/bin/env bash

panic () {
    echo "ERROR: $1" >&2
    exit 1
}

testdir=$(dirname "$0")
echo "===== Test directory is $testdir"
cd "$testdir" || panic "Unable to open test directory $testdir"

# For each .test testcase in the current directory...
for x in *.test ; do
    echo "===== Running testcase: $x"
    # Use the other script to run the test
    bash RunSingleTest.sh "$x" || panic "Testcase $x failed!"
done
exit 0
