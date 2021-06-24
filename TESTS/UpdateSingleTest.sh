#!/usr/bin/env bash

function panic {
    echo "ERROR: $1" >&2
    exit 1
}

# Get first command-line argument as test file
testfile="$1"
[ -f "$testfile" ] || panic "Unable to open testfile!"

# Parse some data from testfile and the corresponding input file
echo "=== Reading testfile '$testfile'..."
diskfit_input_file=$(sed -Ee "5{s/'([^']+)'/\\1/;q};d" < "$testfile")
echo "=== Using DiskFit input  file '$diskfit_input_file'..."
[ -f "$diskfit_input_file" ] || panic "DiskFit input file '$diskfit_input_file' doesn't exist!"
diskfit_output_file=$(sed -Ee "8{s/'([^']+)'\s*#?.*/\\1/;q};d" < "$diskfit_input_file")
echo "=== Using DiskFit output file '$diskfit_output_file'..."
diskfit_output_file_no_ext="${diskfit_output_file%\.out}"

# Clear output directory
diskfit_output_directory=$(dirname "$diskfit_output_file")
echo "=== Cleaning up output directory '$diskfit_output_directory'..."
mkdir -vp "$diskfit_output_directory"
rm -vf "$diskfit_output_directory"/*

# Double-check expected directory
expected_directory=$(sed -Ee "3q;d" < "$testfile")
echo "=== Using expected outputs directory '$expected_directory'..."
mkdir -vp "$expected_directory"

# Run DiskFit and save output to a temporary file in the output directory
echo "=== Running DiskFit..."
# sed -e '/^-- start input/,/^-- end input/{/^--/d; p};d' < "$testfile" | ./DiskFit 2>&1 | tee "$diskfit_output_file_no_ext".console_output
tail -n +5 "$testfile" | ../CODE/DiskFit 2>&1 | tee "$diskfit_output_file_no_ext".console_output

# Copy files to update expected
echo "=== Updating expected results..."
cp -rv "$diskfit_output_directory"/* "$expected_directory"/

# Run diff to compare output to expected output
echo "=== Comparing results..."
diff -q "$expected_directory" "$diskfit_output_directory" || panic "Test failed due to output miscompare!"
echo "Test passed!"
exit 0
