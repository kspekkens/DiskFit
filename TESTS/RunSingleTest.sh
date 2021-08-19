#!/usr/bin/env bash

function panic {
    echo "ERROR: $1" >&2
    exit 1
}

# Assert script was run from its directory (bash ./script)
starting_dir=$(pwd)
cd "$(dirname $0)"
script_dir=$(pwd)
[ "$starting_dir" == "$script_dir" ] || panic "Please run '$0' from its directory: bash `basename $0`"

# Get first command-line argument as test file
testfile="$1"
[ -f "$testfile" ] || panic "Unable to open testfile!"

# Copy DiskFit executable if it's not here already
diskfit="./DiskFit"
diskfit_built_executable="../CODE/DiskFit"
df_valid=0
if [ ! -e "$diskfit" ]; then
    # It's not here. Look at ../CODE/DiskFit
    if [ -e "$diskfit_built_executable" ]; then
        echo "Copying DiskFit executable from CODE directory..."
        cp "$diskfit_built_executable" "$diskfit"
    else
        echo "DiskFit doesn't exist in TESTS directory or in CODE directory."
        echo "Please compile or provide the DiskFit executable."
    fi
fi
if [ -e "$diskfit" ]; then
    # It exists
    if [ -d "$diskfit" ]; then
        # It exists but it is a directory
        echo "DiskFit folder exists in the current directory."
        echo "Please delete the folder named DiskFit, then"
        echo "re-run this script."
    elif [ -x "./DiskFit" ]; then
        # It exists, is not a directory, and is executable
        # This is the passing case
        df_valid=1
    else
        # It exists, but is either a directory or not executable.
        echo "DiskFit exists but is not executable."
        echo "Please ensure you have the correct DiskFit"
        echo "executable file, then run 'chmod +x DiskFit'"
        echo "to allow it to run."
    fi
fi
[ "$df_valid" -eq 1 ] || panic "Invalid DiskFit executable!"

# Copy EXAMPLE directory if it doesn't exist
example_dir="./EXAMPLE"
if [ ! -e "$example_dir" ]; then
    # It's not here. Look at ../EXAMPLE
    if [ -e "../EXAMPLE" ]; then
        echo "Copying EXAMPLE directory..."
        cp -r "../EXAMPLE" "$example_dir"
    fi
fi
[ -d "$example_dir" ] || panic "Example directory doesn't exist!"

# Determine what DiskFit input and output files are being used
echo "=== Reading testfile '$testfile'..."
diskfit_input_file=$(sed -e "5{" -e "s/'\([^']*\)'/\1/" -e "q" -e "}" -e "d" < "$testfile")
echo "=== Using DiskFit input  file '$diskfit_input_file'..."
[ -f "$diskfit_input_file" ] || panic "DiskFit input file '$diskfit_input_file' doesn't exist!"
diskfit_output_file=$(sed -e "8{" -e "s/'\([^']*\)'[[:space:]]*#*.*/\1/" -e "q" -e "}" -e "d" < "$diskfit_input_file")
[ -z "$diskfit_output_file" ] && panic "Unable to read output file from input file!"
echo "=== Using DiskFit output file '$diskfit_output_file'..."
diskfit_output_file_no_ext="${diskfit_output_file%\.out}"

# Clear output directory
diskfit_output_directory=$(dirname "$diskfit_output_file")
[ -z "$diskfit_output_directory" -o "$diskfit_output_directory" = "." ] && panic "Invalid output directory '$diskfit_output_directory'"
echo "=== Cleaning up output directory '$diskfit_output_directory'..."
mkdir -p "$diskfit_output_directory"
rm -f "$diskfit_output_directory"/*

# Double-check expected directory
expected_directory=$(sed -e "3q" -e "d" < "$testfile")
echo "=== Using expected outputs directory '$expected_directory'..."
[ -d "$expected_directory" ] || echo "WARNING: Expected outputs directory '$expected_directory' doesn't exist!"

# Run DiskFit and save output to a temporary file in the output directory
echo "=== Running DiskFit..."
# sed -e '/^-- start input/,/^-- end input/{/^--/d; p};d' < "$testfile" | ./DiskFit 2>&1 | tee "$diskfit_output_file_no_ext".console_output
tail -n +5 "$testfile" | "$diskfit" 2>&1 >"$diskfit_output_file_no_ext".console_output

# Run diff to compare output to expected output
echo "=== Comparing results..."
diff -q "$expected_directory" "$diskfit_output_directory" || panic "Test failed due to output miscompare!"
echo "Test passed!"
exit 0
