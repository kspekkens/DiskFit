#!/bin/bash
function extract() {
	file="$1"
	prefix="$2"
	cat "${file}" | sed -f NiceDiffOutput.sed | sed -e "s/^/${prefix} /"
}
function gen_equals() {
	count="$1"
	i=0
	while [ "$i" -lt "$count" ]; do 
		i=$(($i+1))
		echo "================================================================================"
	done
}
if [ ! -f "$1" ] || [ ! -f "$2" ] ; then
	echo "Usage: bash $0 <file a> <file b>"
	echo "Prints out an easy to read comparison of each output parameter"
	echo "for two DiskFit runs."
	echo "Example: bash $0 <expected .out file> <actual .out file>"
	exit 1
fi
a=$(extract "$1" '<')
b=$(extract "$2" '>')
wca=$(wc -l <<< "$a")
wcb=$(wc -l <<< "$b")
wcm="$wcb"
[ "$wca" -gt "$wcb" ] && wcm="$wca"
echo "< $1"
echo "> $2"
paste -d '\n' <( gen_equals "$wcm" ) <( cat <<< "$a" ) <( cat <<< "$b" )
