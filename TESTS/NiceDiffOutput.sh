#!/bin/bash
function extract() {
	file="$1"
	prefix="$2"
	cat "${file}" | sed -n -f <( cat <<-SEDSCR
		# If we match the "start" string, branch forward.
		/Best fitting values/b mloop
		# Otherwise, delete line and restart.
		d

		:mloop
		# Read in new line
		n

		# If we have found the end string, quit.
		/^Degrees of freedom in fit/{q}

		# If line matches any number of things we don't want to print,
		# ignore it and go back to start of loop.
		/^-----/{ z ; b mloop }
		/^# points Dn used in fit/{ z ; b mloop }
		/^Minimization Details/{ z ; b mloop }
		# Skip blank lines
		/^\s*$/{ z ; b mloop }

		# Otherwise, print string and loop
		s|^|${prefix} |
		p
		b mloop
SEDSCR
	)
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
