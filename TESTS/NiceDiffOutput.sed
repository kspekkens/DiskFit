#n
# If we match the "start" string, branch forward.
/Best fitting values/b mloop
# Otherwise, delete line and restart.
d

:mloop
# Read in new line
n

# If we have found the end string, quit.
/^Degrees of freedom in fit/{
    q
}

# If line matches any number of things we don't want to print,
# ignore it and go back to start of loop.
/^-----/{
    b mloop
}
/^# points Dn used in fit/{
    b mloop
}
/^Minimization Details/{
    b mloop
}
# Skip blank lines
/^\s*$/{
    b mloop
}

# Otherwise, print string and loop
p
b mloop
