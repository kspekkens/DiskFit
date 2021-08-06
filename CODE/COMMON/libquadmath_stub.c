// libquadmath_stub.c
// Cooper Harasyn, 2021-08-05
// This file is required since we aren't allowed to link against
// libquadmath as that would violate the LGPL that libquadmath is
// licensed under.
// Since we don't use any functionality from that library, we have
// no issues omitting it functionally. However, libgfortran tries to
// link against certain functions in that library, so we need to
// provide something for it to link against. Hence, the `lqm_panic`
// function, which just prints out an error message and then quits the
// program. This is a good response since the actual libquadmath
// function isn't linked in, and it should never be called.
#include <stdlib.h>
#include <stdio.h>
const char * const panic_str =
"==================================================================\n"
"= Call to a libquadmath function - unexpected coding error!      =\n"
"= This should never happen! Quitting...                          =\n"
"==================================================================\n";
void lqm_panic(void) {
    fputs(panic_str, stderr);
    exit(127);
}
// Add all undefined libquadmath functions here:
void quadmath_snprintf() __attribute__((weak, alias ("lqm_panic")));
void strtoflt128()       __attribute__((weak, alias ("lqm_panic")));
