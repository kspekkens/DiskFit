#!/bin/sh
panic () {
    echo "ERROR: $1" >&2
    exit 1
}
resolve_lib () {
    lib="${1#-l}"
    oldIFS="$IFS"
    IFS=":"
    file=""
    for dir in $gfort_library_path ; do
        [ -d "$dir" ] || continue
        file="$dir/lib$lib.a" ; [ -f "$file" ] && break
        file="$dir/lib$lib.o" ; [ -f "$file" ] && break
        file=""
    done
    IFS="$oldIFS"
    if [ -f "$file" ] ; then
        resolved_lib="$file"
        return 0
    else
        resolved_lib=""
        return 1
    fi
}
fortest="./fortest.f"
fortest_exe="./fortest"
[ -f "${fortest}" ] && panic "${fortest} already exists"
[ -f "${fortest_exe}" ] && panic "${fortest_exe} already exists"

cat >"${fortest}" <<END
      end
END
gfort_library_path="$(gfortran -### -o "${fortest_exe}" "${fortest}" 2>&1 | grep "LIBRARY_PATH" | sed -e s/LIBRARY_PATH=//)"
gfort_libraries="$(gfortran -### --static -o "${fortest_exe}" "${fortest}" 2>&1 | grep "collect" | fmt -1 | grep -e '^ *-l' -e '.a$' | grep -v -e "-lcrt0" -e "-lSystem" -e "-lm" -e "-lquadmath" | fmt -999)"

lib_to_use=""
for lib in ${gfort_libraries} ; do
    if ( echo "${lib}" | grep "^-l" >/dev/null ) ; then
        resolve_lib "${lib}" || panic "Unable to resolve library ${lib}"
        lib="$resolved_lib"
    fi
    lib_to_use="${lib_to_use} ${lib}"
done
echo "${lib_to_use}" "${lib_to_use}"
rm "${fortest}"
