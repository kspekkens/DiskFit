FROM alpine:3.14
RUN apk add --no-cache musl-dev gfortran make patch coreutils
COPY CODE /build/CODE
RUN make -C /build/CODE docker_build=1 clean-all all
