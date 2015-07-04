#!/bin/bash
# build.sh [--prefix=PREFIX] [--build=BUILDTYPE]    \
#            [--search=SEARCHPATH] [--test] [CMake arguments]
#  
#  Build the project with the supplied configurations,
#    or else default values.
#
#   PREFIX gives the location to install.
#   BUILDTYPE is either Debug or Release 
#     (or some other CMake recognizable build type)
#   SEARCHPATH is an optional location to search for headers 
#     and libraries (i.e. SEARCHPATH/include and SEARCHPATH/lib)
#   If --test is supplied then tests will be built.
#   The default setting for PREFIX is /usr/local unless it is not writable
#     in which case it is ~/.local.
#   The default setting for BUILDTYPE is Release
#   The default setting for SEARCHPATH is to be equal to PREFIX
#   Additional arguments will be passed to CMake. Any paths in these arguments
#   should be absolute.

## Parse command line arguments to get
#  PREFIX, SEARCHPATH, BUILDTYPE, TESTS, and MASS
SHELL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SHELL_DIR/parse.sh

## Set CMake arguments
CMAKE_ARGS+=" -DCMAKE_INSTALL_PREFIX=${PREFIX}"
CMAKE_ARGS+=" -DUSER_INCLUDE_PATH=${SEARCHPATH}/include"
CMAKE_ARGS+=" -DUSER_LIBRARY_PATH=${SEARCHPATH}/lib"
CMAKE_ARGS+=" -DCMAKE_BUILD_TYPE=$BUILDTYPE"
CMAKE_ARGS+=$MASS

## Build 
rm -rf build && mkdir build && cd build || exit 1
cmake $CMAKE_ARGS .. && make            || exit 1

## Test
if [[ "$TEST" == "YES" ]]; then
  make test
  if [ ! $? -eq 0 ]; then
    cat Testing/Temporary/LastTest.log
    exit 1
  fi
fi

## Install
make install || exit 1
