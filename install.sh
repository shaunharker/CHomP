#!/bin/bash
# install.sh [--prefix=PREFIX] [--build=BUILDTYPE]    \
#            [--search=SEARCHPATH] [CMake arguments]
#
#  Build the project with the supplied configurations,
#    or else default values.
#
#   PREFIX gives the location to install.
#   BUILDTYPE is either Debug or Release
#     (or some other CMake recognizable build type)
#   SEARCHPATH is an optional location to search for headers
#     and libraries (i.e. SEARCHPATH/include and SEARCHPATH/lib)
#   The default setting for PREFIX is /usr/local unless it is not writable
#     in which case it is ~/.local.
#   The default setting for BUILDTYPE is Release
#   The default setting for SEARCHPATH is to be equal to PREFIX
#   Additional arguments will be passed to CMake. Any paths in these arguments
#   should be absolute.

SRC_ROOT=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
build="$SRC_ROOT/.install/build.sh"
source $SRC_ROOT/.install/parse.sh
cd $SRC_ROOT
$build --prefix=$PREFIX --searchpath=$SEARCHPATH --build=$BUILDTYPE --test $MASS || exit 1
