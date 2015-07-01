#!/bin/bash
# cimg.sh
#   Script to install CImg-1.6.1

## Parse command line arguments to get install PREFIX
SHELL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SHELL_DIR/parse.sh

mkdir -p $PREFIX/include
wget http://downloads.sourceforge.net/project/cimg/CImg_1.6.1.zip || exit 1
unzip CImg_1.6.1.zip || exit 1
mv CImg-1.6.1/CImg.h ${PREFIX}/include/CImg.h || exit 1
