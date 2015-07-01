#!/bin/bash
# boost.sh
#   Script to install Boost-1.58

## Parse command line arguments to get install PREFIX and MASS
SHELL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SHELL_DIR/parse.sh

## Install Boost-1.58
wget http://downloads.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.gz > /dev/null
tar xfz boost_1_58_0.tar.gz > /dev/null
cd boost_1_58_0
./bootstrap.sh --with-libraries=$MASS > /dev/null
./b2 --prefix=$PREFIX install > /dev/null
