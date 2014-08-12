#!/bin/bash
CUR_DIR=`pwd`
if [ $# -ge 1 ]; then
    # The user supplied an argument
    PREFIX=${1}
    # Get absolute path name of install directory
    mkdir -p "${PREFIX}" 2> /dev/null
    cd "${PREFIX}" > /dev/null 2>&1
    if [ $? != 0 ] ; then
        echo "ERROR: directory '${PREFIX}' does not exist nor could be created."
        echo "Please choose another directory."
        exit 1
    else
        PREFIX=`pwd -P`
    fi
    echo "CHomP will be installed in '${PREFIX}'"    
    ARGUMENT=-DCMAKE_INSTALL_PREFIX=${CHOMP_INSTALL_PREFIX}
fi

cd ${CUR_DIR}
rm -rf build
mkdir build
cd build
cmake $ARGUMENT ..
make
make install
