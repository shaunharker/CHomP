#!/bin/bash
## Parse command line arguments
## !! --> Also, will mkdir the install location if it does not exist!

for i in "$@"; do case $i in
    -p=*|--prefix=*) PREFIX="${i#*=}";         shift;;
    -b=*|--build=*) BUILDTYPE="${i#*=}";       shift;;
    -s=*|--searchpath=*) SEARCHPATH="${i#*=}"; shift;;
    -t|--test) TEST="YES";                     shift;;
    *)               MASS+="${i#*=} ";         shift;;
esac; done
absolute() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
if [ ! -z "$PREFIX" ]; then 
  mkdir -p $PREFIX || ( echo Permission to mkdir $PREFIX denied && exit 1 ); 
  PREFIX=$(absolute $PREFIX); 
elif [ -w /usr/local ]; then 
  PREFIX=/usr/local; 
elif [ -d ~/.local ]; then
  PREFIX=$(absolute ~/.local ); 
else
  echo Run with admin privileges or choose a non-system install path with --prefix && exit 1; 
fi;
if [ ! -w $PREFIX ]; then echo Permission to write to $PREFIX denied && exit 1; fi
if [ -z "$BUILDTYPE" ]; then BUILDTYPE=Release; fi
if [ -z "$SEARCHPATH" ]; then SEARCHPATH=$PREFIX; fi
if [ -d "$SEARCHPATH" ]; then SEARCHPATH=$(absolute $SEARCHPATH); else exit 1; fi
