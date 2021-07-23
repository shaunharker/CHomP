# CHomP -- Computational Homology Project software

[![Build Status](https://travis-ci.org/shaunharker/CHomP.svg?branch=lite)](https://travis-ci.org/shaunharker/CHomP)[![Coverage Status](https://coveralls.io/repos/shaunharker/CHomP/badge.svg?branch=lite)](https://coveralls.io/r/shaunharker/CHomP?branch=lite)


## Overview
This software is a tool for computational homology for both cellular complexes and induced homology on maps. It is based on discrete Morse theory. It is used to compute Conley indices in the [Conley-Morse-Database](https://github.com/shaunharker/conley-morse-database) project.

See "LICENSE" for license details.
See "INSTALL" for more installation notes.

## Installation

### Dependencies

* cmake
* boost

### Example installation on Ubuntu 20.04 with conda installed

```
sudo apt-get install cmake
conda install boost
git clone https://github.com/shaunharker/CHomP.git
cd CHomP
git checkout lite
./install.sh # doesn't install to system, just to ./build/bin
./build/bin/chomp-simplicial-z2 ./examples/torus.simp
./build/bin/chomp-simplicial-z3 ./examples/torus.simp
```

### Example installation on macOS with homebrew

```
brew install cmake
brew install boost
git clone https://github.com/shaunharker/CHomP.git
cd CHomP
git checkout lite
./install.sh # doesn't install to system, just to ./build/bin
./build/bin/chomp-simplicial-z2 ./examples/torus.simp
./build/bin/chomp-simplicial-z3 ./examples/torus.simp
```

## More information

Documentation can be found [here](http://chomp.rutgers.edu/Projects/Databases_for_the_Global_Dynamics/software/LorentzCenterAugust2014.pdf).

Please e-mail Shaun Harker sharker81@gmail.com
about problems or feature requests.
