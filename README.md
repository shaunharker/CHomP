# CHomP -- Computational Homology Project software

[![Build Status](https://travis-ci.org/shaunharker/CHomP.svg?branch=master)](https://travis-ci.org/shaunharker/CHomP)[![Coverage Status](https://coveralls.io/repos/shaunharker/CHomP/badge.svg?branch=master)](https://coveralls.io/r/shaunharker/CHomP?branch=master)



## Overview 
This software is a tool for computational homology for both cellular complexes and induced homology on maps. It is based on discrete Morse theory. It is used to compute Conley indices in the [Conley-Morse-Database](https://github.com/shaunharker/conley-morse-database) project.

See "LICENSE" for license details.
See "INSTALL" for more installation notes.


## Installation

### Dependencies

* cmake
* boost
* GMP
* CImg
* X11

### General Instructions

* First install dependencies.
* Then, either:

#### Option 1. Install in default location.

To install into the system, i.e. `/usr/local` install with

    ./install.sh

In the event of an error (e.g. the installer says `Run with admin privileges or choose a non-system install path with —prefix`) then see the troubleshooting steps below.

#### Option 2. Install in a custom location.

Install via

    ./install.sh --prefix=/path/to/install/folder

### Specific instructions for macOS:

For macOS the easiest way to install dependencies is with "homebrew" <http://brew.sh>.

Install homebrew via

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Install dependencies via

```
brew install cmake
brew install gmp
brew install boost
brew install cimg
```

To install X11, follow the instructions on <https://www.xquartz.org>

To install CHomP:

```
git clone https://github.com/shaunharker/CHomP.git
cd CHomP
./install.sh
```

To test an example:

```
chomp-simplicial ./examples/simplex.simp
```


## Examples

### Cubical complexes

**Example:**

    chomp-cubical ./examples/square.cub

File format: a list of new-line separated d-tuples, where d is the dimension of the cubical complex, and the tuples represent d-dimensional cubes which will be added to the complex (along with all of their faces):

Example: _(Annulus-like shape made out of squares)_
```
(0, 1)
(0, 2)
(1, 0)
(1, 2)
(2, 0)
(2, 1)
(2, 2)
```

### Simplicial complexes

    chomp-simplicial ./examples/simplex.simp

File format: each line represents a simplex in the complex given by its vertices. For each simplex encountered, the simplex and all of its faces are generated. For instance `0 1 2` constructs a 2-cell, 3 1-cells, and 3 0-cells. The line `4` would create a single 0-cell corresponding to the vertex 4. Because subfaces are automatically created, it is sufficient to include maximal simplices (though allowed to include all)

Example: _(Tetrahedron along with an isolated vertex)_
```
0 1 2
0 1 3
0 2 3
1 2 3
4
```


### General chain complexes (given by matrices for boundary map)

    chomp-matrix ./examples/torus.mat

File format: The format is in terms of a sparse matrix representation of the boundary maps $d_n : C_{n+1} \to C_{n}$ for $n = 0 \cdots D-1$. The file has $D$ sections, corresponding to each map $d_0$, $d_1$, in turn. Each section begins with a line `n`, where $n$ is the index of the boundary map. Then follow lines indicating the non-zero entries for the incidence matrix for $d_n$. In particular each line is a triple `i j k` indicating that $d_{n}_{i,j} = k$, i.e. the jth cell in $C_{n+1}$ has in its boundary $k$ times the cell $i$ in $C_{n}$. It is allowed to include zero entries in the matrix explicitly, and it necessary for the software to learn of any cells which are not in the image of a boundary operator. 

Example: _(Torus; 1 2d-cell, 2 1d-cells, and 1 1d cell)_

```
0
0 0 0
0 1 0
1
0 0 0
1 0 0
```

## Troubleshooting

**Question.** The installer fails with the message

    Run with admin privileges or choose a non-system install path with —prefix

What do I do?

**Answer.** The default install location is `/usr/local` and typing `./install.sh` is equivalent to typing `./install.sh --prefix=/usr/local`. The error indicates that `/usr/local` is owned by the system and is not writable. Thus it doesn't do the install unless you give it a different location to install or otherwise change the situation.

Three explicit solutions are given below. The first two are for macOS users.

**Solution 1.** Install `homebrew`.

On macOS the package manager `homebrew` is widely used and sets the permissions on `/usr/local` to be writable. `CHomP` expects this is installed (as it is the simplest way to install the dependencies) and has made `/usr/local` writable. Therefore, if `homebrew` is not installed then installing it will likely solve the problem as it will fix the permissions on `/usr/local`.

To check if `homebrew` is installed, type

    brew help

If this says `command not found` then `homebrew` is not installed. Instructions for installing homebrew can be found on <http://brew.sh>. Currently, it is as simple as running the following line on the terminal and following instructions:

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

**Solution 2.** Change the permissions on `/usr/local`

On macOS permissions on `/usr/local` are typically made writable by user (thus preventing the above problem) when `homebrew` is installed.

It is possible that an OS update has mucked up the permissions on `/usr/local` subsequent to `homebrew` being installed. In particular the El Capitan update was known to do this.

Thus, the permissions on `/usr/local` are not suitable and need to be fixed. A known solution to this problem (<http://stackoverflow.com/questions/16432071/how-to-fix-homebrew-permissions>) is to type the following lines on the command line:

    sudo chown -R "$USER":admin /usr/local
    sudo chown -R "$USER":admin /Library/Caches/Homebrew

**Solution 3.** Make a custom installation which is not to /usr/local

_Note. This solution is not recommended to macOS users, but it IS recommended for Linux users who do not have the privileges to change the permissions on `/usr/local` or otherwise want a custom install location._

To install to a custom location, use

    ./install.sh --prefix=/path/to/where/you/want/it/installed

This will install the executables to

    /path/to/where/you/want/it/installed/bin

Now update `~/.bashrc` to add this location to the `PATH` environmental variable, so that the command line can find the executables (e.g. `chomp-simplicial`). That is, add the line

    export PATH=/path/to/where/you/want/it/installed/bin/:$PATH

to the end of your `~/.bashrc` file (create one if it does not exist).

## More information

Documentation can be found [here](http://chomp.rutgers.edu/Projects/Databases_for_the_Global_Dynamics/software/LorentzCenterAugust2014.pdf).

Please e-mail Shaun Harker sharker81@gmail.com
about problems or feature requests.
