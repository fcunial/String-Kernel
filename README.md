<img align="right" src="./logo.png" width="395" height="160"/>

# Bwtmaw

Tools for computing minimal absent words (MAWs) and minimal rare words (MRWs) in small space. The tool allows computing several scores for M\*Ws, and to output just those with high score.

References
------------

The theory behind this code is described in the following paper:

* D. Belazzougui, and F. Cunial (2017). [A framework for space-efficient string kernels](https://link.springer.com/article/10.1007/s00453-017-0286-4). Algorithmica 79.3 (2017): 857-883.

Requirements
------------

* A C compiler with support for OpenMP, such as [g++](https://gcc.gnu.org).
* A 64-bit operating system. The code was tested on both Mac OS X and Linux.

Building
------------

First of all, clone the git and checkout this branch
```
git clone https://github.com/fcunial/Bwtman.git
cd Bwtman
git checkout v2_saad
```

Then we need to install `libdivsufsort64`
```
git clone https://github.com/y-256/libdivsufsort.git
cd libdivsufsort
mkdir build
cd build
```
Change the line in CmakeLists.txt
```
option(BUILD_DIVSUFSORT64 "Build libdivsufsort64" OFF) 
```
to `ON` and Run `Cmake`
```
cmake -DCMAKE_BUILD_TYPE="Release"  -DCMAKE_INSTALL_PREFIX="/usr/local" ..
make
sudo make install
cd ../..
```
Then get the path of libdivsufsort64.a. If you follow the instructions above, then the path will be "/usr/local/lib/libdivsufsort64.a"

install `jansson`
```
wget http://digip.org/jansson/releases/jansson-2.13.1.tar.gz
bunzip2 -c jansson-2.13.1.tar.bz2 | tar xf -
cd jansson-2.13.1
./configure
make
make check
make install
cd ..
```
Then get the path of libjansson.a. In make file (Makefile) change the following:
```
CC="/usr/bin/gcc"
```
to be the path of your gcc compiler, use `which gcc` to find the path.

After that change the path of the libraries installed earlier as follow:
```
DIVSUFSORT_OBJS=$(LIB_PATH)/libdivsufsort64.a
JANSSON=$(LIB_PATH)/libjansson.a
```
Compiling
------------
To compile the framework use 
```
make
```
or use specific target, like
```
make buildIndex 
make run_MAWs_single 
make run_MRWs_single 
mkdir data
```

Running
------------
To run use
```
./run_MAWs_single cases/config-maws-countOnly.json
```
 the other cases can be found in cases/.

Be aware that the .fna / .fasta or .txt files should be in data/ folder.

<!-- 
Related code
---------

The following software computes MAWs as well:

* [CST-based language model](https://github.com/eehsan/cstlm): implements interpolated Markov models with Kneser-Ney smoothing using a similar setup of data structures as in this project. -->
