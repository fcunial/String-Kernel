<img align="right" src="./logo.png" width="395" height="160"/>

# Bwtmaw

Tools for computing minimal absent words (MAWs) and minimal rare words (MRWs) in small space. The tool allows computing several scores for M\*Ws, and to output just those with high score.

References
------------

The theory behind this code is described in the following paper:

* D. Belazzougui, and F. Cunial (2017). [A framework for space-efficient string kernels](https://link.springer.com/article/10.1007/s00453-017-0286-4). Algorithmica 79.3 (2017): 857-883.

Requirements
------------

* A modern, C++11 ready compiler such as [g++](https://gcc.gnu.org) version 4.9 or higher, or [clang](https://clang.llvm.org) version 3.2 or higher.
* A 64-bit operating system. The code was tested on both Mac OS X and Linux.

Installing and testing
------------

Compile the rest with `make`:

```
make tests
make optimized
```

The `tests` executable runs the test suite.

Huge pages
----------

The iterator is approx. 15% faster if huge pages are used. In Linux, one can list the page sizes supported by the system with:

```
$ pagesize --all
4096
2097152
1073741824
```

To use huge pages, one must allocate a pool of them. To check which huge page pools are already available, do:

```
$ hugeadm --pool-list
       Size  Minimum  Current  Maximum  Default
   2097152        0        0        0        *       
1073741824        0        0        0        
```

