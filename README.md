<img align="right" src="./logo.png" width="395" height="197"/>

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

Related code
---------

The following software computes MAWs as well:

<!-- * [CST-based language model](https://github.com/eehsan/cstlm): implements interpolated Markov models with Kneser-Ney smoothing using a similar setup of data structures as in this project. -->
