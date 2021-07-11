<img align="right" src="./logo.png" width="395" height="160"/>

# Bwtmaw

Tools for computing minimal absent words (MAWs) and minimal rare words (MRWs) in small space. The tool allows computing several scores for M\*Ws, and to output just those with high score.

## References

The theory behind this code is described in the following paper:

* D. Belazzougui, and F. Cunial (2017). [A framework for space-efficient string kernels](https://link.springer.com/article/10.1007/s00453-017-0286-4). Algorithmica 79.3 (2017): 857-883.

## Requirements

* A modern, C++11 ready compiler such as [g++](https://gcc.gnu.org) version 4.9 or higher, or [clang](https://clang.llvm.org) version 3.2 or higher.
* A 64-bit operating system. The code was tested on both Mac OS X and Linux.

## Installing and testing

Compile the rest with `make`:

```
make tests
make optimized
```

The `tests` executable runs the test suite.


## Huge pages

Like most data structures built on the BWT, the iterator accesses the underlying index without any spatial or temporal locality. This means that most accesses induce a cache miss and, symmetrically, a TLB miss. The TLB is typically small and the BWT index is typically large: if one uses small memory pages, just a few page translation entries can fit in the TLB, and every query to the TLB is very likely a miss. Conversely, if one uses large memory pages, a large fraction of the page translation entries needed by the iterator can fit into the TLB, and misses become rare. If the TLB is implemented as part of one or more levels of CPU cache, fewer translation entries might also mean more cache space for data. The time spent serving a TLB miss can be a significant fraction of the total runtime of the iterator, since the iterator is memory-bound.

Running the count-only variant of the one-string iterator with 2MB pages instead of the default 4KB pages gives a speedup of approx 15% on a real 22GB DNA text (index size in RAM 6.5 GB) using 24 threads running on 24 physical cores. Using 1GB pages saves just 3.6% of the time with respect to 2MB pages.

### Configuration

Please refer to more specialized documentation for details on huge pages (for example, the [libhugetlbfs HOWTO](https://github.com/libhugetlbfs/libhugetlbfs/blob/master/HOWTO), the [Huge Pages LWN article](https://lwn.net/Articles/376606), the [Red Hat Performance Tuning Guide](https://access.redhat.com/documentation/en-us/red_hat_enterprise_linux/7/html/performance_tuning_guide/sect-red_hat_enterprise_linux-performance_tuning_guide-memory-configuring-huge-pages)). Here we just present a quick walkthrough. In Linux, one can list the page sizes supported by the system with:
```
$ pagesize --all
4096
2097152
1073741824
```
The default page size used by the kernel is 4KB. To use huge pages, one has to allocate a pool of them. To check which huge page pools are already available, type:
```
$ hugeadm --pool-list
       Size  Minimum  Current  Maximum  Default
   2097152        0        0        0        *       
1073741824        0        0        0        

$ grep Huge /proc/meminfo 
AnonHugePages:         0 kB
HugePages_Total:       0
HugePages_Free:        0
HugePages_Rsvd:        0
HugePages_Surp:        0
Hugepagesize:       2048 kB
```
The `Default` flag means that the currently running system will be using 2MB huge pages. To make it use 1GB huge pages instead, one has to add the following line to the kernel boot command line, and reboot the machine:
```
hugepagesz=1GB default_hugepagesz=1GB
```
To allocate 50 thousand huge pages of 2 MB each at runtime, do the following (as root):
```
$ hugeadm -v --add-temp-swap --pool-pages-min 2MB:50000
hugeadm:INFO: page_size<2MB> adjust<50000> counter<0>
hugeadm:INFO: 1, 50000 -> 50000, 50000
Setting up swapspace version 1, size = 10236 KiB
no label, UUID=1014c55a-d37c-46ad-804e-481861124bf6
hugeadm:INFO: setting HUGEPAGES_TOTAL to 50000
hugeadm:INFO: setting HUGEPAGES_OC to 0
```
Then check if allocation succeeded with:
```
$ hugeadm --pool-list
      Size  Minimum  Current  Maximum  Default
   2097152    50000    50000    50000        *
1073741824        0        0        0         
```
Allocation might fail if the memory is highly fragmented (especially when using 1GB pages). This appears as the following message:
```
hugeadm:WARNING: failed to set pool minimum to 100 became 0
```
If this happens, rebooting the machine will solve the problem. One could also fix a specific configuration of huge pages at boot time, by adding the following line to the kernel boot command line:
```
hugepagesz=1GB default_hugepagesz=1GB hugepages=10
```
The last thing that needs to be configured is which group of users can use huge pages. You can set it by typing (as root):
```
hugeadm --set-shm-group=groupID
```
Now every user in `groupID` can run an already-compiled version of the iterator, so that all its `malloc()` calls are automatically backed by huge pages:
```
$ LD_PRELOAD=libhugetlbfs.so HUGETLB_MORECORE=yes iterator
```
If errors like the following appear, huge pages are not backing the program:
```
libhugetlbfs: WARNING: Hugepage size 2097152 unavailable
libhugetlbfs: WARNING: Hugepage size 1073741824 unavailable
```


