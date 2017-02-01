# kmer-counter

Compilation
-----------

From Altius hosts:

```
$ module add gcc
$ module add glibc
$ make
```

In theory, this should compile okay under OS X, but that is not tested.

Usage
-----

Provide a four-column BED file with the interval's genomic sequence in the fourth column (*i.e.*, ID field).

```
$ ./kmer-counter --k=6 --offset=12195 foo.bed4
```

Notes
-----

I am using a [hash table]<https://en.wikipedia.org/wiki/Hash_table> implementation from [emilib]<https://github.com/emilk/emilib/blob/master/emilib/hash_map.hpp>. A discussion about performance characteristics compared with the C++ STL `std::map` is [available from the author]<http://www.ilikebigbits.com/blog/2016/8/28/designing-a-fast-hash-table>.