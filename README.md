# kmer-counter

Compilation
-----------

Run `make` to build the `kmer-counter` binary.

This has been compiled under Ubuntu 18.04.4, Cygwin 3.1.4, and Mac OS X 10.15.3, using concurrent GCC/glibc and Clang toolkits.

Usage
-----

Run `kmer-counter --help` for a list of options.

There are a couple ways to use this.

1. You can provide a single-line FASTA input and write counts to standard output, *e.g.*:

```
$ ./kmer-counter --fasta --k=6 sequences.fa
>foo    CGTTAA:1 TTAACG:1
>bar    TTCTTA:1 TAGGGC:1 AAATTC:1 GTGGAA:1 AACTTC:1 ...
...
```

2. For a more complex use case, you can provide a four-column BED file with the interval's genomic sequence in the fourth column (*i.e.*, ID field), along with the number *k* for the k-mers you want to count, an *offset* value for mer-keys (explained below), and a *results directory* to write results, *e.g.*:

```
$ ./kmer-counter --bed --k=6 --offset=12195 --results-dir="6mers" intervals.bed4
```

The above example generates 6-mers of the sequences from the file `intervals.bed4`.

The results are stored in a folder called `6mers`, which contains two files `count.bed` and `map.txt`.

The first file `count.bed` contains a BED file of intervals from `intervals.bed4`, where the fourth column contains a space-delimited pair of "mer"-keys and the number of times that key is seen. Mer-keys are numbers which begin at the `offset` value provided on the command-line.

The second file `map.txt` contains a tab-delimited pairing of mers and their mer-key, as found in `count.bed`.

Notes
-----

I am using a [hash table](https://en.wikipedia.org/wiki/Hash_table) implementation from [Emil Ernerfeldt](https://github.com/emilk/emilib/blob/master/emilib/hash_map.hpp). A discussion about performance characteristics compared with the C++ STL `std::unordered_map` is [available from the author](http://www.ilikebigbits.com/blog/2016/8/28/designing-a-fast-hash-table).
