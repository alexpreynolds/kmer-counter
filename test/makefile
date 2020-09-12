SHELL           := /bin/bash
PWD             := $(shell pwd)
BIN              = ../kmer-counter

all: 2mer

2mer:
	cd .. && $(MAKE) clean && $(MAKE) && cd $(PWD)
	$(BIN) --bed --k=2 --offset=100 --results-dir="2mer" 2mer.bed
	diff -s 2mer/count.bed 2mer-expected/count.bed
	diff -s 2mer/map.txt 2mer-expected/map.txt

2mer_valgrind:
	valgrind --tool=memcheck --track-origins=yes --track-fds=yes --leak-check=yes $(BIN) --k=2 --offset=100 --results-dir="2mer" 2mer.bed

4mer_fasta:
	cd .. && $(MAKE) clean && $(MAKE) && cd $(PWD)
	./generate-random-sequences.py 5 1000 123 | fold -w 80 -s > 4mer-test.fa
	$(BIN) --fasta --k=4 4mer-test.fa > 4mer-observed.txt
	diff -s 4mer-observed.txt 4mer-expected.txt

clean:
	rm -rf 2mer
	rm -rf *~
	rm -f 4mer-observed.txt
	rm -f 4mer-test.fa
	cd .. && $(MAKE) clean && cd $(PWD)
