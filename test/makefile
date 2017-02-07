SHELL           := /bin/bash
PWD             := $(shell pwd)
BIN              = ../kmer-counter

all: 2mer

2mer:
	cd .. && $(MAKE) clean && $(MAKE) && cd $(PWD)
	$(BIN) --k=2 --offset=100 --results-dir="2mer" 2mer.bed
	diff -s 2mer/count.bed 2mer-expected/count.bed
	diff -s 2mer/map.txt 2mer-expected/map.txt

clean:
	rm -rf 2mer
	rm -rf *~
