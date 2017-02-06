SHELL           := /bin/bash
PWD             := $(shell pwd)
BIN              = ../kmer-counter

all: 2mer

2mer:
	cd .. && $(MAKE) clean && $(MAKE) && cd $(PWD)
	$(BIN) --k=2 --offset=100 --results-dir="2mer" 2mer.bed

clean:
	rm -rf 2mer
	rm -rf *~
