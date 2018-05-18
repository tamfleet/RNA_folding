rna_fold: rna_fold.o
	g++ -Wall -pedantic -g -o rna_fold rna_fold.o

rna_fold.o: rna_fold.cpp
	g++ -Wall -pedantic -g -std=c++11 -c rna_fold.cpp

clean:
	rm -rf rna_fold.o rna_fold
