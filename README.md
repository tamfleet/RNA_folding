# RNA Folding
> Implements RNA Secondary Structure Prediction via Base Pair Maximization using Nussinov's Algorithm

## Background

This method to predict the secondary structure of RNA sequences finds the optimal set of base pairs following the logic that pairing represent bonds and the more bonds there are, the more stable the structure is. However, this method ignores the stabilizing energies for stacked base pairs or destabilizing energies of loops (unpaired) regions.

### Problem

Given a string over {A, T, C, G}, find the optimal set of base pairs i, j such that
* There must be at least 4 positions between i and j (i < j – 4)
* Only A can pair with T and only C can pair with U
* A position is in at most one pairing
* If (i, j) and (i’, j’) are paired, then i < i’ < j < j’ is not permitted (no pseudoknots)


## Compile & Run Commands

Compile using `make` -or-
  `g++ -Wall -pedantic -g -std=c++11 -o rna_fold rna_fold.cpp`

Run using `./rna_fold <read_file>`
`<read-file>`` File containing an RNA sequence consisting of bases 'A', 'U', 'G', and 'C'

## Output

Outputs the optimum structure in dot-parentheses format. This structure can be inputted with the structure into the _forna_ RNA visualization tool found [here](http://rna.tbi.univie.ac.at/forna/). Without the dot-parentheses, _forna_ predicts the structure using a minimized free energy (MFE) method.

## Source
