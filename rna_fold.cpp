/**
 * Tamara Fleet
 * CSCI 582 Bioinformatics
 * RNA Secondary Structure Prediction
 *
 * Updated:   5/17/2018
 * Purpose:   Implements RNA Secondary Structure Prediction via Base Pair
 *              Maximization using Nussinov's Algorithm
 * Compile:   g++ -Wall -pedantic -g -std=c++11 -o rna_fold rna_fold.cpp
 * Run:       ./rna_fold <read_file>
 **/


#include <iostream>
#include <fstream>
#include <vector>

const int MIN_LOOP_LENGTH = 4; // # of positions required between pairs

int getRead(std::ifstream &file, std::vector<char> *read_vec);
char* nussinovAlgorithm(std::vector<char> *read_vec);
void initMatrix(int** matrix, int size);
int OPT(int i, int j, std::vector<char> *sequence, int** N);
int isComplementary(char a, char b);
void traceback(int i, int j, int** N, char* structure, std::vector<char> *sequence);
void printMatrix(int** matrix, int n);

int main(int argc, char *argv[])
{
  //---------------------------------------------------------------------------
  // Error check number of command line arguments
  if (argc != 2)
  {
    std::cerr <<"Incorrect number of command line arguments. Correct format is:"
      << " ./run <read_file>" << std::endl;
    return 1;
  }

  //---------------------------------------------------------------------------
  // Open read file and get read
  std::ifstream infile(argv[1], std::ios::in);
  if (!infile)
  {
    std::cerr << "Could not open read input file <" << argv[1] << ">." <<
      std::endl;
    return 1;
  }

  std::vector<char> read;
  if (getRead(infile, &read))
  {
    std::cerr << "Error reading in input file <" << argv[1] << ">." <<
      std::endl;
    return 1;
  }

  //---------------------------------------------------------------------------
  // Run Nussinov's Algorithm & print optimum structure in dot-parentheses format
  char* structure = nussinovAlgorithm(&read);
  for (unsigned int i=0; i<read.size(); i++)
  {
    std::cout << structure[i];
  }
  std::cout << std::endl;

  return 0;
}

/**
 * Reads in string from file stream
 *
 * file:      Input file stream associated with the read file
 * read_vec:  Pointer to the vector that holds the RNA sequence
 */
int getRead(std::ifstream &file, std::vector<char> *read_vec)
{
  char value;
  while (file >> value)
  {
    if (file.fail()) {
      return 1;
    }
    //std::cout << value;
    if (value == 'A' || value == 'U' || value == 'G' || value == 'C')
      read_vec->push_back(value);
    else
      return 1;
  }
  return 0;
}


/**
 * Base Pair Maximization using Nussinov's Algorithm
 *
 * read_vec:  Pointer to the vector that holds the RNA sequence
 */
char* nussinovAlgorithm(std::vector<char> *read_vec)
{
  int n = read_vec->size();  // n = sequence length

  //---------------------------------------------------------------------------
  // Allocate and initialize matrix
  int** nussinov = new int*[n];
  for (int i=0; i<n; i++)
    nussinov[i] = new int[n]();
  initMatrix(nussinov, n);

  //---------------------------------------------------------------------------
  // Recursively fill out Nussinov's matrix with OPT scores starting
  //  with i = 0, j = n-1
  nussinov[0][n-1] = OPT(0, n-1, read_vec, nussinov);

  // uncomment the following line to print out final matrix
  //printMatrix(nussinov, n);

  //---------------------------------------------------------------------------
  // Traceback through Nussinov's matrix to determine optimal structure
  char* structure = new char[n]();
  for (int i=0; i<n; i++)
  {
    structure[i] = '.';
  }
  traceback(0, n-1, nussinov, structure, read_vec);
  return structure;
}

/**
 * Initialize matrix with -1s where there can be a pairing
 *
 * matrix:  Nussinov's matrix that holds OPT scores
 * size:    size of sequence/matrix
 */
void initMatrix(int** matrix, int size)
{
  for (int i=0; i<size; i++)
  {
    for (int j=i+1+MIN_LOOP_LENGTH; j<size; j++)
    {
      matrix[i][j] = -1;
    }
  }
}

/**
 * Returns the score of the optimal pairings between indices i and j
 *
 * i:         start index of substring
 * j:         end index of substring
 * sequence:  RNA squence consisting of bases 'A', 'U', 'G', and 'C'
 * N:         Nussinov's matrix that holds optimal scores
 */
int OPT(int i, int j, std::vector<char> *sequence, int** N)
{
  //Base case: i and j are less than MIN_LOOP_LENGTH bases apart; cannot pair
  if (j-i <= MIN_LOOP_LENGTH)
    return 0;
  else if (N[i][j] != -1) // case: if OPT score has already been calculated
  {
    return N[i][j];
  }
  else //i,j are either paired or not paired
  {
    //unpaired
    int unpaired = OPT(i, j-1, sequence, N);

    //paired
    int max_paired = 0;
    for (int t=i; t<j-MIN_LOOP_LENGTH; t++) // all possible pairs
    {
      // check if j can be involved in a pairing with a position t
      if (isComplementary(sequence->at(t), sequence->at(j)))
      {
        int pairing = 1 + OPT(i, t-1, sequence, N) + OPT(t+1, j-1, sequence, N);
        if (pairing > max_paired)
          max_paired = pairing;
      }
    }

    // return max(unpaired, max_paired)
    if (unpaired > max_paired)
    {
      N[i][j] = unpaired; // update OPT matrix with new score
      return unpaired;
    }
    else
    {
      N[i][j] = max_paired; // update OPT matrix with new score
      return max_paired;
    }

  }
}

/**
 * Checks if base a and base b are complementary: AU, UA, CG, or GC
 *
 * a: first base
 * b: second base
 */
int isComplementary(char a, char b)
{
    return ((a == 'A' && b == 'U') || (a == 'U' && b == 'A') || (a == 'C' && b == 'G') || (a == 'G' && b == 'C'));
}


/**
 * Traceback through Nussinov's matrix to determine optimal structure in
 * dot-parentheses format
 *
 * i:         start index of substring
 * j:         end index of substring
 * N:         Nussinov's matrix that holds optimal scores
 * structure: dot-parentheses structure format of the sequence
 * sequence:  RNA squence consisting of bases 'A', 'U', 'G', and 'C'
 */
void traceback(int i, int j, int** N, char* structure, std::vector<char> *sequence)
{
  if (j <= i)
    return;
  else if (N[i][j] == N[i][j-1])
  {
    traceback(i, j-1, N, structure, sequence);
    return;
  }
  else
  {
    for (int k=i; k<j; k++) // find pairing
    {
      if (isComplementary(sequence->at(k), sequence->at(j)))
      {
        int score = 0;
        if (k-1 >= 0) // checks k-1 is in bounds
        {
          score = N[i][k-1];
        }
        if (N[i][j] == 1 + score + N[k+1][j-1]) // k-j pair
        {
          structure[k] = '(';
          structure[j] = ')';
          traceback(i, k-1, N, structure, sequence); //left of k
          traceback(k+1, j-1, N, structure, sequence); // right of k
          return;
        }
      }
    }
  }
}

/**
 * Print int matrix of size n
 *
 * matrix:  2D array to be printed
 * n:       size of matrix
 */
void printMatrix(int** matrix, int n)
{
  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n; j++)
    {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
