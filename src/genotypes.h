#ifndef GENOTYPES_H
#define GENOTYPES_H

/*
Knockoffs for phased haplotypes
*/

#include <vector>
#include <random>
#include "utils.h"

typedef std::vector< std::vector<double> > matrix;
typedef std::vector< std::vector<int> > imatrix;

namespace knockoffs {
  class GenotypeModel {
  public:
    GenotypeModel(const std::vector<double> & r, const matrix & alpha, const matrix & theta, int seed);
    imatrix sample(const imatrix & X);
    std::vector<int> sample(const std::vector<int> & X);
  private:
    void sampleViterbi(const std::vector<int> & X);
    int pair_to_index(int i, int j);
    std::vector<int> single_to_indices(int j);
    double emission_prob(int j, int x, int k, int l);
    double emission_prob(int j, int x, int m);
    imatrix table;
    void knockoffMC(const std::vector<int> & H);
    void emission(const std::vector<int> & Hk);
    matrix theta, a;
    std::vector<double> b;
    matrix beta;
    double beta_const;
    std::vector<double> weights;
    int p, K, nStates;
    std::vector<int> H, Hk, Xk;
    // Partition function for Markov chain knockoffs
    std::vector<double> Z, Z_old, C, C_sums_partial;
    std::vector<int> indices;
    // Random number generation
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    std::mt19937 gen2;
  };
}

#endif
