/*
  This file is part of SNPknock.

    Copyright (C) 2017-2018 Matteo Sesia

    SNPknock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SNPknock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SNPknock.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef HAPLOTYPES_CPP
#define HAPLOTYPES_CPP

#include "haplotypes.h"

using namespace knockoffs;

HaplotypeModel::HaplotypeModel(const std::vector<double> & r, const matrix & alpha, const matrix & _theta,
                               int seed){
  // Store input parameters
  theta = _theta;
  p = alpha.size();
  nStates = alpha[0].size();
  // Process MC parameters
  a = matrix(p, std::vector<double> (nStates));
  b = std::vector<double>(p, 0);
  b[0] = 0;
  for(int l=0; l<nStates; l++) {
    a[0][l] = alpha[0][l];
  }
  for(int j=1; j<p; j++) {
    b[j] = std::exp(-r[j]);
    for(int l=0; l<nStates; l++) {
      a[j][l] = (1.0-b[j])*alpha[j][l];
    }
  }
  // Initialize random number generator
  gen = std::mt19937();
  gen2 = std::mt19937();
  dis = std::uniform_real_distribution<double>(0.0,1.0);
  gen.seed(seed);
  gen2.seed(seed+100000);

  // Initialize variables
  weights = std::vector<double> (nStates);
  beta = matrix(p, std::vector<double> (nStates));
  H = std::vector<int> (p);
  Hk = std::vector<int> (p);
  Xk = std::vector<int> (p);
  // Knockoffs for Markov chains
  Z = std::vector<double> (nStates);
  Z_old = std::vector<double> (nStates);
}

void HaplotypeModel::sampleViterbi(const std::vector<int> & X) {
  // Compute backward weights
  std::fill(beta[p-1].begin(), beta[p-1].end(), 1.0);
  for(int j=p-2; j>=0; j--) {
    beta_const = 0.0;
    for(int l=0; l<nStates; l++) {
      if ( X[j+1] == 1) {
        beta_const += a[j+1][l] * theta[j+1][l] * beta[j+1][l];
      }
      else {
        beta_const += a[j+1][l] * (1.0-theta[j+1][l]) * beta[j+1][l];
      }
    }
    double betaSum = 0.0;
    for(int k=0; k<nStates; k++) {
      if ( X[j+1] == 1) {
        beta[j][k] = beta_const + b[j+1] * theta[j+1][k] * beta[j+1][k];
      }
      else {
        beta[j][k] = beta_const + b[j+1] * (1.0-theta[j+1][k]) * beta[j+1][k];
      }

      betaSum += beta[j][k];
    }
    for(int k=0; k<nStates; k++) {
      beta[j][k] /= betaSum;
    }
  }

  // Forward sampling
  double weights_sum = 0.0;
  for(int k=0; k<nStates; k++) {
    if ( X[0] == 1) {
      weights[k] = a[0][k] * theta[0][k] * beta[0][k];
    }
    else {
      weights[k] = a[0][k] * (1.0-theta[0][k]) * beta[0][k];
    }
    weights_sum += weights[k];
  }
  for(int k=0; k<nStates; k++) {
    weights[k] /= weights_sum;
  }
  H[0] = weighted_choice(dis(gen),weights);

  for(int j=1; j<p; j++) {
    weights_sum = 0.0;
    for(int k=0; k<nStates; k++) {
      if ( X[j] == 1) {
        weights[k] = (a[j][k] + b[j]*(double)(k==H[j-1])) * theta[j][k] * beta[j][k];
      }
      else {
        weights[k] = (a[j][k] + b[j]*(double)(k==H[j-1])) * (1.0-theta[j][k]) * beta[j][k];
      }
      weights_sum += weights[k];
    }
    for(int k=0; k<nStates; k++) {
      weights[k] /= weights_sum;
    }
    H[j] = weighted_choice(dis(gen),weights);
  }
}

void HaplotypeModel::knockoffMC(const std::vector<int> & H) {
  std::fill(Z_old.begin(), Z_old.end(), 1.0);
  double weights_sum, Z_sum;

  for(int j=0; j<p; j++) {
    std::fill(weights.begin(), weights.end(), 1.0);
    std::fill(Z.begin(), Z.end(), 0.0);
    weights_sum = 0;
    Z_sum = 0;

    // Precompute sum for partition function
    for(int k=0; k<nStates; k++) {
      if( j==0 ) {
        Z_sum += a[j][k];
      }
      else {
        Z_sum += (a[j][k] + (double)(k==H[j-1])) * (a[j][k] + (double)(k==Hk[j-1])) / Z_old[k];
      }
    }

    // Compute partition function
    for(int k=0; k<nStates; k++) {
      if(j<p-1) {
        Z[k] += a[j+1][k] * Z_sum;
        if(j==0) {
          Z[k] += b[j+1] * a[j][k];
        }
        else {
          Z[k] += b[j+1] * (a[j][k] + (double)(k==H[j-1])) * (a[j][k] + (double)(k==Hk[j-1])) / Z_old[k];
        }
      }
    }

    // Compute sampling weights
    for(int k=0; k<nStates; k++) {
      if(j < p-1) {
        weights[k] *= (a[j+1][H[j+1]]+b[j+1]*(double)(k==H[j+1]));
      }
      if(j==0) {
        weights[k] *= a[j][k];
      }
      else {
        weights[k] *= (a[j][k]+b[j]*(double)(k==H[j-1]));
        weights[k] *= (a[j][k]+b[j]*(double)(k==Hk[j-1]));
      }
      weights[k] /= Z_old[k];
      Z_old[k] = Z[k];
      weights_sum += weights[k];
    }

    // Normalize weights
    for(int k=0; k<nStates; k++) {
      weights[k] /= weights_sum;
    }

    Hk[j] = weighted_choice(dis(gen2),weights);
  }
}

void HaplotypeModel::emission(const std::vector<int> & Hk) {
  std::vector<double> weights2(2,1.0);
  // A little slower, but random seeds are compatible with older code
  for(int j=0; j<p; j++) {
    weights2[0] = 1.0-theta[j][Hk[j]];
    weights2[1] = theta[j][Hk[j]];
    Xk[j] = weighted_choice(dis(gen),weights2);
  }
}

std::vector<int> HaplotypeModel::sample(const std::vector<int> & X) {
  sampleViterbi(X);
  knockoffMC(H);
  emission(Hk);
  return(Xk);
}

imatrix HaplotypeModel::sample(const std::vector<std::vector<int> > & X) {
  unsigned int n = X.size();
  imatrix XkMatrix(n, std::vector<int>(p));
  for(unsigned int i=0; i<X.size(); i++) {
    XkMatrix[i] = sample(X[i]);
  }
  return(XkMatrix);

}

#endif
