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

#ifndef GENOTYPES_CPP
#define GENOTYPES_CPP

#include "genotypes.h"

using namespace knockoffs;

GenotypeModel::GenotypeModel(const std::vector<double> & r, const matrix & alpha, const matrix & _theta,
                             int seed){
  // Store input parameters
  theta = _theta;
  p = alpha.size();
  K = alpha[0].size();
  nStates = (K*(K+1))/2;

  // Process MC parameters
  a = matrix(p, std::vector<double> (K));
  b = std::vector<double>(p, 0);
  b[0] = 0;
  for(int l=0; l<K; l++) {
    a[0][l] = alpha[0][l];
  }
  for(int j=1; j<p; j++) {
    b[j] = std::exp(-r[j]);
    for(int l=0; l<K; l++) {
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
  table = imatrix(nStates, std::vector<int>(2));
  for(int i=1; i<K; i++) {
    for(int j=0; j<=i; j++) {
      int m = pair_to_index(i,j);
      table[m][0] = i;
      table[m][1] = j;
    }
  }
  C = std::vector<double>(nStates,0);

  Hk = std::vector<int> (p);
  Xk = std::vector<int> (p);
  // Knockoffs for Markov chains
  Z = std::vector<double> (nStates);
  Z_old = std::vector<double> (nStates,1);
  indices = std::vector<int>(K);
  C_sums_partial = std::vector<double>(K);
}

int GenotypeModel::pair_to_index(int k, int l) {
  // Make sure that k >= l
  if(k<l) {
    int tmp = l;
    l = k;
    k = tmp;
  }
  return((k*(k+1))/2+l);
}

std::vector<int> GenotypeModel::single_to_indices(int j) {
  std::vector<int> indices(K);
  for(int i=0; i<K; i++) {
    indices[i] = pair_to_index(i,j);
  }
  return(indices);
}

double GenotypeModel::emission_prob(int j, int x, int m) {
  return emission_prob(j, x, table[m][0], table[m][1]);
}

double GenotypeModel::emission_prob(int j, int x, int k, int l) {
  double prob = 0;
  if(x==0) {
    prob = (1-theta[j][k])*(1-theta[j][l]);
  }
  else if(x==1) {
    prob = theta[j][k]*(1-theta[j][l]) + theta[j][l]*(1-theta[j][k]);
  }
  else if(x==2) {
    prob = theta[j][k]*theta[j][l];
  }
  return(prob);
}

void GenotypeModel::sampleViterbi(const std::vector<int> & X) {
  // Compute backward weights
  std::fill(beta[p-1].begin(), beta[p-1].end(), 1.0);

  for(int j=p-2; j>=0; j--) {
    // Precompute scalar sum
    double beta_pre_scalar = 0;
    for(int k=0; k<K; k++) {
      for(int l=0; l<K; l++) {
        int m = pair_to_index(k,l);
        beta_pre_scalar += a[j+1][k] * a[j+1][l] * emission_prob(j+1, X[j+1], m) * beta[j+1][m];
      }
    }
    // Precompute vector sum
    std::vector<double> beta_pre_vec(K,0);
    for(int k=0; k<K; k++) {
      std::vector<int> indices = single_to_indices(k);
      for(int l=0; l<K; l++) {
        int k1 = table[indices[l]][0];
        int k2 = table[indices[l]][1];
        double factor = 0;
        if(k1==k2) {
          factor = a[j+1][k1];
        }
        else {
          factor = a[j+1][k1] * (k2==k) + a[j+1][k2] * (k1==k);
        }
        beta_pre_vec[k] += factor * emission_prob(j+1, X[j+1], indices[l]) * beta[j+1][indices[l]];
      }
    }
    // Assemble backward weights
    double betaSum = 0.0;
    for(int k=0; k<K; k++) {
      for(int l=0; l<=k; l++) {
        int m = pair_to_index(k,l);
        double pEmit = emission_prob(j+1, X[j+1], m);
        beta[j][m]  = beta_pre_scalar;
        beta[j][m] += b[j+1] * (beta_pre_vec[k] + beta_pre_vec[l]);
        beta[j][m] += b[j+1] * b[j+1] * pEmit * beta[j+1][m];
        betaSum += beta[j][m];
      }
    }

    // Normalize the backward weights
    for(int m=0; m<nStates; m++) {
      beta[j][m] /= betaSum;
    }

  }

  // Forward sampling
  double weights_sum = 0.0;
  for(int k=0; k<K; k++) {
    for(int l=0; l<=k; l++) {
      int m = pair_to_index(k,l);
      if(k==l) {
        weights[m] = a[0][k] * a[0][l];
      }
      else {
        weights[m] = 2.0 * a[0][k] * a[0][l];
      }
      weights[m] *= emission_prob(0, X[0], m) * beta[0][m];
      weights_sum += weights[m];
    }
  }

  for(int k=0; k<nStates; k++) {
    weights[k] /= weights_sum;
  }
  H[0] = weighted_choice(dis(gen),weights);

  for(int j=1; j<p; j++) {
    weights_sum = 0.0;
    for(int k=0; k<K; k++) {
      for(int l=0; l<=k; l++) {
        int m = pair_to_index(k,l);     // New state
        int k0  = table[H[j-1]][0];
        int l0  = table[H[j-1]][1];
        weights[m]  = (a[j][k]+b[j] * (double)(k==k0)) * (a[j][l]+b[j] * (double)(l==l0));
        if(k!=l) {
          weights[m] += (a[j][k]+b[j] * (double)(k==l0)) * (a[j][l]+b[j] * (double)(l==k0));
        }
        weights[m] *= emission_prob(j, X[j], m) * beta[j][m];
        weights_sum += weights[m];
      }
    }

    for(int k=0; k<nStates; k++) {
      weights[k] /= weights_sum;
    }
    H[j] = weighted_choice(dis(gen),weights);
  }
}

void GenotypeModel::knockoffMC(const std::vector<int> & H) {
  std::fill(Z_old.begin(), Z_old.end(), 1.0);
  double weights_sum;
  for(int j=0; j<p; j++) {
    std::fill(weights.begin(), weights.end(), 1.0);
    std::fill(Z.begin(), Z.end(), 0.0);
    std::fill(C.begin(), C.end(), 1.0);
    std::fill(C_sums_partial.begin(), C_sums_partial.end(), 0.0);
    weights_sum = 0;

    if(j<p-1) {
      // Precompute vector C
      double C_sum = 0;
      for(int k1=0; k1<K; k1++) {
        for(int k2=0; k2<=k1; k2++) {
          int k12 = pair_to_index(k1,k2);
          if(j==0) {
            C[k12] = a[j][k1] * a[j][k2];
            if(k1!=k2) {
              C[k12] += a[j][k1] * a[j][k2];
            }
          }
          else {
            int l1  = table[H[j-1]][0];
            int l2  = table[H[j-1]][1];
            int l1k = table[Hk[j-1]][0];
            int l2k = table[Hk[j-1]][1];
            if(k1==k2) {
              C[k12]  = (a[j][k1]+b[j]*(k1==l1)) * (a[j][k2]+b[j]*(k2==l2));
              C[k12] *= (a[j][k1]+b[j]*(k1==l1k)) * (a[j][k2]+b[j]*(k2==l2k));
            }
            else {
              double tmp1 = (a[j][k1]+b[j]*(k1==l1)) * (a[j][k2]+b[j]*(k2==l2));
              tmp1 += (a[j][k2]+b[j]*(k2==l1)) * (a[j][k1]+b[j]*(k1==l2));
              double tmp2 = (a[j][k1]+b[j]*(k1==l1k)) * (a[j][k2]+b[j]*(k2==l2k));
              tmp2 += (a[j][k2]+b[j]*(k2==l1k)) * (a[j][k1]+b[j]*(k1==l2k));
              C[k12] = tmp1 * tmp2;
            }
          }
          C[k12] /= Z_old[k12];
          // Precompute sum of C
          C_sum += C[k12];
          // Precompute partial sums of C
          C_sums_partial[k1] += C[k12];
          C_sums_partial[k2] += C[k12];
        }
      }

      // Compute partition function
      for(int k1=0; k1<K; k1++) {
        for(int k2=0; k2<=k1; k2++) {
          int k12 = pair_to_index(k1,k2);
          Z[k12]  = a[j+1][k1] * a[j+1][k2] * C_sum;
          Z[k12] += b[j+1] * b[j+1] * C[k12];
          Z[k12] += a[j+1][k1] * b[j+1] * C_sums_partial[k2];
          if(k1!=k2) {
            Z[k12] += a[j+1][k1] * a[j+1][k2] * C_sum;
            Z[k12] += a[j+1][k2] * b[j+1] * C_sums_partial[k1];
          }
        }
      }
    }

    // Compute weights
    for(int k1=0; k1<K; k1++) {
      for(int k2=0; k2<=k1; k2++) {
        int k12 = pair_to_index(k1,k2);
        if(j>0) {
          int l1  = table[H[j-1]][0];
          int l2  = table[H[j-1]][1];
          int l1k = table[Hk[j-1]][0];
          int l2k = table[Hk[j-1]][1];
          double weight_left = 1.0;
          if(k1==k2) {
            weight_left = (a[j][k1]+b[j]*(k1==l1)) * (a[j][k2]+b[j]*(k2==l2));
            weight_left *= (a[j][k1]+b[j]*(k1==l1k)) * (a[j][k2]+b[j]*(k2==l2k));
          }
          else {
            double tmp1 = (a[j][k1]+b[j]*(k1==l1)) * (a[j][k2]+b[j]*(k2==l2));
            tmp1 += (a[j][k2]+b[j]*(k2==l1)) * (a[j][k1]+b[j]*(k1==l2));
            double tmp2 = (a[j][k1]+b[j]*(k1==l1k)) * (a[j][k2]+b[j]*(k2==l2k));
            tmp2 += (a[j][k2]+b[j]*(k2==l1k)) * (a[j][k1]+b[j]*(k1==l2k));
            weight_left = tmp1 * tmp2;
          }          
          weights[k12] *= weight_left;
        }
        if(j<p-1) {
          int l1  = table[H[j+1]][0];
          int l2  = table[H[j+1]][1];
          double weight_right = (a[j+1][l1]+b[j+1]*(k1==l1)) * (a[j+1][l2]+b[j+1]*(k2==l2));
          if(l1!=l2) {
            weight_right += (a[j+1][l1]+b[j+1]*(k2==l1)) * (a[j+1][l2]+b[j+1]*(k1==l2));
          }
          weights[k12] *= weight_right;
        }
        weights[k12] /= Z_old[k12];
        Z_old[k12] = Z[k12];
        weights_sum += weights[k12];
      }
    }

    // Normalize weights
    for(int k12=0; k12<nStates; k12++) {
      weights[k12] /= weights_sum;
    }

    // Sample knockoff copy at this site
    Hk[j] = weighted_choice(dis(gen2),weights);
  }
}

void GenotypeModel::emission(const std::vector<int> & Hk) {
  std::vector<double> weights3(3,1.0);
  for(int j=0; j<p; j++) {
    weights3[0] = emission_prob(j, 0, Hk[j]);
    weights3[1] = emission_prob(j, 1, Hk[j]);
    weights3[2] = emission_prob(j, 2, Hk[j]);
    Xk[j] = weighted_choice(dis(gen),weights3);
  }
}

std::vector<int> GenotypeModel::sample(const std::vector<int> & X) {
  sampleViterbi(X);
  knockoffMC(H);
  emission(Hk);
  return(Xk);
}

imatrix GenotypeModel::sample(const imatrix & X) {
  unsigned int n = X.size();
  imatrix XkMatrix(n, std::vector<int>(p));
  for(unsigned int i=0; i<X.size(); i++) {
    XkMatrix[i] = sample(X[i]);
  }
  return(XkMatrix);

}

#endif
