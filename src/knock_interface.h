#ifndef KNOCK_INTERFACE_H
#define KNOCK_INTERFACE_H

#include <vector>
#include <iostream>
#include <Rcpp.h>
#include "dmc_knock.h"
#include "hmm_knock.h"

typedef std::vector< double > vector;
typedef std::vector< std::vector<double> > vector2;
typedef std::vector< vector2 > vector3 ;
typedef std::vector< int > ivector;
typedef std::vector< std::vector<int> > ivector2;

vector  numToVec(const Rcpp::NumericVector &);
vector2 numToVec2(const Rcpp::NumericVector &, int);
vector3 numToVec3(const Rcpp::NumericVector &, int, int);
ivector2 numToIntVec2(const Rcpp::NumericVector &, int);

#endif
