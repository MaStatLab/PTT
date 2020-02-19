#ifndef HELPERS_H
#define HELPERS_H

#define MAXVAR 15  // down no more than 14 levels
#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;
using namespace arma;
using namespace std;

#ifndef __USE_MISC
#define __USE_MISC
typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;
#endif

union INDEX_TYPE_t {
  unsigned short var[MAXVAR+1]; //var[MAXVAR] gives the bits that represent the left and right children
  unsigned long long index;
};
typedef INDEX_TYPE_t INDEX_TYPE;

double unifRand();

double log_exp_x_plus_exp_y(double x, double y);

unsigned long long pow2(int k);

void print_index(INDEX_TYPE & I, int level);

void print_index_2(INDEX_TYPE & I, int level, int k);

unsigned int Choose(int n, int k);

unsigned int get_node_index(INDEX_TYPE& I,int level);

INDEX_TYPE init_index(int n,int level);

INDEX_TYPE make_child_index(INDEX_TYPE& I, unsigned short part_dim, int level, ushort which);

INDEX_TYPE get_next_node(INDEX_TYPE& I, int p, int level);

uint convert_to_inverse_base_2(double x, int k);

uint convert_to_base_2(double x, int k);

#endif
