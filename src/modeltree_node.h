#ifndef MODELTREE_NODE_H
#define MODELTREE_NODE_H
#define MAXVAR 12

union INDEX_TYPE_t {
  unsigned short var[MAXVAR+1]; //var[MAXVAR] gives the bits that represent the left and right children
  unsigned long long index;
};
typedef INDEX_TYPE_t INDEX_TYPE;

#endif
