#include "stdafx.h"
#include "modeltree_node.h"
#include "gbt.h"

unsigned long long pow2(int k) {

  unsigned long long res = (unsigned long long) 1 << k;
  return res;
}

void print_index(INDEX_TYPE & I, int level) {
  ushort x_curr = 0;
  ushort index_prev_var = 0;

  for (int i = 0; i < level; i++) {
    x_curr += I.var[i] - index_prev_var - 1;
    cout << "X" << x_curr << "=" << ((I.var[MAXVAR] >> i) & (ushort) 1) << ",";
    index_prev_var = I.var[i];
  }
  cout << I.var[MAXVAR];
}

void print_index_2(INDEX_TYPE & I, int level, int k) {
  ushort x_curr = 0;
  ushort index_prev_var = 0;
  ushort lower = 0;
  ushort x_curr_count = -1;

  for (int i = 0; i < level; i++) {
    if ( I.var[i] - index_prev_var - 1 > 0 ) { // next variable
      if (x_curr > 0) cout << ";";
      cout << "X" << x_curr << "_l=" << lower << ", X" << x_curr << "_u=" << lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;
 
      lower = 0;
      x_curr_count = 0;
    }  
    else { 
      x_curr_count++;
    }
      
    x_curr += I.var[i] - index_prev_var - 1;
    lower |= (((I.var[MAXVAR] >> i) & (ushort) 1)) << (k-x_curr_count-1);  
    index_prev_var = I.var[i];
  }
  
  if (level > 0) {
    if (x_curr >0) cout << ";";
    cout << "X" << x_curr << "_l=" << lower << ", X" << x_curr << "_u=" << lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;
  }

  if (level == 0) {

  }
}

unsigned int Choose(int n, int k) {
    unsigned int c = 1;
    unsigned int d = 1;
    for (int i = 0; i <k; i++) {
        c *= (n-i); d *= (i+1);
    }
    return c  / d;
}

unsigned int get_node_index(INDEX_TYPE& I,int level){
     unsigned long long  r = 0;
 
     unsigned long long numerator = 1;
     unsigned long long denominator = 1;
     
     for (int i = 0; i < level; i++) {
       numerator = 1;
       denominator *= (i+1);     
       for (int j = 1; j <= i+1; j++) {
     	 numerator *= I.var[i]-j;
      }
       r += numerator / denominator;
     }

     return NUMNODEVAR*(r*pow2(level) + (unsigned long long) I.var[MAXVAR]);
}
inline INDEX_TYPE init_index(int n,int level) {
    INDEX_TYPE init;

    for (int i=0; i < level; i++) {init.var[i] = i+1;}    
    for (int i=level; i <= MAXVAR; i++) {init.var[i] = 0;}
    
    return init;
}

inline INDEX_TYPE make_child_index(INDEX_TYPE& I, unsigned short part_dim, int level, ushort which) {
    INDEX_TYPE child_index = I;
    unsigned short data = part_dim+1; 
    int i;
    int j;
    int x_curr;
    int child_index_var_prev;

    if (level == 0) {
      x_curr=1;
      i=0;
      child_index_var_prev = 0;
    }
    
    else {
       x_curr = child_index.var[0]; // current dimension
       child_index_var_prev = child_index.var[0];
       i = 1;
    }

    while (i<MAXVAR) {
        while (child_index.var[i] >0 && data >= x_curr ) {
	  x_curr += child_index.var[i] - child_index_var_prev - 1;
	  child_index_var_prev = child_index.var[i];
	  i++;

        }

        if (child_index.var[i] == 0 && data >= x_curr) {
	    child_index.var[i] = data - x_curr + 1 + child_index_var_prev;
	    j=i;
	    i=MAXVAR;
        } else { // this corresponds to the first i such that data < x_curr
	    for (int h = level; h >= i; h--) {
	      child_index.var[h] = child_index.var[h-1]+1;
	    }
	    
	    child_index.var[i-1] = child_index.var[i] - (x_curr - data + 1); 
	    
	    j=i-1; 
	    i=MAXVAR;
        }
    }
    child_index.var[MAXVAR] = ((I.var[MAXVAR] << 1) & ((~((ushort)0)) << (j+1) )) | ((ushort) which << j) | ( I.var[MAXVAR] & ~((~((ushort)0)) << j) ); // update the bits of the child

    return child_index;
}


inline INDEX_TYPE get_next_node(INDEX_TYPE& I, int p, int level) {
  
    INDEX_TYPE node = I;
    int i = level-1; int j = p + level -2;
    while (i>=0 && node.var[i] == j+1) {i--;j--;}
    if (i < 0) { //reach the end of nodes
      for (int h = 0; h <= MAXVAR; h++) {
	node.var[h]=0; //invalid node
      }
    } else {
        node.var[i] += 1;
        for (j=i+1;j<level;j++) {
            node.var[j] = node.var[i]+ j-i;
        }
    }
    node.var[MAXVAR] = 0;
    
    return node;
}

inline uint convert_to_inverse_base_2(double x, int k) {
  
  uint x_base_2 = (uint) floor (x * pow2(k) );
  uint x_base_inverse_2 = 0;

  for (int i = 0; i < k; i++) {
    
    x_base_inverse_2 |= ((x_base_2 >> (k - 1 -i)) & 1) << i;
  }

  return x_base_inverse_2;
}

inline uint convert_to_base_2(double x, int k) {
  
  uint x_base_2 = (uint) floor (x * pow2(k) );

  return x_base_2;
}


