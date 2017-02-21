#include "cgbt.h"

#include <vector>
using namespace std;


CondGBT::CondGBT(Mat< unsigned int > X, Mat< unsigned int> Y, int k_X, int k_Y, int p_X, int p_Y,
                 double rho0_X, int rho0_mode_X, double rho0_Y, int rho0_mode_Y, int tran_mode_Y,
                 double lognu_lowerbound_Y, double lognu_upperbound_Y, int n_grid_Y = 1, int n_s_Y=4, double beta_Y=0.1):
              GBT(X,k_X,p_X,rho0_X,rho0_mode_X,tran_mode_Y,lognu_lowerbound_Y,lognu_upperbound_Y,n_grid_Y,1,beta_Y),
              p_Y(p_Y),k_Y(k_Y),rho0_Y(rho0_Y), rho0_mode_Y(rho0_mode_Y),n_s_Y(n_s_Y) {

         nobs = X.n_rows;
         init(X,Y);
}

CondGBT::~CondGBT() {
  clear();
}

double CondGBT::get_root_logphi() {

  return models[0][2+n_s];
  // return GBT::get_root_logphi() - models[0][0] * k * log(2.0);

}

double CondGBT::get_root_logrho() {

  // return log(make_rho0(0)) + models[0][2] - get_root_logphi();
  return GBT::get_root_logrho();
}

void CondGBT::init(Mat< unsigned int > X, Mat< unsigned int > Y) {
  unsigned long long i, j;
  unsigned long long l;

  models = new double*[k+1];
  modelscount = new unsigned int[k+1];
  gbt_ptrs = new GBT **[k+1];


  for (i=0; i <= k; i++) {
    modelscount[i] = Choose(p + i - 1, i);
    models[i] = new double[(unsigned long long) (modelscount[i]*numnodevar) << i];
    gbt_ptrs[i] = new GBT *[(unsigned long long) modelscount[i] << i];

    for (j=0; j < modelscount[i] ; j++) {
      for (l=0; l < pow2(i); l++) {
        models[i][(j*pow2(i)+l)*numnodevar] = 0; // Initialize NODE_CURR[0] to 0
        gbt_ptrs[i][j*pow2(i)+l] = NULL;
      }
    }
  }

  // step 1: add data into the models
  INDEX_TYPE I_root = init_index(p,0);

  for (i=0; i < nobs; i++) { //since we do it only once, do this for now
    add_data_to_subtree(I_root,0,1,0,X.row(i).t(),Y.row(i).t()); // code the data matrix such that each observation is in contiguity
  }

  // step 2: initialize prior parameters
  init_prior();

}

void CondGBT::clear() {

  for (int i =0; i <= k; i++) {

    for (unsigned long long j = 0; j < ((unsigned long long) modelscount[i] << i); j++) {
      if (gbt_ptrs[i][j] != NULL) {
         gbt_ptrs[i][j]->clear();
      }
    }

    delete [] gbt_ptrs[i];
  }


/*
  delete [] models; models = NULL;
  delete [] modelscount; modelscount= NULL;

  for (int s = 0; s < n_s; s++) {
    delete [] logrho_mat[s];
  }

  delete [] logrho_mat;
  delete [] logrho_vec;
  logrho_mat=NULL;
*/
}

GBT ** CondGBT::get_node_gbt_ptr(INDEX_TYPE & I, int level) {
  return &gbt_ptrs[level][get_node_index(I,level)/numnodevar];
}

void CondGBT::sample_subtree(INDEX_TYPE& I,int level,int s) {

  int i,t;
  double u, v, cum_prob_curr;

  double *CHILD_0;
  double *CHILD_1;

  make_posterior_logrho_vec(logrho_vec,I,level,s);

  double log_prob = 0;

  // first sample state
  u = unifRand();

  if (level == k || log(u) <= logrho_vec[n_s]) { // pruned?

    t = n_s;
    sample_nodes.push_back(make_pair(make_pair(I,level),make_pair(1.0/0,0)));

  }

  else {
    t = 0; // not stopped

    cum_prob_curr = exp(logrho_vec[n_s]);

    // next sample partition
    v = unifRand();
    cum_prob_curr = 0;

    for(i = 0; i < p && cum_prob_curr < v; i++  ) {
      cum_prob_curr += get_add_prob(I,i,0,level);
    }

    // i - 1 is the dimension to divide

    INDEX_TYPE child_index_0 = make_child_index(I, i-1, level, 0);
    INDEX_TYPE child_index_1 = make_child_index(I, i-1, level, 1);

    CHILD_0 = get_node(child_index_0,level+1);
    CHILD_1 = get_node(child_index_1,level+1);

    sample_nodes.push_back(make_pair(make_pair(I,level),make_pair(0,0)));
    sample_subtree(child_index_0,level+1,0);
    sample_subtree(child_index_1,level+1,0);
  }
}


void CondGBT::find_hmap(int print) {
  INDEX_TYPE I_root = init_index(p,0);
  hmap_nodes.clear();
  find_hmap_subtree(I_root,0);   // the sampled nodes are stored in sample_nodes
  if (print) { // print the sampled nodes
    vector< pair<INDEX_TYPE, pair<int,int> > >::iterator it;
    for (it = hmap_nodes.begin(); it < hmap_nodes.end(); it++) {
      print_index(it->first, it->second.first);
      cout << endl;
      print_index_2(it->first, it->second.first, k);
      cout << endl;
    }
  }
}

void CondGBT::find_hmap_subtree(INDEX_TYPE& I, int level) {

  int j,t;
  double *node = get_node(I,level);
  double logrho0 = log(make_rho0(level));

  int curr_max_state, curr_max_dim;
  double curr_phi_max, curr_Zj_max, curr_Zj;

  if (level == k) { // stopping?
    hmap_nodes.push_back(make_pair(I,make_pair(level,n_s+1)));
  }

  else {
    curr_max_state = n_s;
    curr_phi_max = node[2+n_s+n_s];

    for(t = 0; t < n_s; t++) {
      if (node[2+n_s+t] > curr_phi_max) {
        curr_max_state = t;
        curr_phi_max = node[2+n_s+t];
      }
    }

    // so now the hmap state is curr_max_state
    hmap_nodes.push_back(make_pair(I, make_pair(level,curr_max_state+1)));

    if (curr_max_state < n_s) {

      curr_max_dim = 0;
      curr_Zj_max = get_add_prob(I,0,curr_max_state,level);

      for (j = 1; j < p; j++) {
        curr_Zj = get_add_prob(I,j,curr_max_state,level);

        if (curr_Zj > curr_Zj_max) {
          curr_max_dim = j;
          curr_Zj_max = curr_Zj;
        }
      }


      INDEX_TYPE child_index_0 = make_child_index(I, curr_max_dim, level, 0);
      INDEX_TYPE child_index_1 = make_child_index(I, curr_max_dim, level, 1);

      find_hmap_subtree(child_index_0,level+1);
      find_hmap_subtree(child_index_1,level+1);
    }
  }

}

pair< vector< vector< ushort > >, vector< double > > CondGBT::find_hmap_part() {
  vector< vector< ushort > > part_points;
  vector< double > hmap_rhos;

  vector< pair<INDEX_TYPE, pair<int,int> > >::iterator it;
  int i;
  vector< ushort > v(2*p+2);

  INDEX_TYPE I;
  int level;
  int state;

  ushort x_curr = 0;
  ushort index_prev_var = 0;
  ushort lower = 0;
  ushort x_curr_count = -1;

  find_hmap(0);

  for (it = hmap_nodes.begin(); it < hmap_nodes.end(); it++) {

    I = it->first;
    level = it->second.first;
    state = it->second.second;
    v[2*p] = level;
    v[2*p+1] = state;

    for (i = 0; i < p; i++) {
      v[2*i] = 0;
      v[2*i+1] = ((ushort) 1 << k) - 1;
    }

    x_curr = 0;
    index_prev_var = 0;
    lower = 0;
    x_curr_count = -1;

    for (i = 0; i < level; i++) {
      if ( I.var[i] - index_prev_var - 1 > 0 ) { // next variable

        v[2*x_curr] = lower;
        v[2*x_curr+1] = lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;

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

    v[2*x_curr] = lower;
    v[2*x_curr+1] = lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;

    part_points.push_back(v);

    make_posterior_logrho_vec(logrho_vec,I,level,0);
    hmap_rhos.push_back(exp(logrho_vec[n_s]));
  }

  return make_pair(part_points,hmap_rhos);
}


double CondGBT::get_add_prob(INDEX_TYPE& I,int i, int t, int level) { //get the splitting probability of dimension i
  double *node = get_node(I,level);
  double loglambda0 = (-1.0) * log((double) p);
  double *CHILD_0 = get_child(I,i,level,0);
  double *CHILD_1 = get_child(I,i,level,1);

  return exp( loglambda0 + CHILD_0[2+n_s] + CHILD_1[2+n_s] - node[1]);
}


/*
int CondGBT::find_sample_part(vector< vector< ushort > > &part_points) {

  vector< pair< INDEX_TYPE, int > >::iterator it;
  int i;
  vector<ushort> v(2*p+1);

  INDEX_TYPE I;
  int level;

  ushort x_curr = 0;
  ushort index_prev_var = 0;
  ushort lower = 0;
  ushort x_curr_count = -1;

  part_points.clear();

  for (it = sample_nodes.begin(); it < sample_nodes.end(); it++) {


    I = it->first;
    level = it->second;
    v[2*p] = level;

    for (i = 0; i < p; i++) {
      v[2*i] = 0;
      v[2*i+1] = ((ushort) 1 << k) - 1;
    }

    x_curr = 0;
    index_prev_var = 0;
    lower = 0;
    x_curr_count = -1;

    for (i = 0; i < level; i++) {
      if ( I.var[i] - index_prev_var - 1 > 0 ) { // next variable

        v[2*x_curr] = lower;
        v[2*x_curr+1] = lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;

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

    v[2*x_curr] = lower;
    v[2*x_curr+1] = lower + ((ushort) 1 << (k-x_curr_count-1)) - 1;

    part_points.push_back(v);

  }

  return part_points.size();
}
*/


void CondGBT::sample() {
  INDEX_TYPE I_root = init_index(p,0);
  sample_nodes.clear();

  sample_subtree(I_root,0,0);   // the sampled nodes are stored in sample_nodes
}

int CondGBT::update() {

  double *NODE_CURR;
  INDEX_TYPE I;

  for (int level=k; level>=0;level--) { //do it from the largest models;

    make_prior_logrho_mat(level);

    unsigned count = 0;
    I = init_index(p,level);



    while (count < modelscount[level]) {

      NODE_CURR = get_node(I,level);

      for (uint j = 0; j < pow2(level) ; j++) {

        I.var[MAXVAR] = j;
        update_node(NODE_CURR,level,I);
        NODE_CURR += numnodevar;

      }

      I = get_next_node(I,p,level); count++;
    }
  }

  return 0;
}


int CondGBT::update_node(double *NODE_CURR, int level, INDEX_TYPE I) {

  double *CHILD_0;
  double *CHILD_1;
  GBT **NODE_CURR_GBT_PTR;

  double Z[2]; // Z is a vector of length n_s_x + 1, that stores log Z(A,t,x) for the current node A
  double phi[2]; // phi is a vector of length n_s_x, that stores log phi(A,s,x) for the current node A for s=0,1,2,...,n_s_x
  // note that for s = stopping state, by design phi(A,s,x) = Z(A,s,x), which is stored in Z[n_s_x].

  double Zi;

  NODE_CURR_GBT_PTR = get_node_gbt_ptr(I,level);
  if (NODE_CURR[0] > 0) {
    if (NODE_CURR_GBT_PTR[0] == NULL) cout << "Error: Null pointer for non-empty GBT. This shouldn't happen!" << endl;

    NODE_CURR_GBT_PTR[0]->update(); // if there is at least a data point in the node

    if (level == k) { // deepest level allowed, these nodes have prior rho=1

      NODE_CURR[2] = phi[0] = phi[1] = Z[1] = NODE_CURR_GBT_PTR[0]->get_root_logphi(); // = 0
      NODE_CURR[1] = Z[0] = 0; // -DBL_MAX; // NODE_CURR[1+t] stores log Z(A,t,x) for t = 0,1,2,...,n_s_x
      NODE_CURR[3] = Z[1];	    // NODE_CURR[2+n_s_x+s] stores log phi(A,s,x) for s = 0,1,2,...,n_s_x-1

    } else if (NODE_CURR[0] <= 1) { // if the node contains no more than 1 data point

      if (NODE_CURR[0] == 1) { // if it contains one data point
        NODE_CURR[2] = phi[0] = phi[1] = Z[1] = NODE_CURR_GBT_PTR[0]->get_root_logphi();
      }

      else NODE_CURR[2] = Z[1] = 0; // if it contains no data point

      NODE_CURR[4] = NODE_CURR[3] = NODE_CURR[1] = Z[0] = Z[1];


    } else { // other models need recursion to compute phi and rho

      NODE_CURR[2] = Z[1] = NODE_CURR_GBT_PTR[0]->get_root_logphi(); // Z(A,infty,x) can be computed directly

      Z[0] = -DBL_MAX;

      for (int i=0; i < p; i++) { // for each dimension i

        CHILD_0 = get_child(I,i,level,0);
        CHILD_1 = get_child(I,i,level,1);

        Zi = loglambda0 + CHILD_0[3] + CHILD_1[3];

        if (Z[0] == -DBL_MAX) {
          Z[0] = Zi;
        } else {
           Z[0] = log_exp_x_plus_exp_y (Z[0], Zi);
        }

      }

      // Now Z[t] stores Z(A,t,x) for all states t
      NODE_CURR[1] = Z[0]; // NODE[1+t] stores Z(A,t,x)

      // Next we compute the phi(A,s,x)
      phi[0] = logrho_mat[0][0] + Z[0];
      phi[0] = log_exp_x_plus_exp_y(phi[0],logrho_mat[0][1]+Z[1]);
      phi[1] = Z[1]; // phi(A,infty,x) = Z(A,infty,x)
      NODE_CURR[3] = phi[0];

    }

    NODE_CURR[3] = phi[0];
    NODE_CURR[4] = phi[1];
  }

  else { // if NODE_CURR contains no data, then everything is zero
    NODE_CURR[2] = phi[0] = phi[1] = Z[1] = 0;
    NODE_CURR[1] = Z[0] = 0;
    NODE_CURR[3] = Z[1];
    NODE_CURR[4] = 0;
  }


  return 0;
}




/* The following implementation of CondGBT::update_node should match the results above */
/*
int CondGBT::update_node(double *NODE_CURR, int level, INDEX_TYPE I) {

  double *CHILD_0;
  double *CHILD_1;
  GBT ** NODE_CURR_GBT_PTR;

  NODE_CURR_GBT_PTR = get_node_gbt_ptr(I,level);
  NODE_CURR_GBT_PTR[0]->update();

  double Z[n_s+1]; // Z is a vector of length n_s + 1, that stores log Z(A,t,x) for the current node A
  double phi[n_s+1]; // phi is a vector of length n_s, that stores log phi(A,s,x) for the current node A for s=0,1,2,...,n_s
  // note that for s = stopping state, by design phi(A,s,x) = Z(A,s,x), which is stored in Z[n_s].

  double Zi;
  int t,s; // nonstopping state

  if (NODE_CURR[0] == 0) {
    for (s=1; s<= 2+2*n_s; s++) {
      NODE_CURR[s] = 0;
    }
  }

  else {
    if (NODE_CURR_GBT_PTR[0] == NULL) cout << "Error: Null pointer for non-empty GBT. This shouldn't happen!" << endl;

    if (level == k) { // deepest level allowed, these nodes have prior rho=1

      NODE_CURR[1+n_s] = phi[n_s] = Z[n_s] = NODE_CURR_GBT_PTR[0]->get_root_logphi(); // = 0
      // NODE_CURR[1+n_s] = phi[n_s] = Z[n_s] = NODE_CURR[0]*(level-k)*log(2.0); // = 0
      // NODE_CURR[1+t] stores log Z(A,t,x) for t = 0,1,2,...,n_s
      // NODE_CURR[2+n_s+s] stores log phi(A,s,x) for s = 0,1,2,...,n_s-1

      for (s=0; s<n_s; s++) {
        t = s;
        NODE_CURR[1+t] = Z[t] = 0; // -DBL_MAX; // NODE_CURR[1+t] stores log Z(A,t,x) for t = 0,1,2,...,n_s
        // In this case Z(A,t,x) for t nonstopping is -inf
        NODE_CURR[2+n_s+s] = Z[n_s];	    // NODE_CURR[2+n_s+s] stores log phi(A,s,x) for s = 0,1,2,...,n_s-1
      }

    } else if (NODE_CURR[0] <= 1) { // if the node contains no more than 1 data point


        if (NODE_CURR[0] == 1) { // if it contains one data point
          NODE_CURR[1+n_s] = Z[n_s] = NODE_CURR_GBT_PTR[0]->get_root_logphi();
        }

        else NODE_CURR[1+n_s] = Z[n_s] = 0; // if it contains no data point

        for (t=0; t< n_s; t++) {
          // in either case, phi(A,s,x) and Z(A,t,x) all correspond to the uniform likelihood
          s=t;
          NODE_CURR[2+n_s+s] = NODE_CURR[1+t] = Z[t] = Z[n_s];
        }

    } else { // other models need recursion to compute phi and rho

        NODE_CURR[1+n_s] = Z[n_s] = NODE_CURR_GBT_PTR[0]->get_root_logphi(); // Z(A,infty,x) can be computed directly

        for (t=0; t < n_s; t++) { // non-stopping states

          Z[t] = -DBL_MAX;

          for (int i=0; i < p; i++) { // for each dimension i

            CHILD_0 = get_child(I,i,level,0);
            CHILD_1 = get_child(I,i,level,1);


            Zi = loglambda0 + CHILD_0[2+n_s+t] + CHILD_1[2+n_s+t];


            if (Z[t] == -DBL_MAX) {
              Z[t] = Zi;
            } else {
              Z[t] = log_exp_x_plus_exp_y (Z[t], Zi);
            }
          }


          // Now Z[t] stores Z(A,t,x) for all states t
          NODE_CURR[1+t] = Z[t]; // NODE[1+t] stores Z(A,t,x)


          // Next we compute the phi(A,s,x)

          for (s = 0; s < n_s; s++) { // for the non-stopping states

            if (t == 0) {
              phi[s] = logrho_mat[s][t] + Z[t];
            }

            else {
              phi[s] = log_exp_x_plus_exp_y(phi[s], logrho_mat[s][t] + Z[t]);
            }
          }
        }

        for (s=0; s < n_s; s++) {
          phi[s] = log_exp_x_plus_exp_y(phi[s],logrho_mat[s][n_s]+Z[n_s]);
        }

        phi[n_s] = Z[n_s]; // phi(A,infty,x) = Z(A,infty,x)

        for (s=0; s<=n_s; s++) {

          NODE_CURR[2+n_s+s] = phi[s];
        }
    }
  }

  return 0;
}
*/

void CondGBT::add_data_to_subtree(INDEX_TYPE I, int level, int x_curr, int part_count,
                                  Col< unsigned int > obs_X, Col< unsigned int > obs_Y) {

  double *NODE_CURR;
  GBT ** NODE_CURR_GBT_PTR;

  INDEX_TYPE I_child;


  NODE_CURR = get_node(I,level);
  NODE_CURR_GBT_PTR = get_node_gbt_ptr(I,level);

  NODE_CURR[0] += 1;

  INDEX_TYPE I_GBT_root = init_index(p_Y,0);

  if (NODE_CURR[0] == 1) {

    NODE_CURR_GBT_PTR[0] = new GBT(conv_to< Mat< unsigned int> >::from(obs_Y.t()),k_Y,p_Y,rho0_Y,rho0_mode_Y,tran_mode,
                                   lognu_lowerbound,lognu_upperbound,n_grid,n_s,beta); // need to specify initial values
  }

  else {

    NODE_CURR_GBT_PTR[0]->add_data_to_subtree(I_GBT_root,0,1,0,obs_Y);
    (NODE_CURR_GBT_PTR[0]->nobs)++;

  }

  int i = 0;

  if (level < k) { // not maximum (leaf) node

    i = x_curr - 1;
    I_child = make_child_index(I, i, level, (obs_X(i) >> part_count) & 1);
    add_data_to_subtree(I_child, level+1, x_curr, part_count+1, obs_X, obs_Y);


    for (i = x_curr; i < p; i++) {
      I_child = make_child_index(I,i,level,obs_X(i) & 1);
      add_data_to_subtree(I_child, level+1, i+1, 1, obs_X, obs_Y);
    }
  }
}



void CondGBT::remove_data_from_subtree(INDEX_TYPE I, int level, int x_curr, int part_count, Col< unsigned int > obs_X, Col< unsigned int > obs_Y) {
  double *NODE_CURR;
  GBT ** NODE_CURR_GBT_PTR;

  INDEX_TYPE I_child;
  INDEX_TYPE I_GBT = init_index(p_Y,0);

  NODE_CURR = get_node(I,level);
  NODE_CURR[0] -= 1;

  NODE_CURR_GBT_PTR = get_node_gbt_ptr(I,level);
  NODE_CURR_GBT_PTR[0]->remove_data_from_subtree(I_GBT,0,1,0,obs_Y);

  int i = 0;

  if (level < k) { // not maximum (leaf) node

    i = x_curr - 1;
    I_child = make_child_index(I, i, level, (obs_X(i) >> part_count) & 1);
    remove_data_from_subtree(I_child, level+1, x_curr, part_count+1, obs_X, obs_Y);


    for (i = x_curr; i < p; i++) {
      I_child = make_child_index(I,i,level,obs_X(i) & 1);
      remove_data_from_subtree(I_child, level+1, i+1, 1, obs_X, obs_Y);
    }
  }
}

int CondGBT::update_subtree_add_new_data(INDEX_TYPE I, int level, int x_curr, int part_count, Col< unsigned int > new_obs_X, Col< unsigned int > new_obs_Y) {
  double *NODE_CURR;
  GBT ** NODE_CURR_GBT_PTR;

  INDEX_TYPE I_child;
  INDEX_TYPE I_GBT = init_index(p_Y,0);

  NODE_CURR = get_node(I,level);
  NODE_CURR[0] += 1;

  NODE_CURR_GBT_PTR = get_node_gbt_ptr(I,level);
  if (NODE_CURR_GBT_PTR[0] != NULL) {
    NODE_CURR_GBT_PTR[0]->update_subtree_add_new_data(I_GBT,0,1,0,new_obs_Y);
  }
  else {
    NODE_CURR_GBT_PTR[0] = new GBT(conv_to< Mat< unsigned int> >::from(new_obs_Y.t()),k_Y,p_Y,rho0_Y,rho0_mode_Y,tran_mode,
                                   lognu_lowerbound,lognu_upperbound,n_grid,n_s,beta);
    NODE_CURR_GBT_PTR[0]->update();
  }

  int i;

  if (level < k) { // not maximum (leaf) node

    i = x_curr - 1;
    I_child = make_child_index(I, i, level, (new_obs_X(i) >> part_count) & 1);
    update_subtree_add_new_data(I_child, level+1, x_curr, part_count+1, new_obs_X, new_obs_Y);


    for (i = x_curr; i < p; i++) {

      I_child = make_child_index(I,i,level,new_obs_X(i) & 1);
      update_subtree_add_new_data(I_child, level+1, i+1, 1, new_obs_X, new_obs_Y);
    }

  }

  make_prior_logrho_mat(level);
  update_node(NODE_CURR,level,I);

  return 0;
}


int CondGBT::update_subtree_remove_new_data(INDEX_TYPE I, int level, int x_curr, int part_count, Col< unsigned int > new_obs_X, Col< unsigned int > new_obs_Y) {
  double *NODE_CURR;
  GBT ** NODE_CURR_GBT_PTR;

  INDEX_TYPE I_child;
  INDEX_TYPE I_GBT = init_index(p_Y,0);

  NODE_CURR = get_node(I,level);
  NODE_CURR[0] -= 1;

  NODE_CURR_GBT_PTR = get_node_gbt_ptr(I,level);
  NODE_CURR_GBT_PTR[0]->update_subtree_remove_new_data(I_GBT,0,1,0,new_obs_Y);

  int i;

  if (level < k) { // not maximum (leaf) node

    i = x_curr - 1;
    I_child = make_child_index(I, i, level, (new_obs_X(i) >> part_count) & 1);
    update_subtree_remove_new_data(I_child, level+1, x_curr, part_count+1, new_obs_X, new_obs_Y);


    for (i = x_curr; i < p; i++) {

      I_child = make_child_index(I,i,level,new_obs_X(i) & 1);
      update_subtree_remove_new_data(I_child, level+1, i+1, 1, new_obs_X, new_obs_Y);
    }

  }

  make_prior_logrho_mat(level);
  update_node(NODE_CURR,level,I);

  return 0;

}


vector<double> CondGBT::compute_predictive_density(Mat< unsigned int >Xnew, Mat< unsigned int > Ynew) {

  int nobs_new = Xnew.n_rows;

  INDEX_TYPE I_root = init_index(p,0);
  vector<double> predictive_density(nobs_new);
  double new_root_logphi;
  double original_root_logphi;

  original_root_logphi = get_root_logphi();

  for (int i=0; i < nobs_new; i++) { //since we do it only once, do this for now

    update_subtree_add_new_data(I_root,0,1,0,Xnew.row(i).t(),Ynew.row(i).t());
    new_root_logphi = get_root_logphi();

    predictive_density[i] = exp(new_root_logphi - original_root_logphi);
    update_subtree_remove_new_data(I_root,0,1,0,Xnew.row(i).t(),Ynew.row(i).t());
  }

  return predictive_density;

}
