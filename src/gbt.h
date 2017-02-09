#ifndef GBT_H
#define GBT_H
#define NUMNODEVAR 50
class GBT {  
public:
    //data dimmensions
    int nobs;   //number of observations
    int p;      //total number of predictors
    
    //model parameters
    int k; // maximum model size
    double rho0; // stopping probability at level 0 (sample space)
    int n_s;
    int rho_mode,tran_mode,alpha_mode;
    double beta; // transition probabilities
    double lognu_lowerbound,lognu_upperbound;
    double b; // split point on lognu for weak vs strong shrinkage

    int n_grid;
    double loglambda0;
    double **logrho_mat;
    double *logrho_vec;

    //constructors
    GBT(uint *Xwork, int nobs, int k, int p, double rho0, int rho_mode, int tran_mode, double lognu_lowerbound, double lognu_upperbound, int n_grid, int n_s, double beta);
    ~GBT();
    
    //main functions
    int update_node(double *, int, INDEX_TYPE);
    int update();

    //getters
    double get_root_logrho();
    double get_root_logphi();
    double get_log_Ma(double theta0,int n_0,int n_1,int t);
    double make_alpha0(int t, int level);
    double make_rho0(int level);
    void make_prior_logrho_mat(int level);
    void make_prior_logrho_vec(double *logrho_vec,int level, int s);
    void make_posterior_logrho_vec(double *logrho_vec,INDEX_TYPE& I, int level, int s);
    double get_node_logphi(vector<double> &node_INFO);

    void sample();
    void find_hmap(int print);
    vector< vector<ushort> > find_part();
    int find_sample_part(vector< vector< ushort> > &part_points, vector< vector< double> > &nu_and_probs);
    
    vector<double> compute_predictive_density(uint *Xnew, int nobs_new);

private:
    //data
     
    //working variables
    double **models;
    unsigned int *modelscount;
    
    vector< pair< pair< INDEX_TYPE, int >, pair< double,double > > > sample_nodes;
    vector< pair< INDEX_TYPE, pair<int,int> > > hmap_nodes;

    void init(uint *Xwork);
    void clear();
    
    void sample_subtree(INDEX_TYPE& I,int level,int s, double prob);
    void find_hmap_subtree(INDEX_TYPE& I, int level);

    double get_add_prob(INDEX_TYPE& I,int i, int t, int level);
    double *get_child(INDEX_TYPE& I, int i,int level,ushort which);
    double *get_node(INDEX_TYPE& I, int level);
    
    int convert_obs_to_uint(int *obs,INDEX_TYPE I,int level);
    void add_data_to_subtree(INDEX_TYPE I, int level, int x_curr, int part_count, uint *obs);
    void remove_data_from_subtree(INDEX_TYPE I, int level, int x_curr, int part_count, uint *obs);
    int update_subtree_add_new_data(INDEX_TYPE I, int level, int x_curr, int part_count, uint *new_obs);
    int update_subtree_remove_new_data(INDEX_TYPE I, int level, int x_curr, int part_count, uint *new_obs);
  

};
#endif
