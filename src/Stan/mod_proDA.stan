data {
  int nsamples;// number of samples
  int nribosamples;// number of samples
  int G;// number of proteins
  int T;// info on the number of conditions
  int totalmissing;// info on th total amount of missing data
  matrix[G,T] lMSmu;
  matrix[G,T] lSeqmu;
  matrix[G,T] lMSsigma[T];
  matrix[G ,T] lSeqsigma[T];
  matrix[T ,T+1] mybs; // spline basis for the protein trajectory
  matrix[T,T] mydbs; // differentiated spline basis for the RNA (piecewise linear right now)
  real l_st_priorsd;
  real l_ribo_priorsd;
  real l_ribo_priornu;
  real l_pihalf_priormu;
  real l_pihalf_priorsd;
}

parameters {
  vector<lower=-10,upper=10>[G] l_st; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lsynths;  // log vector of fold changes due to synthesis
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector[G] prot0; // initial LOG amount of protein
}

transformed parameters{
    matrix [G,T] prot; // amounts of protein
    matrix [G,T] seq; // amounts of mRNA
    vector[G] Kd; // the degred
    vector[G] lKd; // the log degrad
    vector[G] lKs; // the (log) synthesis constant
    matrix [G,T] synth; // vector of fold changes due to synthesis
    matrix [G,T] fcs; // vector of fold changes due to synthesis

    //get Kd
    lKd = log2(log(2)) -  l_pihalf;
    Kd = exp(log(2)*lKd);
    
    //get Ks
    lKs = l_st -  lKd;

    //get the fold changes due to synthesis (which must be positive)
    synth = exp(log(2)*lsynths);
    
    fcs = synth - rep_matrix( Kd,T);
    
    prot = append_col(prot0, fcs ) * (mybs)'; // get the full protein trajectory

    seq = prot + log2(synth * (mydbs)' ) - rep_matrix(lKs,T); // get the mRNA trajectory this implies

}

model {
  // this needs to become steady state I think
  // for(c in 1:ncond){
    // l_st[,c] ~ normal(mu0, sqrt(sigma20));

  l_st ~ normal(0,l_st_priorsd);
  l_pihalf ~ normal(l_pihalf_priormu,l_pihalf_priorsd);
  
  for(t in 1:T){
    lsynths[,t] ~ student_t(l_ribo_priornu,0,l_ribo_priorsd);
  }// put a prior on fold changes

  for(g in 1:G){
    seq[g,] ~ multi_normal(lSeqmu[g,],lSeqsigma[g]);
    prot[g,] ~ multi_normal(lMSmu[g,],lMSsigma[g]);
  }
  
}


