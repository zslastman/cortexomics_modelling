data {
  int nsamples;// number of samples
  int nribosamples;// number of samples
  int G;// number of proteins
  int T;// info on the number of conditions
  int totalmissing;// info on th total amount of missing data
  matrix[G, nsamples] lMS;// data
  matrix[G, nribosamples] lribo;// data
  matrix[G, nribosamples] ribo_sigma;  
  int experimental_design[nsamples];// indicator variable matching the samples to conditions
  int experimental_design_r[nribosamples];
  real zeta[nsamples];// the spread of the dropout point for a library, gets combined with the varianc per protein
  real rho[nsamples];// rho the location of the dropout point for a given library
  real mu0;// fit by proDD the mean of means
  real sigma20;// fit by proDD - the variance in means
  real eta;// fit by proDD - th evariance in protein variances
  real nu;// fit by proDD - the mean of protein variances
  matrix[T,T+1] mybs; // spline basis for the protein trajectory
  matrix[T,T] mydbs; // differentiated spline basis for the RNA (piecewise linear right now)
  real l_st_priorsd;
  real l_ribo_priorsd;
  real l_ribo_priornu;
  real l_pihalf_priormu;
  real l_pihalf_priorsd;
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector[T] ribnorm; // normalization factor for the riboseq
  vector[T] protnorm;  // normalization factor for the protein
}

parameters {
  real<lower=0> sigma2[G];// the variance for that protein
  vector<lower=-10,upper=10>[G] l_st; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lsynths;  // log vector of fold changes due to synthesis
  vector[G] prot0; // initial LOG amount of protein
}

transformed parameters{
    matrix [G,T] prot; // amounts of protein
    matrix [G,T] mRNA; // amounts of mRNA
    vector[G] Kd; // the degred
    vector[G] lKd; // the log degrad
    vector[G] lKs; // the (log) synthesis constant
    real zetastar[totalmissing];
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
    
    prot = append_col(prot0, fcs ) * (mybs)' - (rep_matrix(protnorm,G)'); // get the full protein trajectory

    mRNA = prot + log2(synth * (mydbs)' ) - rep_matrix(lKs,T) - (rep_matrix(ribnorm,G)'); // get the mRNA trajectory this implies


  {
    int counter = 1;
    for(i in 1:G){
      for(j in 1:nsamples){
        if(is_inf(lMS[i, j])){
          zetastar[counter] = zeta[j] * sqrt(1 + sigma2[i]/zeta[j]^2);
          counter = counter + 1;
        }
      }
    }
  }
}

model {
  // this needs to become steady state I think
  // for(c in 1:ncond){
    // l_st[,c] ~ normal(mu0, sqrt(sigma20));
  prot0 ~ normal(mu0, sqrt(sigma20));
  sigma2 ~ scaled_inv_chi_square(nu, sqrt(eta));

  l_st ~ normal(0,l_st_priorsd);
  l_pihalf ~ normal(l_pihalf_priormu,l_pihalf_priorsd);
  for(t in 1:T){
    lsynths[,t] ~ student_t(l_ribo_priornu,0,l_ribo_priorsd);
  }// put a prior on fold changes

  // }
  #prior disribution on the differences in protein 
  ribnorm ~ normal(0,3);
  protnorm ~ normal(0,3);
  
  {
    int counter = 1;
    for(i in 1:G){
      for(jr in 1:nribosamples){
        lribo[i,jr] ~ normal(mRNA[i,experimental_design_r[jr]], ribo_sigma[i,jr]); // use confidence intervals to weight the observations
      }
      for(j in 1:nsamples){
        if(is_inf(lMS[i, j])){
          target += normal_lccdf(prot[i, experimental_design[j]] | rho[j], fabs(zetastar[counter]));
          counter += 1;
        }else{
          lMS[i, j] ~ normal(prot[i, experimental_design[j]], sqrt(sigma2[i]));
        }
      }
    }
  }
}
