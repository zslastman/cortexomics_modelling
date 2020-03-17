data {
  int nsamples;// number of samples
  int nribosamples;// number of samples
  int G;// number of proteins
  int T;// info on the number of conditions
  int totalmissing;// info on th total amount of missing data
  matrix[G, nsamples] lMS;// data
  matrix[G, nribosamples] lribo;// data
  matrix[G, nribosamples] voom_sigma;  
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
  real l_pihalf_priorsd;
}

parameters {
  real<lower=0> sigma2[G];// the variance for that protein
  vector<lower=-10,upper=10>[G] l_st; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lsynths;  // log vector of fold changes due to synthesis
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
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

    prot[0] = prot0
    
    for(i in 2:T){
      #prot[1] = prot[i-1]*exp(-Kd *t) - (s1 / Kd^2) + ( (s1 *t) / Kd) + (s0/Kd)
      #prot[1] = prot[i-1]*exp(-Kd) - (s1 / Kd^2) + ( (s1) / Kd) + (s0/Kd)#eliminate T
      log(prot[1]) = log(prot[i-1]*exp(-Kd) - (s1 / Kd^2) +  (s1 / Kd) + (s0/Kd) )#eliminate T
      #or in terms of timeopints
      log(prot[1]) = log(prot[i-1])*exp(-Kd) + (snew/Kd) - ( snew - sold / Kd^2)   #eliminate T
      
      #we can break this into two pieces for logsumExp
      log(prot[i-1]*exp(-Kd))#which becomes
      log(prot[i-1]) -Kd
      #and
      lKs - ((mRNA[i]-mRNA[i-1])/Kd) + mRNA[i]+mRNA[i-1]  - lKd
      
      
    }





    
    prot[i] = log_sum_exp(
      prot[i-1]-lKd, #The amount left under degredation
      lKs - ((mRNA[i]-mRNA[i-1])/Kd) + mRNA[i]+mRNA[i-1]  - lKd #the amount synthesized
    )
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
  l_pihalf ~ normal(0,l_pihalf_priorsd);
  for(t in 1:T){
    lsynths[,t] ~ normal(0,l_ribo_priorsd);
  }// put a prior on fold changes

  // }


  {
    int counter = 1;
    for(i in 1:G){
      for(jr in 1:nribosamples){
        lribo[i,jr] ~ normal(mRNA[i,experimental_design_r[jr]], voom_sigma[i,jr]); // use confidence intervals to weight the observations
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

