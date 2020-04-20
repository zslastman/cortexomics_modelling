functions {
  vector log_minus_exp(vector a, vector b) {
   return log(exp(a - b) -1 )+ b;
  }
  vector log_sum_exp(vector a, vector b){
  return log(exp(a - b) + 1) + b;
  }

}


data {
  int G;// number of proteins
  int T;// info on the number of conditions
  matrix[G,T] lMSmu;
  matrix[G,T] lSeqmu;
  matrix[T,T] lMSsigma[G];
  matrix[T,T] lSeqsigma[G];
  real l_st_priorsd;
  real l_ribo_priorsd;
  real l_ribo_priornu;
  real l_pihalf_priormu;
  real l_pihalf_priorsd;
}

parameters {
  vector<lower=-10,upper=10>[G] l_st; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lribo;  // log vector of ribo-seq levels
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector[G] prot0; // initial LOG amount of protein
}

transformed parameters{
    matrix [G,T] prot; // amounts of protein
    vector[G] Kd; // the degred
    vector[G] lKd; // the log degrad
    vector[G] lKs; // the (log) synthesis constant
    //get Kd
    lKd = log(log(2)) -  l_pihalf;
  
    //get Ks
    lKs = l_st -  lKd;
  
    prot[,1] = prot0;
    for(i in 2:T){
      // put thisinto wolfram  exp(  - d_0  T  )  * (P0 +   Integrate [  ( s_0 + (s_1 - s_0) t)  * exp( d_0 (t) )  dt , from 0 to T ] )
      //this is correct, but stans functions don't wory this way, we need to accumulate
     // prot[,i] = logSumMinus(
      //    log_sum_exp(prot[,i-1] - exp(lKd)  ,  lribo[,i-1] -2 *lKd  ,  lribo[,i]-lKd  ,-exp(lKd)+lribo[,i]-2*lKd ),
     //     log_sum_exp(lribo[,i]-2*lKd , -exp(lKd)+lribo[,i-1]-lKd , -exp(lKd)+lribo[,i-1]-2*lKd)
     // )
      prot[,i] = prot[,i-1] - exp(lKd);//1
      prot[,i] = log_sum_exp(prot[,i],lKs+ lribo[,i]-lKd);//5
      prot[,i] = log_sum_exp(prot[,i],lKs+-exp(lKd)+lribo[,i]-2*lKd );//6
    
      prot[,i] = log_sum_exp(prot[,i],lKs+lribo[,i-1]-lKd);//7
      prot[,i] = log_sum_exp(prot[,i],lKs+lribo[,i-1] -2 *lKd );//8
      prot[,i] = log_minus_exp(prot[,i],lKs+lribo[,i]-2*lKd );//9
      prot[,i] = log_minus_exp(prot[,i],lKs-exp(lKd)+lribo[,i-1] - lKd );//2
      prot[,i] = log_minus_exp(prot[,i],lKs+lribo[,i-1] - lKd );//3
      prot[,i] = log_minus_exp(prot[,i],lKs+-exp(lKd)+lribo[,i-1] - 2*lKd );//4

    }
}

model {
  // this needs to become steady state I think
  // for(c in 1:ncond){
    // l_st[,c] ~ normal(mu0, sqrt(sigma20));

  l_st ~ normal(0,l_st_priorsd);
  l_pihalf ~ normal(l_pihalf_priormu,l_pihalf_priorsd);
  
  //for(t in 1:T){
  //  lribo[,t] ~ student_t(l_ribo_priornu,0,l_ribo_priorsd);
  //}// put a prior on fold changes

  for(g in 1:G){
  
    lSeqmu[g] ~ multi_normal(lribo[g],lSeqsigma[g]);

    lMSmu[g]  ~ multi_normal(prot[g],lMSsigma[g]);
    
  }
  
}



generated quantities{
    matrix[G,T] lMSdev;
    lMSdev = prot - lMSmu;
}

