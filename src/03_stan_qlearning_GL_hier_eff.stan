functions {
  real qlearningpn_log(int[] ACC, real alphaG, real alphaL, real beta, int ntrials, int[] pair, int[] reward) {
    vector[6] Q;
    real LL;
    real logPt;
    int a;
    int A;
    int B;
    real pe; // prediction error
    real offset; // for logsumexp trick
      
    Q <- rep_vector(0, 6);
    LL <- 0.0;
    for( t in 1:ntrials ){
      if(ACC[t]>=0){ // missing values are skipped
        A <- 2*(pair[t])-1;
        B <- A+1;
        if(ACC[t]==1){
          a <- A;
        } else {
          a <- B;
        }
        pe <- (reward[t]-Q[a]);
        Q[a] <- Q[a] + (if_else(pe>0,alphaG,alphaL)*pe);

        offset <- if_else( (Q[A]/beta)>(Q[B]/beta), (Q[A]/beta), (Q[B]/beta) );
        
        logPt <- Q[a]/beta - (log( exp( Q[A]/beta - offset) + exp( Q[B]/beta - offset)) + offset);
        LL <- LL+logPt;
      }
    }
    //print("a=", alpha, " b=", beta, " ntrials=", ntrials, " Q=", Q, " LL=", LL);

    return LL;
  }
}
data {
  int<lower=0> ntrials;
  int<lower=0> nsubj;

  int<lower=-1,upper=1> ACCA[nsubj,ntrials];
  int<lower=1,upper=3> pairA[nsubj,ntrials];
  int<lower=0,upper=1> rewardA[nsubj,ntrials];

  int<lower=-1,upper=1> ACCB[nsubj,ntrials];
  int<lower=1,upper=3> pairB[nsubj,ntrials];
  int<lower=0,upper=1> rewardB[nsubj,ntrials];

}

parameters {
  //real<lower=0,upper=1> alpha[nsubj]; # individual QL parameters
  real alphaGi[nsubj]; // alphai keeps track of the group-level-derived individual logit(alpha) pars
  real alphaLi[nsubj];
  real logbeta[nsubj];

  real mu_alphaG; // group-level parameters, mu/sigma
  real mu_alphaL;
  real<lower=0> sig_alpha; 
  
  real mu_beta; // group-level parameters, mu/sigma
  real<lower=0> sig_beta; 
  
  // effects
  real effBalphaG;
  real effBalphaL;
  real effBbeta;
}
transformed parameters {
  real<lower=0,upper=1> alphaGA[nsubj]; // inv-logit transform the alphai parameters 
  real<lower=0,upper=1> alphaGB[nsubj]; 
  real<lower=0,upper=1> alphaLA[nsubj]; 
  real<lower=0,upper=1> alphaLB[nsubj]; 
  real<lower=0> betaA[nsubj];
  real<lower=0> betaB[nsubj];

  for(i in 1:nsubj){
    alphaGA[i]  <- inv_logit(alphaGi[i]);
    alphaGB[i] <- inv_logit(alphaGi[i]+effBalphaG);
    alphaLA[i]  <- inv_logit(alphaLi[i]);
    alphaLB[i] <- inv_logit(alphaLi[i]+effBalphaL);
    
    betaA[i] <- exp(logbeta[i]);
    betaB[i] <- exp(logbeta[i] + effBbeta);
  }
  
}
model {

  mu_beta ~ normal( 0, 100 );
  sig_beta ~ uniform( 0, 100 );
  
  mu_alphaG ~ normal( 0, 100 );
  mu_alphaL ~ normal( 0, 100 );
  sig_alpha ~ uniform( 0, 100 );

  effBalphaG ~ normal(0,1);
  effBalphaL ~ normal(0,1);
  effBbeta ~ normal(0,1);

  alphaGi ~ normal( mu_alphaG, sig_alpha);
  alphaLi ~ normal( mu_alphaL, sig_alpha); 
  logbeta ~ normal( mu_beta, sig_beta );
  
  for( i in 1:nsubj ){
    ACCA[i] ~ qlearningpn(alphaGA[i], alphaLA[i], betaA[i], ntrials, pairA[i], rewardA[i]);
    ACCB[i] ~ qlearningpn(alphaGB[i], alphaLB[i], betaB[i], ntrials, pairB[i], rewardB[i]);
  }
}

generated quantities {
  vector[nsubj] log_lik;
  
  for (j in 1:nsubj){
    log_lik[j]<-qlearningpn_log(ACCA[j], alphaGA[j], alphaLA[j], betaA[j], ntrials, pairA[j], rewardA[j]) + 
                qlearningpn_log(ACCB[j], alphaGB[j], alphaLB[j], betaB[j], ntrials, pairB[j], rewardB[j]);
  }
}