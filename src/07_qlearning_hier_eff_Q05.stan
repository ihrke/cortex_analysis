functions {
  real qlearning_log(int[] ACC, real alpha, real beta, int ntrials, int[] pair, int[] reward) {
    vector[6] Q;
    real LL;
    real logPt;
    int a;
    int A;
    int B;
    real pe; // prediction error
    real offset; // for logsumexp trick
      
    Q <- rep_vector(0.5, 6);
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
        Q[a] <- Q[a] + (alpha*pe);

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
  real alphai[nsubj]; // alphai keeps track of the group-level-derived individual logit(alpha) pars
  real logbeta[nsubj];

  real mu_alpha; // group-level parameters, mu/sigma
  real<lower=0> sig_alpha; 
  
  real mu_beta; // group-level parameters, mu/sigma
  real<lower=0> sig_beta; 
  
  real effBalpha;
  real effBbeta;
}
transformed parameters {
  real<lower=0,upper=1> alphaA[nsubj];// inv-logit transform the alphai parameters 
  real<lower=0,upper=1> alphaB[nsubj]; 
  real<lower=0> betaA[nsubj];
  real<lower=0> betaB[nsubj];

  for(i in 1:nsubj){
    alphaA[i] <- inv_logit(alphai[i]);
    alphaB[i] <- inv_logit(alphai[i]+effBalpha); // effect is on the logit-scale
    
    betaA[i] <- exp(logbeta[i]);
    betaB[i] <- exp(logbeta[i] + effBbeta);
  }
  
}
model {
  

  mu_beta ~ normal( 0, 100 );
  sig_beta ~ uniform( 0, 100 );
  
  mu_alpha ~ normal( 0, 100 );
  sig_alpha ~ uniform( 0, 100 );

  effBalpha ~ normal(0,1);
  effBbeta ~ normal(0,1);

  alphai ~ normal( mu_alpha, sig_alpha); 
  logbeta ~ normal( mu_beta, sig_beta );
  
  for( i in 1:nsubj ){
    ACCA[i] ~ qlearning(alphaA[i], betaA[i], ntrials, pairA[i], rewardA[i]);
    ACCB[i] ~ qlearning(alphaB[i], betaB[i], ntrials, pairB[i], rewardB[i]);
  }
}

generated quantities {
  vector[nsubj] log_lik;
  
  for (j in 1:nsubj){
    log_lik[j]<-qlearning_log(ACCA[j], alphaA[j], betaA[j], ntrials, pairA[j], rewardA[j]) + 
                qlearning_log(ACCB[j], alphaB[j], betaB[j], ntrials, pairB[j], rewardB[j]);
  }
}