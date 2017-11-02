##
## model with effects on logit/log scale
##
library(ProjectTemplate)
load.project()
invlogit=arm::invlogit
library(rstan)
theme_set(theme_bw())
rstan_options(auto_write=TRUE) # for parallel

# enabls --force for stallo etc.
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options)
  uncache.all()
bname<-tools::file_path_sans_ext(basename(this.file.name()))
mod.fname=sprintf("./src/%s.stan", bname)

n.chains=2
n.cores=2
n.iter=2000
n.warmup=1000

data=list(
  nsubj=length(unique(learn$subj)),
  ntrials=with(learn, sum(subj==subj[1] & condition==condition[1])),
  ACCA=(learn %>% filter(condition=="A") %>% 
    select(subj, ACC) %>% group_by(subj) %>% mutate(trial=1:n()) %>% 
    ungroup %>% spread( subj, ACC ) %>% select(-trial) %>% as.matrix %>% t),
  pairA=(learn %>% filter(condition=="A") %>% 
           select(subj, pair) %>% group_by(subj) %>% mutate(trial=1:n()) %>% 
           ungroup %>% spread( subj, pair ) %>% select(-trial) %>% as.matrix %>% t),
  rewardA=(learn %>% filter(condition=="A") %>% 
               select(subj, reward) %>% group_by(subj) %>% mutate(trial=1:n()) %>% 
               ungroup %>% spread( subj, reward ) %>% select(-trial) %>% as.matrix %>% t),
  ACCB=(learn %>% filter(condition=="B") %>% 
          select(subj, ACC) %>% group_by(subj) %>% mutate(trial=1:n()) %>% 
          ungroup %>% spread( subj, ACC ) %>% select(-trial) %>% as.matrix %>% t),
  pairB=(learn %>% filter(condition=="B") %>% 
           select(subj, pair) %>% group_by(subj) %>% mutate(trial=1:n()) %>% 
           ungroup %>% spread( subj, pair ) %>% select(-trial) %>% as.matrix %>% t),
  rewardB=(learn %>% filter(condition=="B") %>% 
             select(subj, reward) %>% group_by(subj) %>% mutate(trial=1:n()) %>% 
             ungroup %>% spread( subj, reward ) %>% select(-trial) %>% as.matrix %>% t)
  
)


## init
nsubj=data$nsubj
initfct <- function(chainix){
  nsubj=data$nsubj
  list(alphai=arm::logit(runif(nsubj, 0.01, 0.2)),
       logbeta=log(runif(nsubj,0.01, 0.4)),
       mu_alpha=rnorm(1,0,0.1),
       sig_alpha=0.2,
       mu_beta=rnorm(1,0,0.1),
       sig_beta=0.2,
       effBalpha=rnorm(1,0,0.1),
       effBbeta=rnorm(1,0,0.1)
  )
}
#stop()
if(is.cached.var("fit")){
  printf("WARNING: loading variables from cache\n")
  fit=load.cache.var("fit")
} else {
  start.time = proc.time()
  fit = stan(file=mod.fname, data = data, chains = n.chains, iter = n.iter, warmup = n.warmup, init=initfct, cores=n.cores)
  tellmedone()
  end.time = proc.time()
  show("Time taken for sampling")
  show(end.time-start.time)
  
  cache.var('fit')
}

stop()
print(fit)  

#plot(fit)
m <- as.matrix(fit)
vars=colnames(m)

traceplot(fit, pars=c("alphaA", "alphaB"))
plot(fit, pars=c("alphaA", "alphaB"))
plot(fit, pars=c("betaA"))

plot(fit, pars=vars[str_detect(vars,"eff")])

#waic(fit)




means <- m[,str_detect(vars, 'mu')]
dotplot(means)+ggtitle("means on logit/log scales")#+xlim(0,.3)
means.org <- cbind( alpha=invlogit(means[,'mu_alpha']), beta=exp(means[,'mu_beta']) )
dotplot(means.org)+ggtitle('means on original scales')+xlim(0,.5)

dotplot(m[,str_detect(vars, 'eff')])+ggtitle("Effect on logit scale")+geom_vline(xintercept=0, color='red')

orgscale.effects <- cbind( 
  effBalpha=invlogit(mean(m[,'mu_alpha'])+m[,'effBalpha'])-invlogit(mean(m[,'mu_alpha'])),
  effBbeta=exp(mean(m[,'mu_beta'])+m[,'effBbeta'])-exp(mean(m[,'mu_beta']))
  )
dotplot(orgscale.effects)+geom_vline(xintercept = 0, color='red')+ggtitle("Effect on original scale at the mean")

A=cbind(alpha=invlogit(m[,'mu_alpha']),
        beta=exp(m[,'mu_beta']))
B=cbind(alpha=invlogit(m[,'mu_alpha']+m[,'effBalpha']),
        beta=exp(m[,'mu_beta']+m[,'effBbeta']))
dotplot(A,B)


alphaA<-(extract(fit, 'alphaA')[[1]])
colnames(alphaA) <- 1:nsubj
alphaB<-(extract(fit, 'alphaB')[[1]])
colnames(alphaB) <- 1:nsubj
dotplot(alphaA, alphaB)




betaA<-(extract(fit, 'betaA')[[1]])
colnames(betaA) <- 1:nsubj
betaB<-(extract(fit, 'betaB')[[1]])
colnames(betaB) <- 1:nsubj
dotplot(betaA, betaB)


## paper plots
alpha<-(extract(fit, 'alphaA')[[1]])
colnames(alpha) <- 1:nsubj
p1 <- cbind(subj=1:nsubj,
      mean=apply(alpha, 2, mean),
      hdi=hdi(alpha)) %>% data.frame %>%
  ggplot(aes(x=subj,y=mean,ymin=lower,ymax=upper))+geom_errorbar(width=0.1)+geom_point()+theme_bw()+
  geom_hline(yintercept = invlogit(mean(m[,'mu_alpha'])))+
  geom_hline(yintercept = invlogit(hdi(m[,'mu_alpha'])[1]), linetype=2)+
  geom_hline(yintercept = invlogit(hdi(m[,'mu_alpha'])[2]), linetype=2)+xlab("subject number")+ylab(expression(alpha))

beta<-(extract(fit, 'betaA')[[1]])
colnames(beta) <- 1:nsubj
p2 <- cbind(subj=1:nsubj,
      mean=apply(beta, 2, mean),
      hdi=hdi(beta)) %>% data.frame %>%
  ggplot(aes(x=subj,y=mean,ymin=lower,ymax=upper))+geom_errorbar(width=0.1)+geom_point()+theme_bw()+
  geom_hline(yintercept = invlogit(mean(m[,'mu_beta'])))+
  geom_hline(yintercept = invlogit(hdi(m[,'mu_beta'])[1]), linetype=2)+
  geom_hline(yintercept = invlogit(hdi(m[,'mu_beta'])[2]), linetype=2)+xlab("subject number")+ylab(expression(beta))

a=multiplot(p1,p2,cols=2)
#ggsave(filename=file.path('graphs/model21.pdf'), plot=a)


##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
## look at final Q-values
## return probability correct under the model for each trial
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
Qvalues <- function( alpha, beta, ACC, pair, reward){
  npairs=max(pair)
  ntrials=length(pair)
  Q=rep(0, 2*npairs)
  for( t in 1:ntrials ) {
    if(ACC[t]<0){ # catch missing values
      pa[t]=NA
      next
    }
    
    # a is selected action in trial t 
    # correct action is always first of pair (A,C,E)
    A=2*pair[t]-1
    B=A+1
    a=A+(ACC[t]<1)
    #cat(sprintf("%i, %i, %i\n", A, B, a))
    Q[a] = Q[a] + alpha*(reward[t]-Q[a])
  }
  Q
}

library(Rcpp)
cppFunction('
NumericVector cQvalues(double alpha, double beta, NumericVector ACC, NumericVector pair, NumericVector reward){
  int npairs=max(pair);
  int ntrials=pair.size();
  
  NumericVector Q=NumericVector(2*npairs, 0.0);
  int A,B,a;
  double pe; // prediction error

  for( int t=0; t<ntrials; t++ ){
    if(ACC[t]<0) // missing responses
      continue;
    
    // a is selected action in trial t
    // correct action is always first of pair (A,C,E)
    A=2*(((int)pair[t])-1);
    B=A+1;
    a=A+(ACC[t]<1);
    
    pe=(reward[t]-Q[a]);
    Q[a] = Q[a] + alpha*pe;
  }
  return(Q);

}')
subj="P01"

with(data, Qvalues(0.16, 0.25, ACCA[subj,], pairA[subj,], rewardA[subj,]))
with(data, cQvalues(0.16, 0.25, ACCA[subj,], pairA[subj,], rewardA[subj,]))

m <- as.matrix(fit)
vars=colnames(m)

## translate alpha/beta and data to final Q-values per sample
rbind(
  do.call(rbind, 
  lapply(1:data$nsubj, 
    function(subj){
      do.call(rbind, 
        lapply(cbind( alpha=m[,sprintf('alphaA[%i]',subj)],beta=m[,sprintf('betaA[%i]', subj)]), function(x){
          with(data, cQvalues(x[1],x[2], ACCA[sprintf("P%02i",subj),], pairA[sprintf("P%02i",subj),], rewardA[sprintf("P%02i",subj),]))
        })) %>% data.frame%>% mutate(subj=subj, condition="A")
        })),
  do.call(rbind, 
          lapply(1:data$nsubj, 
                 function(subj){
                   do.call(rbind, 
                           lapply(cbind( alpha=m[,sprintf('alphaB[%i]',subj)],beta=m[,sprintf('betaB[%i]', subj)]), function(x){
                             with(data, cQvalues(x[1],x[2], ACCB[sprintf("P%02i",subj),], pairB[sprintf("P%02i",subj),], rewardB[sprintf("P%02i",subj),]))
                           })) %>% data.frame%>% mutate(subj=subj, condition="B")
                 }))) -> Q

# Q-plot per subj/condition
Q %>% data.frame %>% gather(Q,value, -subj,-condition) %>%
  ggplot(aes(Q,value,color=condition))+
    stat_summary(fun.data=mean_hdi, geom="pointrange",position=position_dodge(width=0.5))+
    facet_wrap(~subj)

# translate Q to p(A) at the end of session
Q %>% data.frame %>% setNames(c("A","B","C","D","E","F", "subj","condition"))


learn %>% filter(condition=="B", ACC>=0) %>% group_by(subj,pair) %>% mutate(trial=1:n()) %>%
  ggplot(aes(x=trial,y=ACC))+geom_point()+facet_grid(subj~pair)

