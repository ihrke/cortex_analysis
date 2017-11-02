##
## estimate Q0 from data
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

n.chains=8
n.cores=8
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
#nsubj=data$nsubj
initfct <- function(chainix){
  nsubj=16
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
  fit = stan(file=mod.fname, data = data, chains = n.chains, iter = n.iter, warmup = n.warmup, init = initfct, cores=n.cores)
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
ggsave(filename=file.path('graphs/model21.pdf'), plot=a)

## effects on original and logit scales
effects=m[,str_detect(vars, 'eff')]
cbind(mean=apply(orgscale.effects, 2, mean),
      sd=apply(orgscale.effects, 2, sd),
      hdi=hdi(orgscale.effects))

cbind(mean=apply(effects, 2, mean),
      sd=apply(effects, 2, sd),
      hdi=hdi(effects))


##
group.pars=extract(fit, c("mu_alpha", "mu_beta"))
cbind(mean=unlist(lapply(group.pars, function(x){mean(invlogit(x))})),
      sd=unlist(lapply(group.pars, function(x){sd(invlogit(x))})),
      lower=unlist(lapply(group.pars, function(x){hdi(as.vector(invlogit(x)))[1]})),
      upper=unlist(lapply(group.pars, function(x){hdi(as.vector(invlogit(x)))[2]}))
        )

group.pars=extract(fit, c("sig_alpha", "sig_beta"))
cbind(mean=unlist(lapply(group.pars, function(x){mean(invlogit(x))})),
      sd=unlist(lapply(group.pars, function(x){sd(invlogit(x))})),
      lower=unlist(lapply(group.pars, function(x){hdi(as.vector(invlogit(x)))[1]})),
      upper=unlist(lapply(group.pars, function(x){hdi(as.vector(invlogit(x)))[2]}))
)

## compare Q0 est-model (08) with Q0=0 (02) -> estimates are identical, igore Q0
load("cache/vars/02_stan_qlearning_hier_eff_fit.RData")
fit02 <- fit
load("cache/vars/08_qlearning_hier_eff_Q0est_fit.RData")
fit08 <- fit
library(loo)
compare(loo(extract_log_lik(fit02)), loo(extract_log_lik(fit08)))

m02 <- as.matrix(fit02)
m08 <- as.matrix(fit08)

mean_hdi <- function(x){
  data.frame(y=mean(x),ymin=hdi(x)[1], ymax=hdi(x)[2]) 
}

rbind(
  data.frame(m02) %>% select(starts_with("alphaA")) %>%mutate(model="02"),
  data.frame(m08) %>% select(starts_with("alphaA")) %>%mutate(model="08")) %>%
  gather(var,samples,-model) %>%
  ggplot(aes(x=var,y=samples, color=model))+stat_summary(fun.data=mean_hdi, geom="pointrange", position=position_dodge(width=0.5))

rbind(
  data.frame(m02) %>% select(starts_with("betaA")) %>%mutate(model="02"),
  data.frame(m08) %>% select(starts_with("betaA")) %>%mutate(model="08")) %>%
  gather(var,samples,-model) %>%
  ggplot(aes(x=var,y=samples, color=model))+stat_summary(fun.data=mean_hdi, geom="pointrange", position=position_dodge(width=0.5))

rbind(
  data.frame(m02) %>% select(effBalpha, effBbeta) %>%mutate(model="02"),
  data.frame(m08) %>% select(effBalpha, effBbeta) %>%mutate(model="08")) %>%
  gather(var,samples,-model) %>%
  ggplot(aes(x=var,y=samples, color=model))+stat_summary(fun.data=mean_hdi, geom="pointrange", position=position_dodge(width=0.5))
