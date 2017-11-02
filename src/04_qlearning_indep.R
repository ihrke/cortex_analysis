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
n.iter=200
n.warmup=100

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
traceplot(fit, pars=c("betaA"))

plot(fit, pars=vars[str_detect(vars,"eff")])

#waic(fit)


##

learn %>% filter(subj=="P08") %>% group_by(condition, pair) %>%
  mutate(trial=1:n()) %>%  mutate(prev_reward=c(NA, reward[-n()])) %>% filter(ACC>=0) %>%
  ggplot(aes(x=trial, y=jitter(ACC,amount=0.05), color=side))+geom_point(size=2)+facet_grid(condition~pair)



learn %>%  group_by(subj,condition, pair) %>%
  mutate(trial=1:n()) %>%  filter(ACC>=0) %>%
  ggplot(aes(x=trial, y=jitter(ACC,amount=0.05)))+geom_point(size=2)+facet_grid(pair~subj)

