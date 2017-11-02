library(ProjectTemplate)
load.project()
theme_set(theme_bw())

# enabls --force for stallo etc.
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options)
  uncache.all()

ntrials=120


library(Rcpp)
cppFunction('double csimandeval(double alpha, double beta, int ntrials, double rewardprobA){
              double ACC=0.0;
              int t;
              double Q[2]={0.5,0.5};
              int A,B,a,reward;
              double pA, offset;
            
              for(t=0; t<ntrials; t++){
                A=0;
                B=1;
                offset=((Q[A]/beta)>(Q[B]/beta) ? (Q[A]/beta) : (Q[B]/beta));
                pA=exp(Q[A]/beta - (log( exp( Q[A]/beta - offset)+exp( Q[B]/beta - offset)) + offset));
                a=(((double)runif(1)[0])<pA) ? A : B;
                reward=(((double)runif(1)[0])<rewardprobA) ? (a==A) : (a==B);
                //Rprintf("pA=%f, a=%i, reward=%i, ACC=%i\\n", pA, a, reward, a==A);
              
                Q[a] = Q[a] + alpha*(reward-Q[a]);
                ACC = ACC+(a==A);
              }
              return ACC/(ntrials);
            }
            ')


rewardprobs=c(0.6, 0.7, 0.8)
alphas=seq(0.01, .5,by=0.005)
betas=seq(0.01,  .5,by=0.005)


## parallel?
library(parallel)
ncpu=8
nrep=10000
pars<-expand.grid(alphas, betas, rewardprobs) %>% setNames(c("alpha", "beta", "prewardA"))

if(is.cached.var('result')){
  result=load.cache.var('result')
} else {
  res=mclapply(data.frame(t(pars)), function(x){ median(replicate(nrep, csimandeval(x[1], x[2], ntrials, x[3]))) }, mc.cores=ncpu)
  result=cbind(pars, result=unlist(res))
  cache.var('result')  
}

result %>% group_by(prewardA) %>%
  #mutate(result=result/max(result)) %>%
  ggplot(aes(alpha, beta, z = result))+geom_raster(aes(fill=result))+
    facet_wrap(~prewardA, scales="free")+
    scale_fill_gradientn(colours = jet.colors(7))

alpha.slice=c(0.03, 0.04, 0.05, 0.08, 0.09, 0.10, 0.12, 0.15, 0.20)
x=alpha.slice[1]

do.call(rbind, 
        lapply(alpha.slice, function(x){
          result %>% filter(as.factor(alpha)==x) %>% mutate(prewardA=as.factor(prewardA)) %>% mutate(alpha.slice=x)
        })) %>%
  ggplot(aes(beta, result,color=prewardA,group=prewardA))+geom_line()+facet_wrap(~alpha.slice)
  

## smoothing spline
do.call(rbind, 
        lapply(alpha.slice, function(x){
          result %>% filter(as.factor(alpha)==x) %>% mutate(prewardA=as.factor(prewardA)) %>% mutate(alpha.slice=x) %>%
            group_by(prewardA) %>%
            mutate(smooth=predict(smooth.spline(beta,result))[['y']]) %>% ungroup
        })) %>%
  ggplot(aes(beta, result,color=prewardA,group=prewardA))+
    geom_point()+geom_line(aes(y=smooth))+
    geom_point(data=(.%>% group_by(alpha.slice,prewardA) %>% 
                       summarise(ymax=max(smooth),xmax=beta[which.max(smooth)])), 
               aes(x=xmax,y=ymax,color=prewardA), shape=9,size=5)+
  
    geom_segment(data=(.%>% group_by(alpha.slice,prewardA) %>% 
                       summarise(ymax=max(smooth),xmax=beta[which.max(smooth)])), 
               aes(x=xmax,y=ymax, xend=xmax,yend=0.5,color=prewardA))+
  
    facet_wrap(~alpha.slice)


## distribution of ACC?
nrep=20000
alpha=0.10
beta=0.05
do.call(rbind,lapply(rewardprobs, function(x){
  data.frame(pA=x,samples=replicate(nrep, csimandeval(alpha, beta, ntrials, x)))})) -> samp


## method of moments for beta function (wikipedia)
mm.beta <- function(x){
  c(shape1=mean(x)*( (mean(x)*(1-mean(x)))/var(x) -1 ), 
    shape2=(1-mean(x))*( (mean(x)*(1-mean(x))/var(x) -1)))
}

samp %>%
  ggplot(aes(x=samples, fill=as.factor(pA)))+
    geom_histogram(aes(y=..density..), position = "identity", binwidth=0.01, alpha=0.4)+
    ggtitle(sprintf("a=%.2f, b=%.2f", alpha, beta))+xlim(0,1)#+
      #stat_function(fun=dbeta, args=with(samp, mm.beta(samples[pA==0.6])))+
      #stat_function(fun=dbeta, args=with(samp, mm.beta(samples[pA==0.7])))+
      #stat_function(fun=dbeta, args=with(samp, mm.beta(samples[pA==0.8])))

lapply(rewardprobs, function(x){with(samp, mm.beta(samples[pA==x]))}) %>% setNames(rewardprobs)


## 
pA.fun <- function(alpha, beta, ntrials, rewardprob){
  Q=c(0.5,0.5)
  A=1;
  B=2;
  d=data.frame(trial=1:ntrials,
               alpha=alpha,
               beta=beta,
               rewardprob=rewardprob,
               ACC=-1,
               reward=-1,
               pA=NA)
  for(t in 1:ntrials){
    offset=ifelse( (Q[A]/beta)>(Q[B]/beta), (Q[A]/beta), (Q[B]/beta) )
    pA=exp(Q[A]/beta - (log( exp( Q[A]/beta - offset)+exp( Q[B]/beta - offset)) + offset))
    a=ifelse(runif(1)<pA, A, B) 
    reward=as.integer(ifelse(runif(1)<rewardprob , (a==A), (a==B) ))
    Q[a] = Q[a] + alpha*(reward-Q[a])
    d$ACC[t]=as.integer(a==A)
    d$reward[t]=reward
    d$pA[t]=pA
  }
  d %>% 
    ggplot(aes(x=trial,y=pA))+geom_line()+geom_point(aes(y=ACC, color=as.factor(reward)))+ylim(0,1)
}

}

