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
              double Q[2]={0.0};
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
  res=mclapply(data.frame(t(pars)), function(x){ mean(replicate(nrep, csimandeval(x[1], x[2], ntrials, x[3]))) }, mc.cores=ncpu)
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


x=with(result, beta[alpha==0.08])
y=with(result, result[alpha==0.08])
f=smooth.spline(x,y)
predict(f)[['x']]
