library(ProjectTemplate)
load.project()

library(lmerTest)
theme_set(theme_bw())

## barplot/SEM ACC/RT
learn %>% filter(ACC>=0) %>% mutate(condition=fct_recode(condition,anodal="A",
                                                         sham="B")) %>% 
  group_by(subj,condition,pair) %>%
  summarise(accuracy=mean(ACC), RT=mean(RT)) %>% gather(var, val, accuracy, RT) %>%
  ggplot(aes(x=pair, y=val, fill=condition))+
  stat_summary(fun.y=mean, geom="bar", position="dodge")+
  stat_summary(fun.data=mean_se, geom="pointrange",  position=position_dodge(width=1))+
  facet_wrap(~var, scales="free")

#---------------------------------------
# Learning: ACC
#---------------------------------------
library(forcats)
## barplot/SEM
learn %>% filter(ACC>=0) %>% mutate(condition=fct_recode(condition,anodal="A",
                                                         sham="B")) %>% 
  group_by(subj,condition,pair) %>%
  summarise(accuracy=mean(ACC)) %>% 
  ggplot(aes(x=pair, y=accuracy, fill=condition))+
    stat_summary(fun.y=mean, geom="bar", position="dodge")+
    stat_summary(fun.data=mean_se, geom="pointrange",  position=position_dodge(width=1))

  
## difference-barplot
learn %>% filter(ACC>=0) %>% group_by(subj,condition,pair) %>%
  summarise(ACC=mean(ACC)) %>% spread(condition, ACC) %>%
  mutate(dACC=A-B) %>%
  ggplot(aes(x=pair, y=dACC))+
    stat_summary(fun.y=mean, geom="bar", position="dodge")+
    stat_summary(fun.data=mean_cl_boot, geom="errorbar",  position="dodge", width=0.1)+
    ylim(-.2, .2)+ggtitle("anodal-sham")

  
## logistic regression for ACC
learn %>% group_by(subj, condition, pair) %>% mutate(trial=1:n()) %>% ungroup %>%
  filter(ACC>=0) %>% data.frame -> d

summary(mod <- glmer( ACC ~ scale(trial) + as.numeric(pair) + condition + (1|subj), 
                      data=d, family=binomial))
summary(mod)

summary(mod2 <- glmer( ACC ~ scale(trial) + as.factor(pair) * condition + (1|subj), 
                      data=d, family=binomial))


# illustratory plot
nsamp=2 # sample nsamp random subjs
use.subjs<-as.vector(sample(unique(d$subj), nsamp))

d %>% mutate(pred=arm::invlogit(predict(mod))) %>% filter(subj %in% use.subjs) %>%
  ggplot(aes(x=trial,y=pred,colour=condition, group=condition))+geom_line()+
  geom_point(aes(y=ACC+0.05*as.numeric(condition)), size=2)+
  facet_grid(subj ~ pair) 

#---------------------------------------
# Learning: RT
#---------------------------------------

## regression for RT/log(RT)
d<-learn %>% group_by(subj, condition) %>% mutate(exptrial=1:n()) %>% ungroup %>% 
  group_by(subj, condition, pair) %>% mutate(trial=1:n()) %>% ungroup %>%
  filter(ACC>=0)


summary(mod <- lmer((RT) ~ scale(trial) + as.numeric(pair) + condition+(1|subj), data=d))
summary(mod2 <- lmer(log(RT) ~ scale(trial) + as.numeric(pair) + condition+(1|subj), data=d))

summary(mod <- lmer((RT) ~ scale(trial) + as.factor(pair) * condition+(1|subj), data=d))
summary(mod2 <- lmer(log(RT) ~ scale(trial) + as.factor(pair) * condition+(1|subj), data=d))

## barplot/SEM
learn %>% filter(ACC>=0) %>% group_by(subj,condition,pair) %>%
  summarise(RT=mean(RT)) %>% 
  ggplot(aes(x=pair, y=RT, fill=condition))+
  stat_summary(fun.y=mean, geom="bar", position="dodge")+
  stat_summary(fun.data=mean_se, geom="errorbar",  position="dodge")

## difference-barplot
learn %>% filter(ACC>=0) %>% group_by(subj,condition,pair) %>%
  summarise(RT=mean(RT)) %>% spread(condition, RT) %>%
  mutate(dRT=A-B) %>%
  ggplot(aes(x=pair, y=dRT))+
  stat_summary(fun.y=mean, geom="bar", position="dodge")+
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",  position="dodge", width=0.1)+ggtitle("anodal-sham")


d %>% ggplot(aes(x=(RT), fill=condition))+geom_histogram(position="identity", alpha=.4)+facet_wrap(~subj)

#---------------------------------------
## Learning: p(stay|win), p(stay|lose): barplot
#---------------------------------------

learn %>% filter(ACC>=0) %>% group_by(subj, condition, pair) %>%
  mutate(prev_reward=c(NA, reward[-n()]), stay=c(NA, diff(ACC)==0)) %>%
  summarise(pstaywin=sum( stay==1 & prev_reward==1, na.rm=T)/(sum(prev_reward==1, na.rm=T)),
            pstaylose=sum( stay==1 & prev_reward==0, na.rm=T)/(sum(prev_reward==0, na.rm=T))) %>%
  gather(winlose, pstay, -subj, -condition, -pair) %>%
  ggplot(aes(x=pair, y=pstay, fill=condition))+
    stat_summary(fun.y=mean, geom="bar", position="dodge")+
    stat_summary(fun.data=mean_se, geom="errorbar",  position="dodge")+
    facet_wrap(~winlose)

## p(stay|win), p(stay|lose): barplot of condition difference
learn %>% filter(ACC>=0) %>% group_by(subj, condition, pair) %>%
  mutate(prev_reward=c(NA, reward[-n()]), stay=c(NA, diff(ACC)==0)) %>%
  summarise(pstaywin=sum( stay==1 & prev_reward==1, na.rm=T)/(sum(prev_reward==1, na.rm=T)),
            pstaylose=sum( stay==1 & prev_reward==0, na.rm=T)/(sum(prev_reward==0, na.rm=T))) %>%
  gather(winlose, pstay, -subj, -condition, -pair) %>%
  spread(condition, pstay) %>%
  mutate(diffpstay=A-B) %>%
    ggplot(aes(x=pair, y=diffpstay))+
      stat_summary(fun.y=mean, geom="bar", position="dodge")+
      stat_summary(fun.data=mean_se, geom="errorbar",  position="dodge")+
      facet_wrap(~winlose)


## pstay: everything together
learn %>% filter(ACC>=0) %>% group_by(subj, condition, pair) %>%
  mutate(prev_reward=c(NA, reward[-n()]), stay=c(NA, diff(ACC)==0)) %>%
  summarise(pstaywin=sum( stay==1 & prev_reward==1, na.rm=T)/(sum(prev_reward==1, na.rm=T)),
            pstaylose=sum( stay==1 & prev_reward==0, na.rm=T)/(sum(prev_reward==0, na.rm=T))) %>% data.frame -> tmp


summary(mod <- lmer( pstaywin ~ condition * as.factor(pair) + (1|subj), 
                      data=tmp))

summary(mod <- lmer( pstaylose ~ condition * as.factor(pair) + (1|subj), 
                     data=tmp))

tmp %>%
  gather(winlose, pstay, -subj, -condition, -pair) -> tmp2
  
summary(mod <- lmer( pstay ~ winlose * condition * as.factor(pair) + (1|subj), 
                     data=tmp2))

#---------------------------------------
## Transfer phase
#---------------------------------------

## choose A
transfer %>% filter(str_detect(type, "A")) %>%
  group_by(subj, condition) %>% summarize(chooseA=sum(chosen=="A")/n()) %>%
  ggplot(aes(x=condition, y=chooseA))+
    stat_summary(fun.y=mean, geom="bar")+
    stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0.1)

## avoid B
transfer %>% filter(str_detect(type, "B")) %>%
  group_by(subj, condition) %>% summarize(avoidB=sum(chosen!="B")/n()) %>%
  ggplot(aes(x=condition, y=avoidB))+
  stat_summary(fun.y=mean, geom="bar")+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", width=0.1)


## only for A
# objects=c("A", "B", "C", "D", "E", "F")
# x="A"
# transfer %>% filter(str_detect(type, x)) %>% 
#   mutate(comp=factor(str_replace(type, x, ""), levels=objects)) %>% 
#   group_by(subj,condition,comp) %>%
#   summarise(choosex=sum(chosen==x)/n()) %>% 
#   ggplot(aes(x=comp, y=choosex, fill=condition))+
#     stat_summary(fun.y=mean, geom="bar", position="dodge")+
#     stat_summary(fun.data=mean_cl_boot, geom="errorbar", position="dodge")+
#     scale_x_discrete(drop=F)


## Barplot: Proportion of choosing X in each combination separately by condition
do.call(rbind, 
  lapply(objects, function(x){
    transfer %>% filter(str_detect(type, x)) %>% 
      mutate(comp=factor(str_replace(type, x, ""), levels=objects)) %>% 
      group_by(subj,condition,comp) %>%
      summarise(choosex=sum(chosen==x)/n()) %>% 
      mutate(baseline=x) %>% data.frame
  })
)%>% ggplot(aes(x=comp, y=choosex, fill=condition))+
  stat_summary(fun.y=mean, geom="bar", position="dodge")+
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", position="dodge")+
  scale_x_discrete(drop=F)+facet_grid(baseline~.)#+coord_flip()

## same but plot mean of difference sham-anodal
do.call(rbind, 
        lapply(objects, function(x){
          transfer %>% filter(str_detect(type, x)) %>% 
            mutate(comp=factor(str_replace(type, x, ""), levels=objects)) %>% 
            group_by(subj,condition,comp) %>%
            summarise(choosex=sum(chosen==x)/n()) %>% 
            mutate(baseline=x) %>% data.frame
        })
) %>% spread( condition, choosex) %>%
  mutate(sham_anodal=B-A) %>%
  ggplot(aes(x=comp, y=sham_anodal))+
    stat_summary(fun.y=mean, geom="bar", position="dodge")+
    stat_summary(fun.data=mean_cl_boot, geom="errorbar", position="dodge")+
    scale_x_discrete(drop=F)+facet_grid(baseline~.)+coord_flip()


### individual data plot
learn %>% group_by(subj, pair, condition) %>%
  mutate(trial=1:n()) %>% ungroup %>%
  filter(ACC>=0) %>%
  group_by(trial, pair, condition) %>%
  mutate(rep.ix=(1:n())/n()) %>% ungroup %>%
  ggplot(aes(x=trial, y=ACC))+
  geom_rect(aes(xmin=trial-.5, xmax=trial+.5, ymin=rep.ix, ymax=rep.ix+0.06, fill=factor(ACC)))+
  stat_summary(fun.data = mean_cl_boot, geom="line", color="red")+
  facet_grid(condition~pair)+scale_fill_manual(values=c("grey","white"))+theme_bw()
  

