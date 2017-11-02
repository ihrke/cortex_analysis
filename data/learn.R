fpath=file.path('data', 'raw', 'Learning')

fnames <- c(list.files( file.path(fpath, "A"), full.names=T ), 
            list.files( file.path(fpath, "B"), full.names=T))
                      

learn <- NULL
for( fname in fnames ){
  condition = strsplit(fname, "/")[[1]][4]
  subj=strsplit(strsplit(fname, "/")[[1]][5], "_")[[1]][1]
  d <- read.table(fname, sep="\t",  header=T, row.names=NULL,
                  col.names=c("block", "pair", "reward", "side", "RT", "ACC"))
  d <- d[,1:6]
  names(d) <- c("block", "pair", "reward", "side", "RT", "ACC")
  d <- within(d, {
    ACC <- as.character(ACC); 
    ACC[ACC=="missed"]<--1; 
    ACC<-as.integer(ACC);
    RT[RT==0] <- NA;
    }) %>% 
    droplevels %>% mutate(subj=subj, condition=condition)
  
  #print(str(d))
  learn <- rbind(learn, d)
}

learn <- within(learn,{
                condition=as.factor(condition);
                subj=as.factor(subj);
                block=as.integer(block);
                })