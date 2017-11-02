fpath=file.path('data', 'raw', 'Transfer')

fnames <- c(list.files( file.path(fpath, "A"), full.names=T ), 
            list.files( file.path(fpath, "B"), full.names=T))

transfer <- NULL
for( fname in fnames ){
  condition = strsplit(fname, "/")[[1]][4]
  subj=strsplit(strsplit(fname, "/")[[1]][5], "_")[[1]][1]
  d <- read.table(fname, sep="\t",  header=T, row.names=NULL)
  names(d) <- c("side", "RT", "ACC", "type", 'tmp')
  d <- within(d, {
    ACC <- as.character(ACC); 
    ACC[ACC=="missed"]<- -1; 
    ACC<-as.integer(ACC);
    RT[RT==0] <- NA;
  }) %>% 
    droplevels %>% mutate(subj=subj, condition=condition) %>% select(-tmp)
  
  transfer <- rbind(transfer, d)
}

# 1:AB
# 2:AC
# 3:AD
# 4:AE
# 5:AF
# 6:BC
# 7:BD
# 8:BE
# 9:BF
# 10:CD
# 11:CE
# 12:CF
# 13:DE
# 14:DF
# 15:EF

transfer <- within(transfer, 
       type<-factor(type, levels=1:15, 
                    labels=c("AB","AC", "AD", "AE", "AF", 
                             "BC", "BD", "BE", "BF", "CD", "CE", "CF", 
                             "DE", "DF", "EF" )))

transfer <- within(transfer,{
  condition=as.factor(condition);
  subj=as.factor(subj);
  side=as.factor(side);
})

## chosen object (ACC variable codes whether the object with higher reward probability would be chosen)
# A=0.8
# B=0.2
# C=0.7
# D=0.3
# E=0.6
# F=0.4
within(transfer,{
  chosen <- rep("miss", length(subj));
  chosen[type=="AB" & ACC==1] <- "A";
  chosen[type=="AB" & ACC==0] <- "B";

  chosen[type=="AC" & ACC==1] <- "A";
  chosen[type=="AC" & ACC==0] <- "C";

  chosen[type=="AD" & ACC==1] <- "A";
  chosen[type=="AD" & ACC==0] <- "D";

  chosen[type=="AE" & ACC==1] <- "A";
  chosen[type=="AE" & ACC==0] <- "E";

  chosen[type=="AF" & ACC==1] <- "A";
  chosen[type=="AF" & ACC==0] <- "F";
  #----------------------------------
  chosen[type=="BC" & ACC==1] <- "C";
  chosen[type=="BC" & ACC==0] <- "B";
  
  chosen[type=="BD" & ACC==1] <- "D";
  chosen[type=="BD" & ACC==0] <- "B";
  
  chosen[type=="BE" & ACC==1] <- "E";
  chosen[type=="BE" & ACC==0] <- "B";
  
  chosen[type=="BF" & ACC==1] <- "F";
  chosen[type=="BF" & ACC==0] <- "B";
  #----------------------------------
  chosen[type=="CD" & ACC==1] <- "C";
  chosen[type=="CD" & ACC==0] <- "D";
  
  chosen[type=="CE" & ACC==1] <- "C";
  chosen[type=="CE" & ACC==0] <- "E";

  chosen[type=="CF" & ACC==1] <- "C";
  chosen[type=="CF" & ACC==0] <- "F";
  #----------------------------------
  chosen[type=="DE" & ACC==1] <- "E";
  chosen[type=="DE" & ACC==0] <- "D";
  
  chosen[type=="DF" & ACC==1] <- "F";
  chosen[type=="DF" & ACC==0] <- "D";
  #----------------------------------
  chosen[type=="EF" & ACC==1] <- "E";
  chosen[type=="EF" & ACC==0] <- "F";
}) -> transfer
