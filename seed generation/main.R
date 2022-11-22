args <- commandArgs(TRUE)
start_0 <- strtoi(args[1], base=10L)
#end_0 <- strtoi(args[2], base=10L)
print(start_0) 

####### generate seed ####### 
iter=1000 # simulation times

# generate random number for seed
seed_0=rep(0,iter) # different parameter set up use the same set of seed
for (i in 1:iter){
  seed_0[i]=sample.int(iter*1000,1)
}
write(t(seed_0),'seed_0.csv',ncol=1,sep=',')





####### delete seeds that failed to generate gs boundary because of No. of events ####### 
iter=1000            # choose seeds that works for case 1, case 2 and case 3 among the first 'iter' seeds in "seed_0.csv"
# increase this number if simulation times becomes bigger


source("case1.R")  # choose seeds that works for case 1 among 'iter' seeds

rm(list=ls(all=TRUE)) # remove all values

iter=1000

source("case23.R")  # choose seeds that works for case 2 and case 3 among the first 'iter' seeds in "seed_0.csv"

rm(list=ls(all=TRUE)) # remove all values

a1=read.table('case1_seed_fail.csv',sep=',')
a1=as.matrix(a1)

#a2=read.table('case23_seed_fail.csv',sep=',')
#a2=as.matrix(a2)
#a=c(a1,a2)

a22 <- try(read.table('case23_seed_fail.csv',sep=','))
if (class(a22)=='try-error'){
  a=a1
} else {
  a2=read.table('case23_seed_fail.csv',sep=',')
  a2=as.matrix(a2)
  a=c(a1,a2)
}

    
b=sort(unique(a))  # choose seeds that works for case 1, case 2 and case 3 among the first 'iter' seeds in "seed_0.csv"

seed_0=read.table('seed_0.csv',sep=',')
seed_0=as.matrix(seed_0)

for (i in length(b):1){
  seed_0=seed_0[-b[i]]
}

#write(t(seed_0),paste0('seed_1_new', ".csv"),ncol=1,sep=',')
start_0=sample.int(100000,1)
write(t(seed_0),paste0('seed_1_new_',start_0,'.csv'),ncol=1,sep=',')

