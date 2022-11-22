source("function.R")
source("iteration.R")

library(simtrial) # generate data
library(gsDesign) # IA and FA boundary
library(survival) # Cox model for hazard ratio
#library(gMCP) # simes.test

# read seed number
seed_0=read.table('seed_1.csv',sep=',')
seed_0=as.matrix(seed_0)
iter=2000              


# weights for p-value combination test
weight_PFS=rbind(c(0.5,0.5,0.5,0.5,0.5,0.5),c(0.5,0.5,0.5,0.5,0.5,0.5),c(0.7,0.3,0.7,0.3,0.7,0.3),c(0,0,0,0,0,0)) # IA1 w1, IA1 w2, IA2 w1, IA2 w2, FA w1, FA w2
weight_OS=rbind(c(0.5,0.5,0.5,0.5,0.5,0.5),c(0.7,0.3,0.7,0.3,0.7,0.3),c(0.7,0.3,0.7,0.3,0.7,0.3),c(0,0,0,0,0,0))  # IA1 w1, IA1 w2, IA2 w1, IA2 w2, FA w1, FA w2


##### study design parameter #####
# parameters to generate data
N_all=924 # total No. of subjects recruited
# enrollment rate
durations=c(2,2,2,27)
xx <- tibble(relativeRate=c( 2.5,5,7.5,10 ),duration=durations)
multiplier <- sum(xx$relativeRate*xx$duration)
enroll_rate=N_all/multiplier*xx$relativeRate

dropout_rate_PFS=rep(0.01,4) # dropout rate(per month): subgroup control, non-subgroup control, subgroup experimental, non-subgroup experimental
dropout_rate_OS=rep(0.001,4)


tau=c(462/924)       # c(0.4,0.5,0.6,0.7)
threshold_futility=rbind(c(0.85,0.83),c(0.9,0.9),c(0.85,0.83))   # c(subgroup, full group)
# get hazard ratio for non-subgroup using hazard ratio of subgroup & full group
hz_PFS=rbind(c(0.7,0.7),c(1,1),c(0.7,3/2.28)) # c(,): hazard ratio for treatment_s, treatment_nonsub
hz_OS=rbind(c(0.7,0.7),c(1,1),c(0.7,5.7/4.37)) # c(,): hazard ratio for treatment_s, treatment_nonsub
# control group median survival time
OS_ns_median=5.7
OS_S_median=10.5
PFS_ns_median=3
PFS_S_median=4
# log(2)/median survival time
lambda_PFS=rbind(c(log(2)/PFS_S_median,log(2)/PFS_ns_median,log(2)/PFS_S_median*hz_PFS[1,1],log(2)/PFS_ns_median*hz_PFS[1,2]),
                 c(log(2)/PFS_S_median,log(2)/PFS_ns_median,log(2)/PFS_S_median*hz_PFS[2,1],log(2)/PFS_ns_median*hz_PFS[2,2]),
                 c(log(2)/PFS_S_median,log(2)/PFS_ns_median,log(2)/PFS_S_median*hz_PFS[3,1],log(2)/PFS_ns_median*hz_PFS[3,2])
) # hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)
lambda_OS=rbind(c(log(2)/OS_S_median,log(2)/OS_ns_median,log(2)/OS_S_median*hz_OS[1,1],log(2)/OS_ns_median*hz_OS[1,2]),
                c(log(2)/OS_S_median,log(2)/OS_ns_median,log(2)/OS_S_median*hz_OS[2,1],log(2)/OS_ns_median*hz_OS[2,2]),
                c(log(2)/OS_S_median,log(2)/OS_ns_median,log(2)/OS_S_median*hz_OS[3,1],log(2)/OS_ns_median*hz_OS[3,2])
) # hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)


# alpha reallocation
# initial alpha for 4 hypotheses: PFS(full group), PFS(subgroup), OS(full group), OS(subgroup).
alpha=c(0.00017,0.01,0.00025,0.01458)

transition=rbind(c(0,0,1,0),c(1,0,0,0),c(1,0,0,0),c(0,0,1,0))
# transition is the matrix with weights of alpha reallocation when any hypothesis is significant
# transition[i,j]: the weights of alpha reallocation from hypothesis i to hypothesis j if hypothesis i is significant 
# adaptive design
alpha_f=c(alpha[1],0,alpha[3],0)  # only full group is selected in adaptive design
alpha_s=c(0,alpha[2],0,alpha[4])  # only subgroup is selected in adaptive design
transition_f=rbind(c(0,0,1,0),c(0,0,0,0),c(1,0,0,0),c(0,0,0,0)) # only full group is selected in adaptive design
transition_s=rbind(c(0,0,0,0),c(0,0,0,0),c(0,0,0,0),c(0,0,0,0)) # only subgroup is selected in adaptive design
### hierarchical test ##
# initialize alpha for 4 hypotheses
alpha_s_ht=c(0,0.0102,0,0.0148) # alpha for testing subgroup
# initialize alpha for 4 hypotheses
alpha_f_ht=c(0.01012,0,0.01488,0) # alpha for testing full group
transition_ht=rbind(c(0,0,1,0),c(0,0,0,1),c(1,0,0,0),c(0,1,0,0))


# No. of OS events(subgroup)
#n_OS_1s=211  # No. of OS events(subgroup) for IA1
n_OS_2s=298  # No. of OS events(subgroup) for IA2 
n_OS_3s=348  # No. of OS events(subgroup) for FA
#n_OS_1s_ht=211
n_OS_2s_ht=297
n_OS_3s_ht=347
# No. of OS events(full group) 
#n_OS_1=502
n_OS_2=688
n_OS_3=779
#n_OS_1_ht=502
n_OS_2_ht=687
n_OS_3_ht=778
# No. of PFS events(subgroup) 
n_PFS_1s=312
n_PFS_2s=395
n_PFS_1s_ht=311
n_PFS_2s_ht=394
# No. of PFS events(full group) 
n_PFS_1=647
n_PFS_2=810
n_PFS_1_ht=647
n_PFS_2_ht=810

n_2=c(n_PFS_2,n_PFS_2s,n_OS_2,n_OS_2s) # prespecified No. of events for IA2
n_3=c(0,0,n_OS_3,n_OS_3s) # prespecified No. of events for FA
n_2_ht=c(n_PFS_2_ht,n_PFS_2s_ht,n_OS_2_ht,n_OS_2s_ht) # prespecified No. of events for IA2 of hierarchical test
n_3_ht=c(0,0,n_OS_3_ht,n_OS_3s_ht) # prespecified No. of events for FA of hierarchical test

# futility analysis (PFS)
n_futility=197    #n_PFS_2*tau[j]  # No. of PFS events for futility


############# three methods ############# 
start_time <- Sys.time()

out=sf_iteration(iter,seed_0,N_all,tau,lambda_PFS,lambda_OS,dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations,
                 weight_PFS,weight_OS,n_futility,threshold_futility,hz_PFS,hz_OS,alpha,transition,alpha_f,alpha_s,
                 transition_f,transition_s,alpha_s_ht,alpha_f_ht,transition_ht,n_PFS_1s,n_PFS_2s,n_OS_3s,n_PFS_1,
                 n_PFS_2,n_OS_3,n_2,n_3,n_PFS_1s_ht,n_PFS_2s_ht,n_OS_3s_ht,n_PFS_1_ht,n_PFS_2_ht,n_OS_3_ht,
                 n_2_ht,n_3_ht)

g_1=out$total
g_2=out$detail
g_3=out$stop
write(t(g_1),paste0('case23_total', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_1),sep=',') 
write(t(g_2),paste0('case23_detail', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_2),sep=',')  
write(t(g_3),paste0('case23_stop', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_3),sep=',')   


end_time <- Sys.time()

print(end_time - start_time)