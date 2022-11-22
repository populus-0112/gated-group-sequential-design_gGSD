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
#weight_PFS=rbind(c(0.2,0.8,0.2,0.8,0.2,0.8),c(0.3,0.7,0.3,0.7,0.3,0.7))
#weight_OS=rbind(c(0.2,0.8,0.2,0.8,0.2,0.8),c(0.3,0.7,0.3,0.7,0.3,0.7))


##### study design parameter #####
# parameters to generate data
N_all=554 # total No. of subjects recruited
# enrollment rate
durations=c(2,2,2,22)
xx <- tibble(relativeRate=c( 2.5,5,7.5,10 ),duration=durations)
multiplier <- sum(xx$relativeRate*xx$duration)
enroll_rate=N_all/multiplier*xx$relativeRate

dropout_rate_PFS=c(0.0083,0.0083,0.0083,0.0083) # dropout rate(per month): subgroup control, non-subgroup control, subgroup experimental, non-subgroup experimental
dropout_rate_OS=c(0.00083,0.00083,0.00083,0.00083)


tau=c(416/554)       # c(0.4,0.5,0.6,0.7)
threshold_futility=rbind(c(0.9,0.85),c(0.9,0.9))   # c(subgroup, full group)
# get hazard ratio for non-subgroup using hazard ratio of subgroup & full group
hz_PFS=rbind(c(0.7,0.7),c(1,1)) # c(,): hazard ratio for treatment_s, treatment_nonsub
hz_OS=rbind(c(0.7,0.7),c(1,1)) # c(,): hazard ratio for treatment_s, treatment_nonsub
# control group median survival time
OS_ns_median=5.7
OS_S_median=10.5
PFS_ns_median=3
PFS_S_median=4
# log(2)/median survival time
lambda_PFS=rbind(c(log(2)/PFS_S_median,log(2)/PFS_ns_median,log(2)/PFS_S_median*hz_PFS[1,1],log(2)/PFS_ns_median*hz_PFS[1,2]),
                 c(log(2)/PFS_S_median,log(2)/PFS_ns_median,log(2)/PFS_S_median*hz_PFS[2,1],log(2)/PFS_ns_median*hz_PFS[2,2])
) # hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)
lambda_OS=rbind(c(log(2)/OS_S_median,log(2)/OS_ns_median,log(2)/OS_S_median*hz_OS[1,1],log(2)/OS_ns_median*hz_OS[1,2]),
                c(log(2)/OS_S_median,log(2)/OS_ns_median,log(2)/OS_S_median*hz_OS[2,1],log(2)/OS_ns_median*hz_OS[2,2])
) # hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)


# alpha reallocation
# initial alpha for 4 hypotheses: PFS(full group), PFS(subgroup), OS(full group), OS(subgroup).
alpha=c(0.00165,0.00835,0.0022,0.0128)

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
alpha_s_ht=c(0,0.00987,0,0.01513) # alpha for testing subgroup
# initialize alpha for 4 hypotheses
alpha_f_ht=c(0.01071,0,0.01429,0) # alpha for testing full group
transition_ht=rbind(c(0,0,1,0),c(0,0,0,1),c(1,0,0,0),c(0,1,0,0))


# No. of OS events(subgroup)
#n_OS_1s=227  # No. of OS events(subgroup) for IA1
n_OS_2s=314  # No. of OS events(subgroup) for IA2 
n_OS_3s=344  # No. of OS events(subgroup) for FA
#n_OS_1s_ht=227
n_OS_2s_ht=313
n_OS_3s_ht=343
# No. of OS events(full group) 
#n_OS_1=332
n_OS_2=446
n_OS_3=482
#n_OS_1_ht=331
n_OS_2_ht=444
n_OS_3_ht=480
# No. of PFS events(subgroup) 
n_PFS_1s=329
n_PFS_2s=378  
n_PFS_1s_ht=329
n_PFS_2s_ht=379
# No. of PFS events(full group) 
n_PFS_1=447
n_PFS_2=508
n_PFS_1_ht=447
n_PFS_2_ht=508

n_2=c(n_PFS_2,n_PFS_2s,n_OS_2,n_OS_2s) # prespecified No. of events for IA2
n_3=c(0,0,n_OS_3,n_OS_3s) # prespecified No. of events for FA
n_2_ht=c(n_PFS_2_ht,n_PFS_2s_ht,n_OS_2_ht,n_OS_2s_ht) # prespecified No. of events for IA2 of hierarchical test
n_3_ht=c(0,0,n_OS_3_ht,n_OS_3s_ht) # prespecified No. of events for FA of hierarchical test

# futility analysis (PFS)
n_futility=189    #n_PFS_2s*tau[j]  # No. of PFS events for futility


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
write(t(g_1),paste0('case1_total', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_1),sep=',') 
write(t(g_2),paste0('case1_detail', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_2),sep=',')   
write(t(g_3),paste0('case1_stop', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_3),sep=',')   


end_time <- Sys.time()

print(end_time - start_time)