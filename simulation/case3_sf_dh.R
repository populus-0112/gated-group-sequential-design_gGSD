source("function.R")
source("iteration_dh.R")

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
threshold_futility=rbind(c(0.85,0.83))   # c(subgroup, full group)
# get hazard ratio for non-subgroup using hazard ratio of subgroup & full group
hz_PFS=rbind(c(0.7,3/2.28)) # c(,): hazard ratio for treatment_s, treatment_nonsub
hz_OS=rbind(c(0.7,5.7/4.37)) # c(,): hazard ratio for treatment_s, treatment_nonsub
# control group median survival time
OS_ns_median=5.7
OS_S_median=10.5
PFS_ns_median=3
PFS_S_median=4
# log(2)/median survival time
lambda_PFS=rbind(c(log(2)/PFS_S_median,log(2)/PFS_ns_median,log(2)/PFS_S_median*hz_PFS[1,1],log(2)/PFS_ns_median*hz_PFS[1,2]))
# hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)
lambda_OS=rbind(c(log(2)/OS_S_median,log(2)/OS_ns_median,log(2)/OS_S_median*hz_OS[1,1],log(2)/OS_ns_median*hz_OS[1,2]))
# hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)


### hierarchical test ##
# initialize alpha for 4 hypotheses
alpha_s_ht=c(0,0.0102,0,0.0148) # alpha for testing subgroup
alpha_s_ht_dh_1=c(0,0.025,0,0) # alpha for testing subgroup when only subgroup selected in futility analysis
alpha_s_ht_dh_2=c(0,0,0,0.025) # alpha for testing subgroup when only subgroup selected in futility analysis
# initialize alpha for 4 hypotheses
alpha_f_ht=c(0.01012,0,0.01488,0) # alpha for testing full group
transition_ht=rbind(c(0,0,1,0),c(0,0,0,1),c(1,0,0,0),c(0,1,0,0))
# transition_ht is the matrix with weights of alpha reallocation when any hypothesis is significant
# transition_ht[i,j]: the weights of alpha reallocation from hypothesis i to hypothesis j if hypothesis i is significant 



# No. of OS events(subgroup)
#n_OS_1s_ht=211
n_OS_2s_ht=297
n_OS_3s_ht=347
n_OS_2s_ht_dh=297
n_OS_3s_ht_dh=347
# No. of OS events(full group) 
#n_OS_1_ht=502
n_OS_2_ht=687
n_OS_3_ht=778
n_OS_2_ht_dh=687
n_OS_3_ht_dh=778
# No. of PFS events(subgroup) 
n_PFS_1s_ht=311
n_PFS_2s_ht=394
n_PFS_1s_ht_dh=311
n_PFS_2s_ht_dh=394
# No. of PFS events(full group) 
n_PFS_1_ht=647
n_PFS_2_ht=810
n_PFS_1_ht_dh=647
n_PFS_2_ht_dh=810

n_2_ht=c(n_PFS_2_ht,n_PFS_2s_ht,n_OS_2_ht,n_OS_2s_ht) # prespecified No. of events for IA2 of hierarchical test
n_3_ht=c(0,0,n_OS_3_ht,n_OS_3s_ht) # prespecified No. of events for FA of hierarchical test
n_2_ht_dh=c(n_PFS_2_ht_dh,n_PFS_2s_ht_dh,n_OS_2_ht_dh,n_OS_2s_ht_dh) # prespecified No. of events for IA2 of hierarchical test
n_3_ht_dh=c(0,0,n_OS_3_ht_dh,n_OS_3s_ht_dh) # prespecified No. of events for FA of hierarchical test

# futility analysis (PFS)
n_futility=197    #n_PFS_2*tau[j]  # No. of PFS events for futility


############# three methods ############# 
start_time <- Sys.time()

out=sf_iteration_dh(iter,seed_0,N_all,tau,lambda_PFS,lambda_OS,dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations,
                    weight_PFS,weight_OS,n_futility,threshold_futility,hz_PFS,hz_OS,alpha_s_ht_dh_1,alpha_s_ht_dh_2,
                    alpha_s_ht,alpha_f_ht,transition_ht,n_PFS_1s_ht,n_PFS_2s_ht,n_OS_3s_ht,n_PFS_1_ht,n_PFS_2_ht,n_OS_3_ht,
                    n_2_ht,n_3_ht,n_PFS_1s_ht_dh,n_PFS_2s_ht_dh,n_OS_3s_ht_dh,n_PFS_1_ht_dh,n_PFS_2_ht_dh,n_OS_3_ht_dh,
                    n_2_ht_dh,n_3_ht_dh)

g_1=out$total
g_2=out$detail
g_3=out$stop
write(t(g_1),paste0('case3(dh)_total', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_1),sep=',') 
write(t(g_2),paste0('case3(dh)_detail', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_2),sep=',')  
write(t(g_3),paste0('case3(dh)_stop', "_", 'sf(',iter,')', ".csv"),ncol=ncol(g_3),sep=',')

end_time <- Sys.time()

print(end_time - start_time)