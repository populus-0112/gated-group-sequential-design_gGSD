source("function.R")

library(simtrial) # generate data
library(gsDesign) # IA and FA boundary
library(survival) # Cox model for hazard ratio
library(gMCP) # simes.test

# read seed number
seed_0=read.table('seed_0.csv',sep=',')
seed_0=as.matrix(seed_0)


# weights for p-value combination test
weight_PFS=rbind(c(0.5,0.5,0.5,0.5,0.5,0.5),c(0.5,0.5,0.5,0.5,0.5,0.5),c(0,0,0,0,0,0)) # IA1 w1, IA1 w2, IA2 w1, IA2 w2, FA w1, FA w2
weight_OS=rbind(c(0.7,0.3,0.7,0.3,0.7,0.3),c(0.3,0.7,0.3,0.7,0.3,0.7),c(0,0,0,0,0,0))  # IA1 w1, IA1 w2, IA2 w1, IA2 w2, FA w1, FA w2


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
threshold_futility=rbind(c(0.9,0.85),c(0.9,0.9))   # subgroup, full group
# get hazard ratio for non-subgroup using hazard ratio of subgroup & full group
hz_PFS=rbind(c(0.7,0.7),c(1,1)) # c(,): hazard ratio for treatment_s, treatment_nonsub
hz_OS=rbind(c(0.7,0.7),c(1,1)) # c(,): hazard ratio for treatment_s, treatment_nonsub
# log(2)/median survival time
OS_ns_median=5.7
OS_S_median=10.5
PFS_ns_median=3
PFS_S_median=4
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
alpha_s_ht=c(0,0.00987,0,0.01513) # when 2 hypotheses of full group are not tested, their alpha are set to be 0
# initialize alpha for 4 hypotheses
alpha_f_ht=c(0.01071,0,0.01429,0) # when 2 hypotheses of subgroup are not tested, their alpha are set to be 0
transition_ht=rbind(c(0,0,1,0),c(0,0,0,1),c(1,0,0,0),c(0,1,0,0))


# No. of OS events(subgroup)
n_OS_1s=227  # No. of OS events(subgroup) for IA1
n_OS_2s=314  # No. of OS events(subgroup) for IA2 
n_OS_3s=344  # No. of OS events(subgroup) for FA
n_OS_1s_ht=227
n_OS_2s_ht=313
n_OS_3s_ht=343
# No. of OS events(full group) 
n_OS_1=332
n_OS_2=446
n_OS_3=482
n_OS_1_ht=331
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
n_2_ht=c(n_PFS_2_ht,n_PFS_2s_ht,n_OS_2_ht,n_OS_2s_ht) # prespecified No. of events for IA2
n_3_ht=c(0,0,n_OS_3_ht,n_OS_3s_ht) # prespecified No. of events for FA

# futility analysis (PFS)
n_futility=253    #n_PFS_2*tau[j]  # No. of PFS events for futility


############# three methods ############# 
start_time <- Sys.time()
seed_fail=NULL # seed not used to analyze
for (j in 1:length(tau)){  # lambda_OS and futility threshold are paired with lambda_PFS
  for (i in 1:nrow(lambda_PFS)){

   
    n=0 # the count of seeds used(seeds used for power analysis + failed seeds due to No. of events > prespecified events) in the simulation
    t=1 # the count of real simulations for power and type I error 
    while (t<=iter){  
      n=n+1
      
      # generate PFS and OS
      seed=seed_0[n]
      whole=PFS_OS(seed,n=N_all,frac=tau[j],lambda_PFS=lambda_PFS[i,],lambda_OS=lambda_OS[i,],
                   dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations)
      PFS=whole[1:N_all,]
      OS=whole[(N_all+1):(2*N_all),]
      
      
      # calculate PFS hazard ratio for early stop
      # x <- cutDataAtCount(PFS,n_futility) # cut by No. of events in full group
      date_cut <- getCutDateForCount(filter(PFS,Stratum=="subgroup"),n_futility) # cut data by No. of subgroup events 
      x <- cutData(PFS,date_cut)
      stage1=nrow(x) # No. of subjects in stage 1 of adaptive design 
      # hazard ratio from Cox model
      g=coxph(Surv(tte, event) ~ Treatment, data = x)
      hz_f=exp(coef(g))
      x_s=x[x['Stratum']=='subgroup',]
      # hazard ratio from Cox model
      g=coxph(Surv(tte, event) ~ Treatment, data = x_s)
      hz_s=exp(coef(g))
      
      k=1
      # futility rule 
      if (hz_s<threshold_futility[i,1] & hz_f<threshold_futility[i,2]){
        select.pop <- "both"
        
        ### group sequential ###
        g_gs=group_sequential(alpha,transition,PFS,OS,n_PFS_1s,n_PFS_2s,n_OS_3s,n_2,n_3)
        
        if (length(g_gs)==1){  # return 0 from group_sequential
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        # ad
        g_ad=adaptive_design(alpha,transition,stage1,select.pop,PFS,OS,
                             n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ad)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        # ht_sf
        # select.pop=full: use full group No. of events to cutoff data, select.pop!=full: use subgroup No. of events to cutoff data
        g_ht=hierarchical_test(alpha_s_ht,alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                               n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                               n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ht)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        
        # ht_fs
        # select.pop=full: use full group No. of events to cutoff data, select.pop!=full: use subgroup No. of events to cutoff data
        g_ht=hierarchical_test(alpha_f_ht,alpha_s_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                               n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1s_ht,
                               n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ht)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
      } else if (hz_s>=threshold_futility[i,1] & hz_f<threshold_futility[i,2]){
        select.pop <- "full"
        
        ### group sequential ###
        g_gs=group_sequential(alpha,transition,PFS,OS,n_PFS_1s,n_PFS_2s,n_OS_3s,n_2,n_3)
        
        if (length(g_gs)==1){  # return 0 from group_sequential
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        # ad
        g_ad=adaptive_design(alpha=alpha_f,transition=transition_f,stage1,select.pop,
                             PFS,OS,n_cut1=n_PFS_1,n_cut2=n_PFS_2,n_cut3=n_OS_3,n_2,n_3,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ad)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        # ht_sf
        g_ht=hierarchical_test(rep(0,4),alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                               n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                               n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ht)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        

        
        # ht_fs
        g_ht=hierarchical_test(alpha_f_ht,rep(0,4),transition=transition_ht,stage1,select.pop,PFS,OS,
                               n_cut1s=n_PFS_1_ht,n_cut2s=n_PFS_2_ht,n_cut3s=n_OS_3_ht,n_cut1f=n_PFS_1s_ht,
                               n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ht)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
      } else if (hz_s<threshold_futility[i,1] & hz_f>=threshold_futility[i,2]){
        select.pop <- "sub"
        
        ### group sequential ###
        g_gs=group_sequential(alpha,transition,PFS,OS,n_PFS_1s,n_PFS_2s,n_OS_3s,n_2,n_3)
        
        if (length(g_gs)==1){  # return 0 from group_sequential
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        # regenerate data after futility analysis when only subgroup is selected
        # If futility analysis chooses subgroup, subjects recruited after futility analysis need to be regenerated using "frac=1"
        # n=N_all(not n_2rs) in order to keep enrollTime the same as original data
        whole_s=PFS_OS(seed,n=N_all,frac=1,lambda_PFS=lambda_PFS[i,],lambda_OS=lambda_OS[i,],
                       dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations)
        
        PFS_1=PFS[1:stage1,]
        n_fut_s=nrow(PFS_1[PFS_1['Stratum']=='subgroup',]) # No. of subjects recruited to subgroup before futility analysis
        
        # PFS and OS dataset for IA and FA
        # total No. of subjects recruited will change when only subgroup is selected
        # and subjects recruited to subgroup should be N_all*tau
        PFS=rbind(PFS[1:stage1,],whole_s[(stage1+1):(stage1+N_all*tau[j]-n_fut_s),]) 
        OS=rbind(OS[1:stage1,],whole_s[(N_all+stage1+1):(N_all+stage1+N_all*tau[j]-n_fut_s),]) 
        
        # ad
        g_ad=adaptive_design(alpha=alpha_s,transition=transition_s,stage1,select.pop,
                             PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ad)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        # ht_sf
        g_ht=hierarchical_test(alpha_s_ht,rep(0,4),transition=transition_ht,stage1,select.pop,PFS,OS,
                               n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                               n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ht)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
        
        # ht_fs
        g_ht=hierarchical_test(rep(0,4),alpha_s_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                               n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1s_ht,
                               n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
        
        if (length(g_ht)==1){  
          seed_fail=rbind(seed_fail,n)
          next
        }
        
      } else {
        select.pop <- "stop"
        
        ### group sequential ###
        g_gs=group_sequential(alpha,transition,PFS,OS,n_PFS_1s,n_PFS_2s,n_OS_3s,n_2,n_3)
        
        if (length(g_gs)==1){  # return 0 from group_sequential
          seed_fail=rbind(seed_fail,n)
          next
        }
        
      }
      
      
      print(t)
      t=t+1
    }
    
    }
  }


if (!is.null(seed_fail)){
  write(t(seed_fail),paste0('case1_seed_fail', ".csv"),ncol=1,sep=',')
}


end_time <- Sys.time()

print(end_time - start_time)