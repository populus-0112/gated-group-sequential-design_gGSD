sf_iteration<-function(iter,seed_0,N_all,tau,lambda_PFS,lambda_OS,dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations,
                       weight_PFS,weight_OS,n_futility,threshold_futility,hz_PFS,hz_OS,alpha,transition,alpha_f,alpha_s,
                       transition_f,transition_s,alpha_s_ht,alpha_f_ht,transition_ht,n_PFS_1s,n_PFS_2s,n_OS_3s,n_PFS_1,
                       n_PFS_2,n_OS_3,n_2,n_3,n_PFS_1s_ht,n_PFS_2s_ht,n_OS_3s_ht,n_PFS_1_ht,n_PFS_2_ht,n_OS_3_ht,
                       n_2_ht,n_3_ht){
  g11=NULL # group sequential
  g22=NULL # adaptive design
  g33=NULL # hierarchical test  
  # detail: use to store significance status of 4 hypotheses of every subject
  detail=NULL
  stop=NULL
  
  for (j in 1:length(tau)){  
    for (i in 1:nrow(lambda_PFS)){  # lambda_OS and futility threshold are paired with lambda_PFS
      seed_fail=NULL # seed not used to analyze
      # initial power array (reject PFS_full, reject PFS_sub, reject OS_full, reject OS_sub)
      g1=matrix(0,1,4) # group sequential, 4 columns: 4 hypotheses
      g2=matrix(0,nrow(weight_PFS),4) # adaptive design
      g3=matrix(0,nrow(weight_PFS),4) # hierarchical test
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
        #plot(survfit(Surv(tte, event) ~ Treatment, data = x))
        hz_f=exp(coef(g))
        x_s=x[x['Stratum']=='subgroup',]
        # hazard ratio from Cox model
        g=coxph(Surv(tte, event) ~ Treatment, data = x_s)
        #plot(survfit(Surv(tte, event) ~ Treatment, data = x_s))
        hz_s=exp(coef(g))
        
        k=1
        while (k<=nrow(weight_PFS)){  # weight_OS are paired with weight_PFS
          
          # futility rule 
          if (hz_s<threshold_futility[i,1] & hz_f<threshold_futility[i,2]){
            select.pop <- "both"
            
            if (k==1){
              ### group sequential ###
              o_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              g_gs=o_gs[1,]
              stop_gs=o_gs[2,]
              # alpha: initial significance level of PFS(F), PFS(S), OS(F), OS(S).
              # transition: transition matrix for alpha reallocation.
              # n_cut1,n_cut2,n_cut3: No. of events used to cutoff data.
              # n_2,n_3: No. of prespecified events at IA2 and FA for gs boundary
              # (each of them is a 4-dimensional vector: No. of events for PFS(F), PFS(S), OS(F), OS(S))
              
              if (length(g_gs)==1){  # return 0 from group_sequential: failed seeds due to No. of events > prespecified events
                seed_fail=rbind(seed_fail,n)
                t=t-1 # combined with 't=t+1' at the end of t loop: t will not add one if No. of events > prespecified events
                # that is: t is the real simulation times except failed seeds
                break # break k loop if No. of events > prespecified events
              }
              
              ### adaptive design ###
              o_ad=adaptive_design(alpha,transition,stage1,select.pop,PFS,OS,
                                   n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              g_ad=o_ad[1,]
              stop_ad=o_ad[2,]
              # n_cut1, n_cut2, n_cut3: No. of events used to cutoff data.
              # weight_PFS[k,],weight_OS[k,]: weights used in inverse normal combination test.
              
              if (length(g_ad)==1){  # return 0 from group_sequential: failed seeds due to No. of events > prespecified events
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              ### hierarchical test ###
              # select.pop=full: use full group No. of events to cutoff data, select.pop!=full: use subgroup No. of events to cutoff data
              o_ht=hierarchical_test(alpha_1=alpha_s_ht,alpha_2=alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              # alpha_1, alpha_2: initial significance level for the first layer hypothesis test and second layer hypothesis test
              # in sf_iteration function, always set subgroup alpha to alpha_1 and full group alpha to alpha_2 as subgroup hypothesis 
              # is tested first
              
              # n_cut1s,n_cut2s,n_cut3s: No. of events used to cutoff data for the first layer hypothesis test
              # n_cut1f,n_cut2f,n_cut3f: No. of events used to cutoff data for the second layer hypothesis test
              # n_2_ht,n_3_ht: No. of prespecified events at IA2 and FA for gs boundary
              # (each of them is a 4-dimensional vector: No. of events for PFS(F), PFS(S), OS(F), OS(S))
              # weight_PFS[k,],weight_OS[k,]: weights used in inverse normal combination test.
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              
              detail=rbind(detail,c('gs_both',g_gs))
              stop=rbind(stop,c('gs_both',stop_gs))
              
              g2[1,]=g2[1,]+g_ad
              
              detail=rbind(detail,c('ad_both',g_ad))
              stop=rbind(stop,c('ad_both',stop_ad))
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_both',g_ht))
              stop=rbind(stop,c('ht_both',stop_ht))
              
              k=k+1
            } else {
              # ad
              o_ad=adaptive_design(alpha,transition,stage1,select.pop,PFS,OS,
                                   n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              g_ad=o_ad[1,]
              stop_ad=o_ad[2,]
              # ht
              o_ht=hierarchical_test(alpha_s_ht,alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              
              g2[k,]=g2[k,]+g_ad
              
              detail=rbind(detail,c('ad_both',g_ad))
              stop=rbind(stop,c('ad_both',stop_ad))
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_both',g_ht))
              stop=rbind(stop,c('ht_both',stop_ht))
              
              k=k+1
            }
            
          } else if (hz_s>=threshold_futility[i,1] & hz_f<threshold_futility[i,2]){
            select.pop <- "full"
            
            if (k==1){
              ### group sequential ###
              o_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              g_gs=o_gs[1,]
              stop_gs=o_gs[2,]
              
              if (length(g_gs)==1){  # return 0 from group_sequential
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # ad
              o_ad=adaptive_design(alpha=alpha_f,transition=transition_f,stage1,select.pop,
                                   PFS,OS,n_cut1=n_PFS_1,n_cut2=n_PFS_2,n_cut3=n_OS_3,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              g_ad=o_ad[1,]
              stop_ad=o_ad[2,]
              # alpha and transition matrix are different when only full group is selected
              
              if (length(g_ad)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # ht
              o_ht=hierarchical_test(alpha_1=rep(0,4),alpha_2=alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              # subgroup is not selected, so set alpha_1=rep(0,4)
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              
              detail=rbind(detail,c('gs_full',g_gs))
              stop=rbind(stop,c('gs_full',stop_gs))
              
              g2[1,]=g2[1,]+g_ad
              detail=rbind(detail,c('ad_full',g_ad))
              stop=rbind(stop,c('ad_full',stop_ad))
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_full',g_ht))
              stop=rbind(stop,c('ht_full',stop_ht))
              
              k=k+1
            } else {
              # ad
              o_ad=adaptive_design(alpha=alpha_f,transition=transition_f,stage1,select.pop,
                                   PFS,OS,n_cut1=n_PFS_1,n_cut2=n_PFS_2,n_cut3=n_OS_3,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              g_ad=o_ad[1,]
              stop_ad=o_ad[2,]
              
              # ht
              o_ht=hierarchical_test(rep(0,4),alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              
              g2[k,]=g2[k,]+g_ad
              detail=rbind(detail,c('ad_full',g_ad))
              stop=rbind(stop,c('ad_full',stop_ad))
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_full',g_ht))
              stop=rbind(stop,c('ht_full',stop_ht))
              
              k=k+1
            }
            
          } else if (hz_s<threshold_futility[i,1] & hz_f>=threshold_futility[i,2]){
            select.pop <- "sub"
            
            if (k==1){
              ### group sequential ###
              o_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              g_gs=o_gs[1,]
              stop_gs=o_gs[2,]
              
              if (length(g_gs)==1){  # return 0 from group_sequential
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # regenerate data after futility analysis when only subgroup is selected
              # If futility analysis chooses subgroup, subjects recruited after futility analysis need to be regenerated using "frac=1"
              # n=N_all in order to keep enrollTime the same as original data
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
              o_ad=adaptive_design(alpha=alpha_s,transition=transition_s,stage1,select.pop,
                                   PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              g_ad=o_ad[1,]
              stop_ad=o_ad[2,]
              # alpha and transition matrix are different when only subgroup is selected
              
              if (length(g_ad)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # ht
              o_ht=hierarchical_test(alpha_1=alpha_s_ht,alpha_2=rep(0,4),transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              # full group is not selected, so set alpha_2=rep(0,4)
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              detail=rbind(detail,c('gs_sub',g_gs))
              stop=rbind(stop,c('gs_sub',stop_gs))
              
              g2[1,]=g2[1,]+g_ad
              detail=rbind(detail,c('ad_sub',g_ad))
              stop=rbind(stop,c('ad_sub',stop_ad))
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_sub',g_ht))
              stop=rbind(stop,c('ht_sub',stop_ht))
              
              k=k+1
            } else {
              # ad
              o_ad=adaptive_design(alpha=alpha_s,transition=transition_s,stage1,select.pop,PFS,OS,
                                   n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              g_ad=o_ad[1,]
              stop_ad=o_ad[2,]
              # ht
              o_ht=hierarchical_test(alpha_s_ht,rep(0,4),transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              
              g2[k,]=g2[k,]+g_ad
              
              detail=rbind(detail,c('ad_sub',g_ad))
              stop=rbind(stop,c('ad_sub',stop_ad))
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_sub',g_ht))
              stop=rbind(stop,c('ht_sub',stop_ht))
              
              k=k+1
            }
            
          } else {
            select.pop <- "stop"
            
            if (k==1){
              ### group sequential ###
              o_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              g_gs=o_gs[1,]
              stop_gs=o_gs[2,]
              
              if (length(g_gs)==1){  # return 0 from group_sequential
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              
              detail=rbind(detail,c('gs_stop',g_gs))
              stop=rbind(stop,c('gs_stop',stop_gs))
              
              detail=rbind(detail,c('ad_stop',rep(0,4)))
              stop=rbind(stop,c('ad_stop',rep(0,4)))
              
              detail=rbind(detail,c('ht_stop',rep(0,4)))
              stop=rbind(stop,c('ht_stop',rep(0,4)))
              
              k=k+1
            } else {
              detail=rbind(detail,c('ad_stop',rep(0,4)))
              stop=rbind(stop,c('ad_stop',rep(0,4)))
              
              detail=rbind(detail,c('ht_stop',rep(0,4)))
              stop=rbind(stop,c('ht_stop',rep(0,4)))
              k=k+1
            } 
            
          } # end if
          
        }
        print(t)
        t=t+1
      }
      
      # generate output table
      
      # group sequential has no weight parameter
      c25=paste(paste(paste(round(g1/iter,4),'(',sep=''),paste(g1,iter,sep="/"),sep=''),')',sep='')
      g=c(hz_PFS[i,],hz_OS[i,],c(NaN,NaN),tau[j],rep(NaN,12),c25)
      g11=rbind(g11,g)
      
      for (ab in 1:nrow(weight_PFS)){
        c24=paste(paste(paste(round(g2[ab,]/iter,4),'(',sep=''),paste(g2[ab,],iter,sep="/"),sep=''),')',sep='')
        g=rbind(c(hz_PFS[i,],hz_OS[i,],threshold_futility[i,],tau[j],weight_PFS[ab,],weight_OS[ab,],c24))
        g22=rbind(g22,g)
        c24=paste(paste(paste(round(g3[ab,]/iter,4),'(',sep=''),paste(g3[ab,],iter,sep="/"),sep=''),')',sep='')
        g=rbind(c(hz_PFS[i,],hz_OS[i,],threshold_futility[i,],tau[j],weight_PFS[ab,],weight_OS[ab,],c24))
        g33=rbind(g33,g)
      }
      if (!is.null(seed_fail)){
        write.table(seed_fail,file="seed_fail.csv",append=T,sep=',',col.names=F,row.names=F)
      }
    }
  }
  
  title=c('method','Hazard ratio(PFS-subgroup)','Hazard ratio(PFS-nonsubgroup)','Hazard ratio(OS-subgroup)',
          'Hazard ratio(OS-nonsubgroup)','Futility stop boundary(Subgroup)','Futility stop boundary(Full group)',
          'Fraction', 
          'weight IA1 PFS stage1','weight IA1 PFS stage2','weight IA2 PFS stage1','weight IA2 PFS stage2',
          'weight FA PFS stage1','weight FA PFS stage2', 'weight IA1 OS stage1','weight IA1 OS stage2',
          'weight IA2 OS stage1','weight IA2 OS stage2','weight FA OS stage1','weight FA OS stage2',
          'Reject PFS-full group','Reject PFS-subgroup','Reject OS-full group','Reject OS-subgroup')
  g11=cbind(rep('gs',nrow(g11)),g11)
  g22=cbind(rep('ad',nrow(g22)),g22)
  g33=cbind(rep('ht',nrow(g33)),g33)
  g_1=rbind(title,g11,g22,g33)
  
  # detail
  title=c('method','Reject PFS-full group',
          'Reject PFS-subgroup','Reject OS-full group','Reject OS-subgroup',
          'Fraction','Hazard ratio(PFS-subgroup)','Hazard ratio(PFS-nonsubgroup)','Hazard ratio(OS-subgroup)',
          'Hazard ratio(OS-nonsubgroup)','Futility stop boundary(Subgroup)','Futility stop boundary(Full group)',
          'weight IA1 PFS stage1','weight IA1 PFS stage2','weight IA2 PFS stage1','weight IA2 PFS stage2',
          'weight FA PFS stage1','weight FA PFS stage2', 'weight IA1 OS stage1','weight IA1 OS stage2',
          'weight IA2 OS stage1','weight IA2 OS stage2','weight FA OS stage1','weight FA OS stage2')
  parameter=rbind(rep(NaN,12),cbind(weight_PFS,weight_OS),cbind(weight_PFS,weight_OS))
  # rep(NaN,12): gs; the first cbind(weight_PFS,weight_OS): ad; the second cbind(weight_PFS,weight_OS): ht.
  parameter=do.call(rbind, replicate(iter, parameter, simplify=FALSE))
  parameter1=NULL
  for (i in 1:nrow(lambda_PFS)){
    parameter1=rbind(parameter1,cbind(do.call(rbind, replicate(nrow(parameter), c(hz_PFS[i,],hz_OS[i,],threshold_futility[i,]), simplify=FALSE))
                                      ,parameter))
  }
  parameter=NULL
  for (j in 1:length(tau)){
    parameter=rbind(parameter,cbind(do.call(rbind, replicate(nrow(parameter1), tau[j], simplify=FALSE))
                                    ,parameter1))
  }
  detail=cbind(detail,parameter)   
  g_2=rbind(title,detail)
  
  output=list("total" = g_1, "detail" = g_2, "stop"=stop)
  return(output)
}









fs_iteration<-function(iter,seed_0,N_all,tau,lambda_PFS,lambda_OS,dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations,
                       weight_PFS,weight_OS,n_futility,threshold_futility,hz_PFS,hz_OS,alpha,transition,alpha_f,alpha_s,
                       transition_f,transition_s,alpha_s_ht,alpha_f_ht,transition_ht,n_PFS_1s,n_PFS_2s,n_OS_3s,n_PFS_1,
                       n_PFS_2,n_OS_3,n_2,n_3,n_PFS_1s_ht,n_PFS_2s_ht,n_OS_3s_ht,n_PFS_1_ht,n_PFS_2_ht,n_OS_3_ht,
                       n_2_ht,n_3_ht){
  g11=NULL # group sequential
  g22=NULL # adaptive design
  g33=NULL # hierarchical test  
  # detail: use to store significance status of 4 hypotheses of every subject
  detail=NULL
  
  for (j in 1:length(tau)){  # lambda_OS and futility threshold are paired with lambda_PFS
    for (i in 1:nrow(lambda_PFS)){
      seed_fail=NULL # seed not used to analyze
      # initial power array (reject PFS_full, reject PFS_sub, reject OS_full, reject OS_sub)
      g1=matrix(0,1,4) # group sequential, 4 columns: 4 hypotheses
      g2=matrix(0,nrow(weight_PFS),4) # adaptive design
      g3=matrix(0,nrow(weight_PFS),4) # hierarchical test
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
        while (k<=nrow(weight_PFS)){  # weight_OS are paired with weight_PFS
          
          # futility rule 
          if (hz_s<threshold_futility[i,1] & hz_f<threshold_futility[i,2]){
            select.pop <- "both"
            
            if (k==1){
              ### group sequential ###
              g_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              # alpha: initial significance level of PFS(F), PFS(S), OS(F), OS(S).
              # transition: transition matrix for alpha reallocation.
              # n_cut1,n_cut2,n_cut3: No. of events used to cutoff data.
              # n_2,n_3: No. of prespecified events at IA2 and FA for gs boundary
              # (each of them is a 4-dimensional vector: No. of events for PFS(F), PFS(S), OS(F), OS(S))
              
              if (length(g_gs)==1){  # return 0 from group_sequential: failed seeds due to No. of events > prespecified events
                seed_fail=rbind(seed_fail,n)
                t=t-1 # combined with 't=t+1' at the end of t loop: t will not add one if No. of events > prespecified events
                # that is: t is the real simulation times except failed seeds
                break # break k loop if No. of events > prespecified events
              }
              
              ### adaptive design ###
              g_ad=adaptive_design(alpha,transition,stage1,select.pop,PFS,OS,
                                   n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              # n_cut1, n_cut2, n_cut3: No. of events used to cutoff data.
              # weight_PFS[k,],weight_OS[k,]: weights used in inverse normal combination test.
              
              if (length(g_ad)==1){  # return 0 from group_sequential: failed seeds due to No. of events > prespecified events
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              ### hierarchical test ###
              # select.pop=full: use full group No. of events to cutoff data, select.pop!=full: use subgroup No. of events to cutoff data
              g_ht=hierarchical_test(alpha_1=alpha_f_ht,alpha_2=alpha_s_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1s_ht,
                                     n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              # alpha_1, alpha_2: initial significance level for the first layer hypothesis test and second layer hypothesis test
              # in fs_iteration function, always set full group alpha to alpha_1 and subgroup alpha to alpha_2 as full group hypothesis 
              # is tested first
              
              # n_cut1s,n_cut2s,n_cut3s: No. of events used to cutoff data for the first layer hypothesis test
              # Since both groups selected, n_cut1s,n_cut2s,n_cut3s are set to be No. of subgroup events.
              # n_cut1f,n_cut2f,n_cut3f: No. of events used to cutoff data for the second layer hypothesis test
              # n_2_ht,n_3_ht: No. of prespecified events at IA2 and FA for gs boundary
              # (each of them is a 4-dimensional vector: No. of events for PFS(F), PFS(S), OS(F), OS(S))
              # weight_PFS[k,],weight_OS[k,]: weights used in inverse normal combination test.
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              
              detail=rbind(detail,c('gs_both',g_gs))
              
              g2[1,]=g2[1,]+g_ad
              
              detail=rbind(detail,c('ad_both',g_ad))
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_both',g_ht))
              
              k=k+1
            } else {
              # ad
              g_ad=adaptive_design(alpha,transition,stage1,select.pop,PFS,OS,
                                   n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              
              # ht
              g_ht=hierarchical_test(alpha_1=alpha_f_ht,alpha_2=alpha_s_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1s_ht,
                                     n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              
              g2[k,]=g2[k,]+g_ad
              
              detail=rbind(detail,c('ad_both',g_ad))
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_both',g_ht))
              
              k=k+1
            }
            
          } else if (hz_s>=threshold_futility[i,1] & hz_f<threshold_futility[i,2]){
            select.pop <- "full"
            
            if (k==1){
              ### group sequential ###
              g_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              
              if (length(g_gs)==1){  # return 0 from group_sequential
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # ad
              g_ad=adaptive_design(alpha=alpha_f,transition=transition_f,stage1,select.pop,
                                   PFS,OS,n_cut1=n_PFS_1,n_cut2=n_PFS_2,n_cut3=n_OS_3,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              # alpha and transition matrix are different when only full group is selected
              
              if (length(g_ad)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # ht
              g_ht=hierarchical_test(alpha_f_ht,rep(0,4),transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1_ht,n_cut2s=n_PFS_2_ht,n_cut3s=n_OS_3_ht,n_cut1f=n_PFS_1s_ht,
                                     n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              # subgroup is not selected, so set alpha_2=rep(0,4)
              # Since only full group selected, n_cut1s,n_cut2s,n_cut3s are set to be No. of full group events.
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              
              detail=rbind(detail,c('gs_full',g_gs))
              
              g2[1,]=g2[1,]+g_ad
              detail=rbind(detail,c('ad_full',g_ad))
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_full',g_ht))
              
              k=k+1
            } else {
              # ad
              g_ad=adaptive_design(alpha=alpha_f,transition=transition_f,stage1,select.pop,
                                   PFS,OS,n_cut1=n_PFS_1,n_cut2=n_PFS_2,n_cut3=n_OS_3,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              
              
              # ht
              g_ht=hierarchical_test(alpha_f_ht,rep(0,4),transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1_ht,n_cut2s=n_PFS_2_ht,n_cut3s=n_OS_3_ht,n_cut1f=n_PFS_1s_ht,
                                     n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              
              g2[k,]=g2[k,]+g_ad
              detail=rbind(detail,c('ad_full',g_ad))
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_full',g_ht))
              
              k=k+1
            }
            
          } else if (hz_s<threshold_futility[i,1] & hz_f>=threshold_futility[i,2]){
            select.pop <- "sub"
            
            if (k==1){
              ### group sequential ###
              g_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              
              if (length(g_gs)==1){  # return 0 from group_sequential
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # regenerate data after futility analysis when only subgroup is selected
              # If futility analysis chooses subgroup, subjects recruited after futility analysis need to be regenerated using "frac=1"
              # n=N_all in order to keep enrollTime the same as original data
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
              # alpha and transition matrix are different when only subgroup is selected
              
              if (length(g_ad)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              # ht
              g_ht=hierarchical_test(rep(0,4),alpha_s_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1s_ht,
                                     n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              # full group is not selected, so set alpha_1=rep(0,4)
              # Since sub group is selected, n_cut1s,n_cut2s,n_cut3s are set to be No. of subgroup events.
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              detail=rbind(detail,c('gs_sub',g_gs))
              
              g2[1,]=g2[1,]+g_ad
              detail=rbind(detail,c('ad_sub',g_ad))
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_sub',g_ht))
              
              k=k+1
            } else {
              # ad
              g_ad=adaptive_design(alpha=alpha_s,transition=transition_s,stage1,select.pop,PFS,OS,
                                   n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3,weight_PFS[k,],weight_OS[k,])
              
              # ht
              g_ht=hierarchical_test(rep(0,4),alpha_s_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1s_ht,
                                     n_cut2f=n_PFS_2s_ht,n_cut3f=n_OS_3s_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              
              g2[k,]=g2[k,]+g_ad
              
              detail=rbind(detail,c('ad_sub',g_ad))
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_sub',g_ht))
              
              k=k+1
            }
            
          } else {
            select.pop <- "stop"
            
            if (k==1){
              ### group sequential ###
              g_gs=group_sequential(alpha,transition,PFS,OS,n_cut1=n_PFS_1s,n_cut2=n_PFS_2s,n_cut3=n_OS_3s,n_2,n_3)
              
              if (length(g_gs)==1){  # return 0 from group_sequential
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g1=g1+g_gs
              
              detail=rbind(detail,c('gs_stop',g_gs))
              
              detail=rbind(detail,c('ad_stop',rep(0,4)))
              
              detail=rbind(detail,c('ht_stop',rep(0,4)))
              k=k+1
            } else {
              detail=rbind(detail,c('ad_stop',rep(0,4)))
              
              detail=rbind(detail,c('ht_stop',rep(0,4)))
              k=k+1
            } 
            
          } # end if
          
        }
        print(t)
        t=t+1
      }
      
      # generate output table
      
      # group sequential has no weight parameter
      c25=paste(paste(paste(round(g1/iter,4),'(',sep=''),paste(g1,iter,sep="/"),sep=''),')',sep='')
      g=c(hz_PFS[i,],hz_OS[i,],c(NaN,NaN),tau[j],rep(NaN,12),c25)
      g11=rbind(g11,g)
      
      for (ab in 1:nrow(weight_PFS)){
        c24=paste(paste(paste(round(g2[ab,]/iter,4),'(',sep=''),paste(g2[ab,],iter,sep="/"),sep=''),')',sep='')
        g=rbind(c(hz_PFS[i,],hz_OS[i,],threshold_futility[i,],tau[j],weight_PFS[ab,],weight_OS[ab,],c24))
        g22=rbind(g22,g)
        c24=paste(paste(paste(round(g3[ab,]/iter,4),'(',sep=''),paste(g3[ab,],iter,sep="/"),sep=''),')',sep='')
        g=rbind(c(hz_PFS[i,],hz_OS[i,],threshold_futility[i,],tau[j],weight_PFS[ab,],weight_OS[ab,],c24))
        g33=rbind(g33,g)
      }
      if (!is.null(seed_fail)){
        write.table(seed_fail,file="seed_fail.csv",append=T,sep=',',col.names=F,row.names=F)
      }
    }
  }
  
  title=c('method','Hazard ratio(PFS-subgroup)','Hazard ratio(PFS-nonsubgroup)','Hazard ratio(OS-subgroup)',
          'Hazard ratio(OS-nonsubgroup)','Futility stop boundary(Subgroup)','Futility stop boundary(Full group)',
          'Fraction', 
          'weight IA1 PFS stage1','weight IA1 PFS stage2','weight IA2 PFS stage1','weight IA2 PFS stage2',
          'weight FA PFS stage1','weight FA PFS stage2', 'weight IA1 OS stage1','weight IA1 OS stage2',
          'weight IA2 OS stage1','weight IA2 OS stage2','weight FA OS stage1','weight FA OS stage2',
          'Reject PFS-full group','Reject PFS-subgroup','Reject OS-full group','Reject OS-subgroup')
  g11=cbind(rep('gs',nrow(g11)),g11)
  g22=cbind(rep('ad',nrow(g22)),g22)
  g33=cbind(rep('ht',nrow(g33)),g33)
  g_1=rbind(title,g11,g22,g33)
  
  # detail
  title=c('method','Reject PFS-full group',
          'Reject PFS-subgroup','Reject OS-full group','Reject OS-subgroup',
          'Fraction','Hazard ratio(PFS-subgroup)','Hazard ratio(PFS-nonsubgroup)','Hazard ratio(OS-subgroup)',
          'Hazard ratio(OS-nonsubgroup)','Futility stop boundary(Subgroup)','Futility stop boundary(Full group)',
          'weight IA1 PFS stage1','weight IA1 PFS stage2','weight IA2 PFS stage1','weight IA2 PFS stage2',
          'weight FA PFS stage1','weight FA PFS stage2', 'weight IA1 OS stage1','weight IA1 OS stage2',
          'weight IA2 OS stage1','weight IA2 OS stage2','weight FA OS stage1','weight FA OS stage2')
  parameter=rbind(rep(NaN,12),cbind(weight_PFS,weight_OS),cbind(weight_PFS,weight_OS)) 
  # rep(NaN,12): gs; the first cbind(weight_PFS,weight_OS): ad; the second cbind(weight_PFS,weight_OS): ht.
  parameter=do.call(rbind, replicate(iter, parameter, simplify=FALSE))
  parameter1=NULL
  for (i in 1:nrow(lambda_PFS)){
    parameter1=rbind(parameter1,cbind(do.call(rbind, replicate(nrow(parameter), c(hz_PFS[i,],hz_OS[i,],threshold_futility[i,]), simplify=FALSE))
                                      ,parameter))
  }
  parameter=NULL
  for (j in 1:length(tau)){
    parameter=rbind(parameter,cbind(do.call(rbind, replicate(nrow(parameter1), tau[j], simplify=FALSE))
                                    ,parameter1))
  }
  detail=cbind(detail,parameter)   
  g_2=rbind(title,detail)
  
  output=list("total" = g_1, "detail" = g_2)
  return(output)
}



