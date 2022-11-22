# parameters in hierarchical_test changed when only subgroup is selected
sf_iteration_dh<-function(iter,seed_0,N_all,tau,lambda_PFS,lambda_OS,dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations,
                       weight_PFS,weight_OS,n_futility,threshold_futility,hz_PFS,hz_OS,alpha_s_ht_dh_1,alpha_s_ht_dh_2,
                       alpha_s_ht,alpha_f_ht,transition_ht,n_PFS_1s_ht,n_PFS_2s_ht,n_OS_3s_ht,n_PFS_1_ht,n_PFS_2_ht,n_OS_3_ht,
                       n_2_ht,n_3_ht,n_PFS_1s_ht_dh,n_PFS_2s_ht_dh,n_OS_3s_ht_dh,n_PFS_1_ht_dh,n_PFS_2_ht_dh,n_OS_3_ht_dh,
                       n_2_ht_dh,n_3_ht_dh){

  g33=NULL # hierarchical test  
  # detail: use to store significance status of 4 hypotheses of every subject
  detail=NULL
  stop=NULL
  
  for (j in 1:length(tau)){  # lambda_OS and futility threshold are paired with lambda_PFS
    for (i in 1:nrow(lambda_PFS)){
      seed_fail=NULL # seed not used to analyze
      # initial power array (reject PFS_full, reject PFS_sub, reject OS_full, reject OS_sub)
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

              ### hierarchical test ###
              # select.pop=full: use full group No. of events to cutoff data, select.pop!=full: use subgroup No. of events to cutoff data
              o_ht=hierarchical_test(alpha_1=alpha_s_ht,alpha_2=alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              # alpha_1, alpha_2: initial significance level for the first layer hypothesis test and second layer hypothesis test
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
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_both',g_ht))
              stop=rbind(stop,c('ht_both',stop_ht))
              
              k=k+1
            } else {
              
              # ht
              o_ht=hierarchical_test(alpha_s_ht,alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
             
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_both',g_ht))
              stop=rbind(stop,c('ht_both',stop_ht))
              
              k=k+1
            }
            
          } else if (hz_s>=threshold_futility[i,1] & hz_f<threshold_futility[i,2]){
            select.pop <- "full"
            
            if (k==1){
              
              # ht
              o_ht=hierarchical_test(rep(0,4),alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_full',g_ht))
              stop=rbind(stop,c('ht_full',stop_ht))
              
              k=k+1
            } else {
              
              # ht
              o_ht=hierarchical_test(rep(0,4),alpha_f_ht,transition=transition_ht,stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht,n_cut2s=n_PFS_2s_ht,n_cut3s=n_OS_3s_ht,n_cut1f=n_PFS_1_ht,
                                     n_cut2f=n_PFS_2_ht,n_cut3f=n_OS_3_ht,n_2_ht,n_3_ht,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_full',g_ht))
              stop=rbind(stop,c('ht_full',stop_ht))
              
              k=k+1
            }
            
          } else if (hz_s<threshold_futility[i,1] & hz_f>=threshold_futility[i,2]){
            select.pop <- "sub"
            
            if (k==1){
              
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
              
              
              # ht
              o_ht=hierarchical_test(alpha_s_ht_dh_1,alpha_s_ht_dh_2,transition=matrix(0,4,4),stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht_dh,n_cut2s=n_PFS_2s_ht_dh,n_cut3s=n_OS_3s_ht_dh,n_cut1f=n_PFS_1s_ht_dh,
                                     n_cut2f=n_PFS_2s_ht_dh,n_cut3f=n_OS_3s_ht_dh,n_2_ht_dh,n_3_ht_dh,weight_PFS[k,],weight_OS[k,])
              # test 2 subgroup hypotheses hierarchically: test PFS in the first layer and test OS in the second layer
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              
              if (length(g_ht)==1){  
                seed_fail=rbind(seed_fail,n)
                t=t-1
                break
              }
              
              g3[1,]=g3[1,]+g_ht
              detail=rbind(detail,c('ht_sub',g_ht))
              stop=rbind(stop,c('ht_sub',stop_ht))
              
              k=k+1
            } else {
              
              # ht
              o_ht=hierarchical_test(alpha_s_ht_dh_1,alpha_s_ht_dh_2,transition=matrix(0,4,4),stage1,select.pop,PFS,OS,
                                     n_cut1s=n_PFS_1s_ht_dh,n_cut2s=n_PFS_2s_ht_dh,n_cut3s=n_OS_3s_ht_dh,n_cut1f=n_PFS_1s_ht_dh,
                                     n_cut2f=n_PFS_2s_ht_dh,n_cut3f=n_OS_3s_ht_dh,n_2_ht_dh,n_3_ht_dh,weight_PFS[k,],weight_OS[k,])
              g_ht=o_ht[1,]
              stop_ht=o_ht[2,]
              
              g3[k,]=g3[k,]+g_ht
              detail=rbind(detail,c('ht_sub',g_ht))
              stop=rbind(stop,c('ht_sub',stop_ht))
              
              k=k+1
            }
            
          } else {
            select.pop <- "stop"
            
            if (k==1){
              
              detail=rbind(detail,c('ht_stop',rep(0,4)))
              stop=rbind(stop,c('ht_stop',rep(0,4)))
              k=k+1
            } else {
              
              detail=rbind(detail,c('ht_stop',rep(0,4)))
              stop=rbind(stop,c('ht_stop',rep(0,4)))
              k=k+1
            } 
            
          } # end if
          
        }
        print(t)
        t=t+1
      }
      
      for (ab in 1:nrow(weight_PFS)){
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
  g33=cbind(rep('ht',nrow(g33)),g33)
  g_1=rbind(title,g33)
  
  # detail
  title=c('method','Reject PFS-full group',
          'Reject PFS-subgroup','Reject OS-full group','Reject OS-subgroup',
          'Fraction','Hazard ratio(PFS-subgroup)','Hazard ratio(PFS-nonsubgroup)','Hazard ratio(OS-subgroup)',
          'Hazard ratio(OS-nonsubgroup)','Futility stop boundary(Subgroup)','Futility stop boundary(Full group)',
          'weight IA1 PFS stage1','weight IA1 PFS stage2','weight IA2 PFS stage1','weight IA2 PFS stage2',
          'weight FA PFS stage1','weight FA PFS stage2', 'weight IA1 OS stage1','weight IA1 OS stage2',
          'weight IA2 OS stage1','weight IA2 OS stage2','weight FA OS stage1','weight FA OS stage2')
  parameter=rbind(cbind(weight_PFS,weight_OS))
  # cbind(weight_PFS,weight_OS): ht.
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