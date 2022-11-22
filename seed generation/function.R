# change seed in each iteration of one simulation(other parameters are the same)
PFS_OS<-function(seed=0110,n,frac,lambda_PFS,lambda_OS,dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations){
  # generate PFS and OS using the same seed(same enrollTime, different failTime)
  # PFS
  set.seed(seed)
  x <-simPWSurv(n,  # total number of subjects enrolled 
                strata=tibble(Stratum=c("subgroup","non-subgroup"), p=c(frac, 1-frac)), # the fraction of subgroup
                block=c("Control","Experimental"),
                enrollRates=tibble(rate=enroll_rate,duration=durations),
                failRates=tibble(Stratum=rep(c("subgroup","non-subgroup"),2),
                                 period=c(1,1,1,1),
                                 Treatment=c(rep("Control",2), rep("Experimental",2)),
                                 duration=c(100,100,100,100),
                                 rate=lambda_PFS),   # hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)
                dropoutRates=tibble(Stratum=rep(c("subgroup","non-subgroup"),2),
                                    period=c(1,1,1,1),
                                    Treatment=c(rep("Control",2), rep("Experimental",2)),
                                    duration=c(100,100,100,100),
                                    rate=dropout_rate_PFS) # censoring rate 
  )
  
  
  # OS (change the median survival time, hazard ratio and dropout rate)
  set.seed(seed)
  y <-simPWSurv(n,  # total number of subjects enrolled 
                strata=tibble(Stratum=c("subgroup","non-subgroup"), p=c(frac, 1-frac)), # the fraction of subgroup
                block=c("Control","Experimental"),
                enrollRates=tibble(rate=enroll_rate,duration=durations), 
                failRates=tibble(Stratum=rep(c("subgroup","non-subgroup"),2),
                                 period=c(1,1,1,1),
                                 Treatment=c(rep("Control",2), rep("Experimental",2)),
                                 duration=c(100,100,100,100),
                                 rate=lambda_OS),   # hazard/failure rate(control_s, control_nonsub, treatment_s, treatment_nonsub)
                dropoutRates=tibble(Stratum=rep(c("subgroup","non-subgroup"),2),
                                    period=c(1,1,1,1),
                                    Treatment=c(rep("Control",2), rep("Experimental",2)),
                                    duration=c(100,100,100,100),
                                    rate=dropout_rate_OS) # censoring rate
  )
  
  return(rbind(x,y)) # x:PFS, y:OS
}


# failTime and dropoutTime follow exponential distribution
#whole=PFS_OS(seed,n=N_all,frac=tau[j],lambda_PFS=lambda_PFS[i,],lambda_OS=lambda_OS[i,],enroll_month,follow_up,dropout_rate_PFS,dropout_rate_OS,enroll_rate,durations)
#N_all=nrow(whole)/2
#PFS=whole[1:N_all,]
#OS=whole[(N_all+1):(2*N_all),]











# used in group_sequential 
gs_logrank<-function(data1,data2,n_cut,step){ 
  
  date_cut <- getCutDateForCount(filter(data1,Stratum=="subgroup"),n_cut) # cut data by No. of subgroup events 
  x <- cutData(data1,date_cut)
  
  x_s=x[x['Stratum']=='subgroup',]
  # calculate No. of real events for significance boundary
  n_data1_s=sum(x_s['event'])
  n_data1=sum(x['event'])
  

  
  # log-rank test of PFS
  x <- cutData(data2,date_cut)
  x_s=x[x['Stratum']=='subgroup',]
  # calculate No. of real events for significance boundary
  n_data2_s=sum(x_s['event'])
  n_data2=sum(x['event'])
  

  
  #p_2=c(pnorm(zPFS_2.F),pnorm(zPFS_2.S),pnorm(zOS_2.F),pnorm(zOS_2.S))
  if (step!=3){

    n_real=c(n_data1,n_cut,n_data2,n_data2_s)
  } else {

    n_real=c(n_data2,n_data2_s,n_data1,n_cut)
  }
 
  return(n_real)
}









########## group sequential ########## 
group_sequential<-function(alpha,transition,PFS,OS,n_PFS_1s,n_PFS_2s,n_OS_3s,n_2,n_3){

  n_2_0=n_2
  n_3_0=n_3
  power_ind=c(0,0,0,0)
  # IA1
  out=gs_logrank(PFS,OS,n_cut=n_PFS_1s,step=1)
  n_1=out

  # IA2
  out=gs_logrank(PFS,OS,n_cut=n_PFS_2s,step=2)
  n_2=out

  # FA
  out=gs_logrank(OS,PFS,n_cut=n_OS_3s,step=3)
  n_3=out
  
  # skip analyze current dataset if No. of events > prespecified events
  if (min(as.numeric(n_1<n_2_0))==0 | min(as.numeric(n_2[3:4]<n_3_0[3:4]))==0){
    return(0)
    break # do not run the following code
  }
  
  return(c(power_ind))
}










# used in adaptive_design & hierarchical_test
ad_logrank<-function(data1,data2,stage1,select.pop,w_data1,w_data2,n_cut,step){

  # cut data1
  if (select.pop!='full'){
    date_cut <- getCutDateForCount(filter(data1,Stratum=="subgroup"),n_cut) # cut data by No. of subgroup events 
    x <- cutData(data1,date_cut)
  } else {
    x <- cutDataAtCount(data1,n_cut) # use No. of events of full group to cut data when only full group is selected 
    date_cut=getCutDateForCount(data1,n_cut) # timepoint when No. of events happens
  }
  
  x_stage1=x[1:stage1,] # subjects recruited in stage 1 
  x_s1=x_stage1[x_stage1['Stratum']=='subgroup',]
  x_stage2=x[(1+stage1):nrow(x),] # subjects recruited in stage 2 
  x_s2=x_stage2[x_stage2['Stratum']=='subgroup',]
  # calculate No. of real events for significance boundary
  n_data1=sum(x_stage1['event'])+sum(x_stage2['event'])
  n_data1_s=sum(x_s1['event'])+sum(x_s2['event'])
  

  

  
  
  # cut data2
  x_stage1=cutData(data2[1:stage1,],date_cut) # subjects recruited in stage 1 
  x_s1=x_stage1[x_stage1['Stratum']=='subgroup',]
  x_stage2=cutData(data2[(1+stage1):nrow(x),],date_cut) # subjects recruited in stage 2 
  x_s2=x_stage2[x_stage2['Stratum']=='subgroup',]
  # calculate No. of real events for significance boundary
  n_data2_s=sum(x_s1['event'])+sum(x_s2['event'])
  n_data2=sum(x_stage1['event'])+sum(x_stage2['event'])
  

  
  
  if (step!=3){ 
    n_real=c(n_data1,n_data1_s,n_data2,n_data2_s)
 
  } else { 
    n_real=c(n_data2,n_data2_s,n_data1,n_data1_s) 
  
  } 
  
  

  
  return(n_real)
}












########## adaptive design ########## 
adaptive_design<-function(alpha,transition,stage1,select.pop,PFS,OS,n_cut1,n_cut2,n_cut3,n_2,n_3,weight_PFS,weight_OS){

  n_2_0=n_2
  n_3_0=n_3
  power_ind=c(0,0,0,0)
  # IA1
  out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[1:2],w_data2=weight_OS[1:2],n_cut=n_cut1,step=1)

  n_1=out

  # IA2
  out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[3:4],w_data2=weight_OS[3:4],n_cut=n_cut2,step=2)

  n_2=out

  # FA
  out=ad_logrank(OS,PFS,stage1,select.pop,w_data1=weight_OS[5:6],w_data2=weight_PFS[5:6],n_cut=n_cut3,step=3)
  n_3=out
  
  # skip analyze current dataset if No. of events > prespecified events
  if (min(as.numeric(n_1<n_2_0))==0 | min(as.numeric(n_2[3:4]<n_3_0[3:4]))==0){
    return(0)
    break # do not run the following code
  }
  
  
  
  return(c(power_ind))
  
}










########## hierarchical test ########## 
hierarchical_test<-function(alpha_1,alpha_2,transition,stage1,select.pop,PFS,OS,
                            n_cut1s,n_cut2s,n_cut3s,n_cut1f,n_cut2f,n_cut3f,n_2,n_3,weight_PFS,weight_OS){
  alpha_1_0=alpha_1
  alpha_2_0=alpha_2
  n_2_0=n_2
  n_3_0=n_3
  power_ind=c(0,0,0,0)
  a1s=a1f=a2s=a2f=a3s=a3f=rep(0,4)
  
  if (sum(alpha_1)>0){
    
    ############ subgroup ############
    # IA1
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[1:2],w_data2=weight_OS[1:2],n_cut=n_cut1s,step=1)
   
    n_1=out
    # IA2
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[3:4],w_data2=weight_OS[3:4],n_cut=n_cut2s,step=2)
    n_2=out
    # FA
    out=ad_logrank(OS,PFS,stage1,select.pop,w_data1=weight_OS[5:6],w_data2=weight_PFS[5:6],n_cut=n_cut3s,step=3)
    n_3=out
    
    # skip analyzing current dataset if No. of events > prespecified events
    if (min(as.numeric(n_1[which(alpha_1>0)]<n_2_0[which(alpha_1>0)]))==0){
      return(0)
      break # do not run the following code
    }
    if (sum(alpha_1[c(3,4)])>0){ # FA only test OS, so only need to check OS events
      if (min(as.numeric(n_2[which(alpha_1[c(3,4)]>0)+2]<n_3_0[which(alpha_1[c(3,4)]>0)+2]))==0){
        return(0)
        break # do not run the following code
      }
    }
    if (sum(alpha_2)>0){
      if (min(as.numeric(n_1[which(alpha_2>0)]<n_2_0[which(alpha_2>0)]))==0){
        return(0)
        break # do not run the following code
      }
      if (sum(alpha_2[c(3,4)])>0){ # FA only test OS, so only need to check OS events
        if (min(as.numeric(n_2[which(alpha_2[c(3,4)]>0)+2]<n_3_0[which(alpha_2[c(3,4)]>0)+2]))==0){
          return(0)
          break # do not run the following code
        }
      }
    }
    # make sure there is no No. of real events> No. of prespecified events when getting intersection boundary
  }
   
  
  if (sum(alpha_1_0)>0 & sum(alpha_2_0)>0){
    n_cut1f=n_cut1s
    n_cut2f=n_cut2s
    n_cut3f=n_cut3s
  }
  
  
  if (sum(alpha_2)>0){
    
    ############ full group ############
    # IA1
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[1:2],w_data2=weight_OS[1:2],n_cut=n_cut1f,step=1)
    n_1=out
    # IA2
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[3:4],w_data2=weight_OS[3:4],n_cut=n_cut2f,step=2)
    n_2=out
    # FA
    out=ad_logrank(OS,PFS,stage1,select.pop,w_data1=weight_OS[5:6],w_data2=weight_PFS[5:6],n_cut=n_cut3f,step=3)
    n_3=out
    
    # skip analyze current dataset if No. of events > prespecified events
    if (min(as.numeric(n_1[which(alpha_2>0)]<n_2_0[which(alpha_2>0)]))==0){
      return(0)
      break # do not run the following code
    }
    if (sum(alpha_2[c(3,4)])>0){ # FA only test OS, so only need to check OS events
      if (min(as.numeric(n_2[which(alpha_2[c(3,4)]>0)+2]<n_3_0[which(alpha_2[c(3,4)]>0)+2]))==0){
        return(0)
        break # do not run the following code
      }
    }
  }
    
  
  return(c(power_ind))
}