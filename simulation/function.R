simes.test<-function(z1,z2){
  p1=pnorm(z1)
  p2=pnorm(z2)
  p=min(2*min(p1,p2),max(p1,p2))
  return(p)
}


# change seed in each iteration (other parameters are the same)
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
  
  
  # OS
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









# sum(a)=0: no alpha reallocation
alpha_reallocation<-function(a,z,alpha,transition){
  d=sum(a[which(alpha>0)]) # No. of significant hypotheses that need to reallocate alpha
  if (d>0){
    delete=NULL
    for (i in 1:d){
      w=(1-pnorm(z))/alpha  # Inf when alpha=0(hypothesis is not tested) 
      w[which(a==0)]=100 # do not pick non-significant hypothesis
      r=which(w==min(w))[1] # pick priority hypothesis to reallocate its alpha
      delete=c(delete,r)
      transition_new=transition
      transition_new[,delete]=0 # deleted hypotheses do not accept reallocated alpha
      b=rep(0,length(w))
      b[r]=1 # the location of significant hypothesis to reallocate its alpha
      alpha=(b*alpha)%*%transition_new+(1-b)*alpha
      # (b*alpha)%*%transition_new: reallocated alpha, (1-b)*alpha: old alpha except the one being reallocated
        for (j in seq(1,nrow(transition))[!seq(1,nrow(transition)) %in% delete]){  
          for (l in seq(1,nrow(transition))[!seq(1,nrow(transition)) %in% c(delete,j)]){
            if (!is.na(transition[j,r]*transition[r,j])){
              if (transition[j,r]*transition[r,j]!=1){
                transition[j,l]=(transition[j,l]+transition[j,r]*transition[r,l])/(1-transition[j,r]*transition[r,j])
              } else if (transition[j,r]*transition[r,j]==1){
                transition[j,l]=0 # alpha doesn't reallocate to others if alpha is only designed to reallocate to r initially 
              }
            }
          }
        }
        transition[r,]=0
        transition[,r]=0
        z[r]=-100 # make sure do not pick up to reallocate alpha next time
    }
  }
  return(rbind(alpha,transition))
}  # alpha should be 0 for significant hypothesis(a=1)






# used in group_sequential 
# For each efficacy analysis, n_cut is related to which survival type(PFS/OS), then set this survival type's dataset    
# to be data1. Dataset of another survival type is data2.
gs_logrank<-function(data1,data2,n_cut,step){ 
  
  date_cut <- getCutDateForCount(filter(data1,Stratum=="subgroup"),n_cut) # cut data by No. of subgroup events 
  x <- cutData(data1,date_cut)
  
  x_s=x[x['Stratum']=='subgroup',]
  # calculate No. of real events for significance boundary
  n_data1_s=sum(x_s['event'])
  n_data1=sum(x['event'])
  
  
  x_f=tensurv(x,txval="Experimental")
  zdata1.F <- with(x_f,sum(OminusE)/sqrt(sum(Var)))
  
  x_s=tensurv(x_s,txval="Experimental")
  zdata1.S <- with(x_s,sum(OminusE)/sqrt(sum(Var)))
  
  # log-rank test for data2
  x <- cutData(data2,date_cut)
  x_s=x[x['Stratum']=='subgroup',]
  # calculate No. of real events for significance boundary
  n_data2_s=sum(x_s['event'])
  n_data2=sum(x['event'])
  
  x_f=tensurv(x,txval="Experimental")
  zdata2.F <- with(x_f,sum(OminusE)/sqrt(sum(Var)))
  
  x_s=tensurv(x_s,txval="Experimental")
  zdata2.S <- with(x_s,sum(OminusE)/sqrt(sum(Var)))
  
  # step=1: IA1, step=2: IA2, step=3: FA
  if (step!=3){  # IA1 and IA2 use No. of PFS events to trigger efficacy analysis
    z=c(zdata1.F,zdata1.S,zdata2.F,zdata2.S)
    z=-z
    n_real=c(n_data1,n_data1_s,n_data2,n_data2_s)
  } else {  # FA uses No. of PFS events to trigger efficacy analysis
    z=c(zdata2.F,zdata2.S,zdata1.F,zdata1.S)
    z=-z
    n_real=c(n_data2,n_data2_s,n_data1,n_data1_s)
  }
  
  return(rbind(z,n_real))
}






# used in group_sequential
alpha_reallocation_iteration<-function(z,n_1,n_2,n_2_0,n_3,n_3_0,alpha,transition,step){
  
  a=rep(0,4)
  
  for (l in c(1,2)){
    if (alpha[l]>0 & step!=3){ # PFS hypotheses are not tested at FA
      # redo test for non-significant hypothesis using reallocated alpha
      g <- gsDesign(k=2, test.type=1, alpha=alpha[l], maxn.IPlan = n_2_0[l], n.I=c(n_1[l],n_2[l]), sfu=sfLDOF)  
      bound=g$upper$bound[step] # z boundary 
      a[l]=as.numeric(z[l]>bound)
    }
  }
  
  for (l in c(3,4)){
    if (alpha[l]>0){ 
      # redo test for non-significant hypothesis using reallocated alpha
      g <- gsDesign(k=3, test.type=1, alpha=alpha[l], maxn.IPlan = n_3_0[l], n.I=c(n_1[l],n_2[l],n_3[l]), sfu=sfLDOF)  
      bound=g$upper$bound[step] # z boundary 
      a[l]=as.numeric(z[l]>bound)
    }
  }
  
  # alpha reallocation
  a_old=a-1
  while(sum(abs(a_old-a))>0 & sum(a)>0){ # stop if no hypothesis becomes significant after using updated alpha
    # sum(a) denotes No. of significant hypotheses which continue to be tested 
    a_old=a
    out=alpha_reallocation(a,z,alpha,transition)
    alpha=out[1,]
    transition=out[2:nrow(out),]
    
    if (step!=1){  # IA2 and FA: the first time alpha reallocation happens, return to IA1 
      break
    }
    
    for (l in c(1,2)){
      if (alpha[l]>0 & step!=3){ # PFS hypotheses are not tested at FA
        # redo test for non-significant hypothesis using reallocated alpha
        g <- gsDesign(k=2, test.type=1, alpha=alpha[l], maxn.IPlan = n_2_0[l], n.I=c(n_1[l],n_2[l]), sfu=sfLDOF)  
        bound=g$upper$bound[step] # z boundary 
        a[l]=as.numeric(z[l]>bound)
      }
    }
    
    for (l in c(3,4)){
      if (alpha[l]>0){ 
        # redo test for non-significant hypothesis using reallocated alpha
        g <- gsDesign(k=3, test.type=1, alpha=alpha[l], maxn.IPlan = n_3_0[l], n.I=c(n_1[l],n_2[l],n_3[l]), sfu=sfLDOF)  
        bound=g$upper$bound[step] # z boundary 
        a[l]=as.numeric(z[l]>bound)
      }
    }
  }
  return(rbind(alpha,a,transition))
}





########## group sequential ########## 
group_sequential<-function(alpha,transition,PFS,OS,n_cut1,n_cut2,n_cut3,n_2,n_3){
  
  n_2_0=n_2
  n_3_0=n_3
  power_ind=c(0,0,0,0)
  stop=c(3,3,3,3)
  # IA1 and IA2 use No. of PFS events to trigger efficacy analysis, FA uses No. of PFS events to trigger efficacy analysis.
  # IA1
  out=gs_logrank(PFS,OS,n_cut=n_cut1,step=1) # step=1: IA1, step=2: IA2, step=3: FA
  z1=out[1,]
  n_1=out[2,]
  # IA2
  out=gs_logrank(PFS,OS,n_cut=n_cut2,step=2)
  z2=out[1,]
  n_2=out[2,]
  # FA
  out=gs_logrank(OS,PFS,n_cut=n_cut3,step=3)
  z3=out[1,]
  n_3=out[2,]
  
  # skip analyzing current dataset if No. of events > prespecified events
  if (min(as.numeric(n_1<n_2_0))==0 | min(as.numeric(n_2[3:4]<n_3_0[3:4]))==0){
    return(0)
    break # do not run the following code
  }
  
  
  a1=a2=a3=rep(0,4) # record significance status of 4 hypotheses in IA1, IA2 and FA
  alpha_old=alpha-1 
  while(sum(alpha)>0 & sum(abs(alpha_old-alpha))>0 & sum(as.numeric(a1>0)+as.numeric(a2>0)+as.numeric(a3>0))<4){
    # IA1
    out=alpha_reallocation_iteration(z=z1,n_1,n_2_0,n_2_0,n_3_0,n_3_0,alpha,transition,step=1) # step=1:IA1, step=2: IA2, step=3: FA
    alpha=out[1,]
    a1=a1+out[2,]
    transition=out[3:nrow(out),]
    for (ii in 1:4){
      if (a1[ii]==1){stop[ii]=1}
    }
    
    alpha_old=alpha
    
    # IA2
    if (sum(alpha)>0){ # sum(alpha)>0: there is hypothesis not significant at IA1
      out=alpha_reallocation_iteration(z=z2,n_1,n_2,n_2_0,n_3_0,n_3_0,alpha,transition,step=2) 
      alpha=out[1,]
      a2=a2+out[2,]
      transition=out[3:nrow(out),]
      for (ii in 1:4){
        if (a2[ii]==1){stop[ii]=2}
      }
      
      if (sum(out[2,])>0){  # alpha reallocated in IA2, then go back to IA1
        next
      }
      
      # FA: only test OS
      if (sum(alpha[c(3,4)])>0){ # if at least 1 of 2 OS hypotheses is non-significant at IA1 and IA2
        out=alpha_reallocation_iteration(z=z3,n_1,n_2,n_2_0,n_3,n_3_0,alpha,transition,step=3)
        alpha=out[1,]
        a3=a3+out[2,]
        transition=out[3:nrow(out),]
        for (ii in 1:4){
          if (a3[ii]==1){stop[ii]=3}
        }
      }
    }
  }
  
  # power for 4 hypotheses respectively
  power_ind=power_ind+as.numeric(a1>0)+as.numeric(a2>0)+as.numeric(a3>0)
  out=rbind(power_ind,stop)
  #return(power_ind)
  return(out)
}










# used in adaptive_design & hierarchical_test
# For each efficacy analysis, n_cut is related to which survival type(PFS/OS), then set this survival type's dataset    
# to be data1. Dataset of another survival type is data2.
ad_logrank<-function(data1,data2,stage1,select.pop,w_data1,w_data2,n_cut,step){
  if (sum(w_data1)==0){  # no prespecified weight
    weight="update"
  } else {
    weight="constant"
  }
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
  
  if (weight=="update" & select.pop=='sub'){ 
    w_data1=c(sum(x_s1['event']),sum(x_s2['event']))/n_data1_s
  } else if (weight=="update" & select.pop!='sub'){
    w_data1=c(sum(x_stage1['event']),sum(x_stage2['event']))/n_data1
  }
  
  # data1 log-rank test, full group
  # stage 1
  xy=tensurv(x_stage1,txval="Experimental")
  zdata1_1.F <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  # stage 2
  xy=tensurv(x_stage2,txval="Experimental")
  zdata1_2.F <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  
  
  # data1 log-rank test, subgroup
  # stage 1
  xy=tensurv(x_s1,txval="Experimental")
  zdata1_1.S <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  # stage 2
  xy=tensurv(x_s2,txval="Experimental")
  zdata1_2.S <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  
  
  # cut data2
  x_stage1=cutData(data2[1:stage1,],date_cut) # subjects recruited in stage 1 
  x_s1=x_stage1[x_stage1['Stratum']=='subgroup',]
  x_stage2=cutData(data2[(1+stage1):nrow(x),],date_cut) # subjects recruited in stage 2 
  x_s2=x_stage2[x_stage2['Stratum']=='subgroup',]
  # calculate No. of real events for significance boundary
  n_data2_s=sum(x_s1['event'])+sum(x_s2['event'])
  n_data2=sum(x_stage1['event'])+sum(x_stage2['event'])
  
  # data2 log-rank test, full group
  # stage 1
  xy=tensurv(x_stage1,txval="Experimental")
  zdata2_1.F <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  # stage 2
  xy=tensurv(x_stage2,txval="Experimental")
  zdata2_2.F <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  
  
  # data2 log-rank test, subgroup
  # stage 1
  xy=tensurv(x_s1,txval="Experimental")
  zdata2_1.S <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  # stage 2
  xy=tensurv(x_s2,txval="Experimental")
  zdata2_2.S <- with(xy,sum(OminusE)/sqrt(sum(Var)))
  
  
  # inverse normal combination test
  
  if (weight=="update"){ 
    if (select.pop=='sub'){ 
      w_data2=c(sum(x_s1['event']),sum(x_s2['event']))/n_data2_s   
    } else if (select.pop!='sub'){
      w_data2=c(sum(x_stage1['event']),sum(x_stage2['event']))/n_data2
    }
  } 
  
  z_1 <- sqrt(w_data1[1]) * qnorm(1 - pnorm(zdata1_1.F)) + sqrt(w_data1[2]) * qnorm(1 - pnorm(zdata1_2.F))
  z_2 <- sqrt(w_data1[1]) * qnorm(1 - pnorm(zdata1_1.S)) + sqrt(w_data1[2]) * qnorm(1 - pnorm(zdata1_2.S))
  z_3 <- sqrt(w_data2[1]) * qnorm(1 - pnorm(zdata2_1.F)) + sqrt(w_data2[2]) * qnorm(1 - pnorm(zdata2_2.F))
  z_4 <- sqrt(w_data2[1]) * qnorm(1 - pnorm(zdata2_1.S)) + sqrt(w_data2[2]) * qnorm(1 - pnorm(zdata2_2.S))
  
  # step=1: IA1, step=2: IA2, step=3: FA
  if (step!=3){ # IA1 and IA2 use No. of PFS events to trigger efficacy analysis
    n_real=c(n_data1,n_data1_s,n_data2,n_data2_s)
    z=c(z_1,z_2,z_3,z_4)
  } else { # FA uses No. of PFS events to trigger efficacy analysis
    n_real=c(n_data2,n_data2_s,n_data1,n_data1_s) 
    z=c(z_3,z_4,z_1,z_2)
  } 
  
  
  if (select.pop=="both"){
    
    # intersection of 2 hypotheses for data1
    p_data1=c(simes.test(zdata1_1.F,zdata1_1.S), # stage 1 intersection
              simes.test(zdata1_2.F,zdata1_2.S)) # stage 2 intersection
    
    # intersection of 2 hypotheses for data2
    p_data2=c(simes.test(zdata2_1.F,zdata2_1.S), # stage 1 intersection
              simes.test(zdata2_2.F,zdata2_2.S)) # stage 2 intersection
    
  } else if (select.pop=="sub"){
    
    # intersection of 2 hypotheses for data1
    p_data1=c(simes.test(zdata1_1.F,zdata1_1.S), # stage 1 intersection
              pnorm(zdata1_2.S)) # stage 2 only has subgroup data
    
    # intersection of 2 hypotheses for data2
    p_data2=c(simes.test(zdata2_1.F,zdata2_1.S), # stage 1 intersection
              pnorm(zdata2_2.S)) # stage 2 only has subgroup data
    
  } else if (select.pop=="full"){
    
    # intersection of 2 hypotheses for data1
    p_data1=c(simes.test(zdata1_1.F,zdata1_1.S), # stage 1 intersection
              pnorm(zdata1_2.F)) # stage 2 only has full group data
    
    # intersection of 2 hypotheses for data2
    p_data2=c(simes.test(zdata2_1.F,zdata2_1.S), # stage 1 intersection
              pnorm(zdata2_2.F)) # stage 2 only has full group data
    
  }
  
  # combination test of 2 stage intersection hypotheses
  z_data1=sqrt(w_data1[1]) * qnorm(1 - p_data1[1]) + sqrt(w_data1[2]) * qnorm(1 - p_data1[2])
  z_data2=sqrt(w_data2[1]) * qnorm(1 - p_data2[1]) + sqrt(w_data2[2]) * qnorm(1 - p_data2[2])
  
  if (step!=3){ # IA1 and IA2 use No. of PFS events to trigger efficacy analysis
    intersection=c(z_data1,z_data2,0,0)
  } else { # FA uses No. of PFS events to trigger efficacy analysis
    intersection=c(z_data2,z_data1,0,0)
  } 
  
  return(rbind(z,n_real,intersection))
}






# used in adaptive_design & hierarchical_test
# intersection(useful at IA1): used to record previous significance status of intersection hypotheses 
# At IA1: skip intersection hypothesis test if the it is already significant previously
ad_alpha_reallocation_iteration<-function(z,z_PFS,z_OS,n_1,n_2,n_2_0,n_3,n_3_0,alpha,transition,step,intersection){  
  
  if (step!=1){ # intersection is only useful at IA1, set intersection to be 0 for IA2 and FA
    intersection=rep(0,4) 
  }
  
  bound=c(0,0,0,0)
  
  a=c(0,0,0,0)
  
  for (l in c(1,2)){
    if (alpha[l]>0 & step!=3){  # PFS are not tested at FA
      # redo test for non-significant hypothesis using reallocated alpha
      g <- gsDesign(k=2, test.type=1, alpha=alpha[l], maxn.IPlan = n_2_0[l], n.I=c(n_1[l],n_2[l]), sfu=sfLDOF)  
      bound[l]=g$upper$bound[step] # z boundary 
      a[l]=as.numeric(z[l]>bound[l])
    }
  }
  
  for (l in c(3,4)){
    if (alpha[l]>0){
      # redo test for non-significant hypothesis using reallocated alpha
      g <- gsDesign(k=3, test.type=1, alpha=alpha[l], maxn.IPlan = n_3_0[l], n.I=c(n_1[l],n_2[l],n_3[l]), sfu=sfLDOF) 
      bound[l]=g$upper$bound[step] 
      a[l]=as.numeric(z[l]>bound[l])
    }
  }
  
  inter_PFS=inter_OS=-100
  # use minimum boundary(subgroup or full group) for intersection hypothesis when subgroup & full group are selected
  if (bound[1]*bound[2]>0){
    inter_PFS=min(bound[2],bound[1])
  } else if (bound[1]==0 & bound[2]>0 & intersection[1]==0){ 
    # bound[i]=0(that is alpha[i]=0) & intersection[i]=1: hypothesis i was tested and significant & intersection test is significant
    # under the above illustrated condition, intersection don't need to be retested (inter_PFS=-100 to avoid retest)
    # IA2 and FA: intersection is set to be 0. That is always test intersection hypothesis at IA2 and FA
    inter_PFS=bound[2]
  } else if (bound[2]==0 & bound[1]>0 & intersection[2]==0){
    inter_PFS=bound[1]
  }
  
  if (bound[3]*bound[4]>0){
    inter_OS=min(bound[3],bound[4])
  } else if (bound[3]==0 & bound[4]>0 & intersection[3]==0){
    inter_OS=bound[4]
  } else if (bound[4]==0 & bound[3]>0 & intersection[4]==0){
    inter_OS=bound[3]
  }
  
  
  if (z_PFS<inter_PFS){ 
    a[1:2]=c(0,0) # PFS hypotheses are not significant if the intersection of PFS hypotheses are not significant 
  }
  if (z_OS<inter_OS){ 
    a[3:4]=c(0,0) # OS hypotheses are not significant if the intersection of OS hypotheses are not significant 
  }

  # prepare for IA1: if PFS/OS intersection test is significant in previous alpha reallocation, then it doesn't need to be retested
  # because z statistic doesn't change and alpha can not decrease for non-significant hypothesis
  if (step==1){  
    intersection[which(a!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
    intersection[c(1,2)]=max(intersection[c(1,2)]) # the significance status of PFS intersection hypothesis should be the same
    intersection[c(3,4)]=max(intersection[c(3,4)]) # the significance status of OS intersection hypothesis should be the same
  }
  
  
  # alpha reallocation
  a_old=a-1
  while(sum(abs(a_old-a))>0 & sum(a)>0){ # stop if no hypothesis becomes significant after using updated alpha
    a_old=a
    out=alpha_reallocation(a,z,alpha,transition)
    alpha=out[1,]
    transition=out[2:nrow(out),]
    
    if (step!=1){  # IA2 and FA: the first time alpha reallocation happens, return to IA1 
      break
    }
    
    for (l in c(1,2)){
      if (alpha[l]>0 & step!=3){ 
        # redo test for non-significant hypothesis using reallocated alpha
        g <- gsDesign(k=2, test.type=1, alpha=alpha[l], maxn.IPlan = n_2_0[l], n.I=c(n_1[l],n_2[l]), sfu=sfLDOF)  
        bound[l]=g$upper$bound[step] # z boundary 
        a[l]=as.numeric(z[l]>bound[l])
      }
    }
    
    for (l in c(3,4)){
      if (alpha[l]>0){
        # redo test for non-significant hypothesis using reallocated alpha
        g <- gsDesign(k=3, test.type=1, alpha=alpha[l], maxn.IPlan = n_3_0[l], n.I=c(n_1[l],n_2[l],n_3[l]), sfu=sfLDOF) 
        bound[l]=g$upper$bound[step] 
        a[l]=as.numeric(z[l]>bound[l])
      }
    }
    
    # retest intersection hypotheses using reallocated alpha
    # make sure significant status is first decided by intersection hypothesis
    # reallocate alpha only when intersection hypothesis & individual hypothesis are both significant
    inter_PFS=inter_OS=-100
    # use minimum boundary(subgroup or full group) for intersection hypothesis when subgroup & full group are selected
    if (bound[1]*bound[2]>0){
      inter_PFS=min(bound[2],bound[1])
    } else if (bound[1]==0 & bound[2]>0 & intersection[1]==0){ 
      # bound[i]=0(that is alpha[i]=0) & intersection[i]=1: hypothesis i was tested and significant & intersection test is significant 
      inter_PFS=bound[2]
    } else if (bound[2]==0 & bound[1]>0 & intersection[2]==0){
      inter_PFS=bound[1]
    }
    
    if (bound[3]*bound[4]>0){
      inter_OS=min(bound[3],bound[4])
    } else if (bound[3]==0 & bound[4]>0 & intersection[3]==0){
      inter_OS=bound[4]
    } else if (bound[4]==0 & bound[3]>0 & intersection[4]==0){
      inter_OS=bound[3]
    }
    
    if (z_PFS<inter_PFS){ 
      a[1:2]=c(0,0) # PFS hypotheses are not significant if the intersection of PFS hypotheses are not significant 
    }
    if (z_OS<inter_OS){ 
      a[3:4]=c(0,0) # OS hypotheses are not significant if the intersection of OS hypotheses are not significant 
    }
    
    if (step==1){  
      intersection[which(a!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
      intersection[c(1,2)]=max(intersection[c(1,2)]) # the significance status of PFS intersection hypothesis should be the same
      intersection[c(3,4)]=max(intersection[c(3,4)]) # the significance status of OS intersection hypothesis should be the same
    }
  }
  # a[i]=1 means significant hypothesis
  
  return(rbind(alpha,a,transition))
}








########## adaptive design ########## 
adaptive_design<-function(alpha,transition,stage1,select.pop,PFS,OS,n_cut1,n_cut2,n_cut3,n_2,n_3,weight_PFS,weight_OS){
  
  n_2_0=n_2
  n_3_0=n_3
  power_ind=c(0,0,0,0)
  stop=c(3,3,3,3)
  # IA1
  out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[1:2],w_data2=weight_OS[1:2],n_cut=n_cut1,step=1)
  z1=out[1,]
  n_1=out[2,]
  z_PFS_IA1=out[3,1]
  z_OS_IA1=out[3,2]
  # IA2
  out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[3:4],w_data2=weight_OS[3:4],n_cut=n_cut2,step=2)
  z2=out[1,]
  n_2=out[2,]
  z_PFS_IA2=out[3,1]
  z_OS_IA2=out[3,2]
  # FA
  out=ad_logrank(OS,PFS,stage1,select.pop,w_data1=weight_OS[5:6],w_data2=weight_PFS[5:6],n_cut=n_cut3,step=3)
  z3=out[1,]
  n_3=out[2,]
  z_PFS_FA=out[3,1]
  z_OS_FA=out[3,2]
  
  # skip analyzing current dataset if No. of events > prespecified events
  if (min(as.numeric(n_1<n_2_0))==0 | min(as.numeric(n_2[3:4]<n_3_0[3:4]))==0){
    return(0)
    break # do not run the following code
  }
  
  a1=a2=a3=rep(0,4)
  
  intersection=rep(0,4)
  alpha_old=alpha-1
  while(sum(alpha)>0 & sum(abs(alpha_old-alpha))>0 & sum(as.numeric(a1>0)+as.numeric(a2>0)+as.numeric(a3>0))<4){
    # IA1
    out=ad_alpha_reallocation_iteration(z=z1,z_PFS=z_PFS_IA1,z_OS=z_OS_IA1,n_1,n_2_0,n_2_0,n_3_0,n_3_0,alpha,transition,step=1,intersection)
    alpha=out[1,]
    a1=a1+out[2,]
    transition=out[3:nrow(out),]
    intersection[which(out[2,]!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
    for (ii in 1:4){
      if (a1[ii]==1){stop[ii]=1}
    }
    
    alpha_old=alpha
    
    # IA2
    if (sum(alpha)>0){
      out=ad_alpha_reallocation_iteration(z=z2,z_PFS=z_PFS_IA2,z_OS=z_OS_IA2,n_1,n_2,n_2_0,n_3_0,n_3_0,alpha,transition,step=2,intersection)
      alpha=out[1,]
      a2=a2+out[2,]
      transition=out[3:nrow(out),]
      for (ii in 1:4){
        if (a2[ii]==1){stop[ii]=2}
      }
      
      if (sum(out[2,])>0){  # alpha reallocated in IA2
        next
      }
      
      # FA
      if (sum(a2[3:4])<2 & sum(alpha[c(3,4)])>0){
        out=ad_alpha_reallocation_iteration(z=z3,z_PFS=z_PFS_FA,z_OS=z_OS_FA,n_1,n_2,n_2_0,n_3,n_3_0,alpha,transition,step=3,intersection)
        alpha=out[1,]
        a3=a3+out[2,]
        transition=out[3:nrow(out),]
      }
    }
  }
  
  # power for 4 hypotheses respectively
  power_ind=power_ind+as.numeric(a1>0)+as.numeric(a2>0)+as.numeric(a3>0)
  out=rbind(power_ind,stop)
  #return(power_ind)
  return(out)
}




# used in hierarchical_test when starting to do test from IA2
# intersection(useful at IA2): used to record previous significance status of intersection hypotheses 
# At IA2: skip intersection hypothesis test if the it is already significant previously
alpha_reallocation_iteration_ht<-function(z,z_PFS,z_OS,n_1,n_2,n_2_0,alpha,transition,step,intersection){  
  
  if (step!=2){ # intersection is only useful at IA2, set intersection to be 0 for FA
    intersection=rep(0,4) 
  }
  
  bound=c(0,0,0,0)
  
  a=c(0,0,0,0)
  
  for (l in c(1,2)){
    if (alpha[l]>0 & step!=3){
      bound[l]=qnorm(1-alpha[l])
      a[l]=as.numeric(z[l]>bound[l])
    }
  }
  
  for (l in c(3,4)){
    if (alpha[l]>0){
      # redo test for non-significant hypothesis using reallocated alpha
      g <- gsDesign(k=2, test.type=1, alpha=alpha[l], maxn.IPlan = n_2_0[l], n.I=c(n_1[l],n_2[l]), sfu=sfLDOF) 
      bound[l]=g$upper$bound[step-1] 
      a[l]=as.numeric(z[l]>bound[l])
    }
  }
  
  inter_PFS=inter_OS=-100
  # use minimum boundary(subgroup or full group) for intersection hypothesis when subgroup & full group are selected
  if (bound[1]*bound[2]>0){
    inter_PFS=min(bound[2],bound[1])
  } else if (bound[1]==0 & bound[2]>0 & intersection[1]==0){
    inter_PFS=bound[2]
  } else if (bound[2]==0 & bound[1]>0 & intersection[2]==0){
    inter_PFS=bound[1]
  }
  
  if (bound[3]*bound[4]>0){
    inter_OS=min(bound[3],bound[4])
  } else if (bound[3]==0 & bound[4]>0 & intersection[3]==0){
    inter_OS=bound[4]
  } else if (bound[4]==0 & bound[3]>0 & intersection[4]==0){
    inter_OS=bound[3]
  }
  
  if (z_PFS<inter_PFS){ 
    a[1:2]=c(0,0) # PFS hypotheses are not significant if the intersection of PFS hypotheses are not significant 
  }
  if (z_OS<inter_OS){ 
    a[3:4]=c(0,0) # OS hypotheses are not significant if the intersection of OS hypotheses are not significant 
  }
  
  # prepare for IA2: if PFS/OS intersection test is significant in previous alpha reallocation, then it doesn't need to be retested
  # because z statistic doesn't change and alpha can not decrease for non-significant hypothesis
  if (step==2){  
    intersection[which(a!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
    intersection[c(1,2)]=max(intersection[c(1,2)]) # the significance status of PFS intersection hypothesis should be the same
    intersection[c(3,4)]=max(intersection[c(3,4)]) # the significance status of OS intersection hypothesis should be the same
  }
  
  
  # alpha reallocation
  a_old=a-1
  while(sum(abs(a_old-a))>0 & sum(a)>0){ # stop if no hypothesis becomes significant after using updated alpha
    a_old=a
    out=alpha_reallocation(a,z,alpha,transition)
    alpha=out[1,]
    transition=out[2:nrow(out),]
    
    if (step!=2){  # FA: the first time alpha reallocation happens, return to IA2 
      break
    }
    
    for (l in c(1,2)){
      if (alpha[l]>0 & step!=3){
        bound[l]=qnorm(1-alpha[l])
        a[l]=as.numeric(z[l]>bound[l])
      }
    }
    
    for (l in c(3,4)){
      if (alpha[l]>0){
        # redo test for non-significant hypothesis using reallocated alpha
        g <- gsDesign(k=2, test.type=1, alpha=alpha[l], maxn.IPlan = n_2_0[l], n.I=c(n_1[l],n_2[l]), sfu=sfLDOF) 
        bound[l]=g$upper$bound[step-1] 
        a[l]=as.numeric(z[l]>bound[l])
      }
    }
    
    # retest intersection hypotheses using reallocated alpha
    # make sure significant status is first decided by intersection hypothesis
    # reallocate alpha only when intersection hypothesis & individual hypothesis are both significant
    inter_PFS=inter_OS=-100
    # use minimum boundary(subgroup or full group) for intersection hypothesis when subgroup & full group are selected
    if (bound[1]*bound[2]>0){
      inter_PFS=min(bound[2],bound[1])
    } else if (bound[1]==0 & bound[2]>0 & intersection[1]==0){
      inter_PFS=bound[2]
    } else if (bound[2]==0 & bound[1]>0 & intersection[2]==0){
      inter_PFS=bound[1]
    }
    
    if (bound[3]*bound[4]>0){
      inter_OS=min(bound[3],bound[4])
    } else if (bound[3]==0 & bound[4]>0 & intersection[3]==0){
      inter_OS=bound[4]
    } else if (bound[4]==0 & bound[3]>0 & intersection[4]==0){
      inter_OS=bound[3]
    }
    
    if (z_PFS<inter_PFS){ 
      a[1:2]=c(0,0) # PFS hypotheses are not significant if the intersection of PFS hypotheses are not significant 
    }
    if (z_OS<inter_OS){ 
      a[3:4]=c(0,0) # OS hypotheses are not significant if the intersection of OS hypotheses are not significant 
    }
    
    # prepare for IA2: if PFS/OS intersection test is significant in previous alpha reallocation, then it doesn't need to be retested
    # because z statistic doesn't change and alpha can not decrease for non-significant hypothesis
    if (step==2){  
      intersection[which(a!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
      intersection[c(1,2)]=max(intersection[c(1,2)]) # the significance status of PFS intersection hypothesis should be the same
      intersection[c(3,4)]=max(intersection[c(3,4)]) # the significance status of OS intersection hypothesis should be the same
    }
    
  }
  # a[i]=1 means significant hypothesis
  
  return(rbind(alpha,a,transition))
}







########## hierarchical test ########## 
hierarchical_test<-function(alpha_1,alpha_2,transition,stage1,select.pop,PFS,OS,
                            n_cut1s,n_cut2s,n_cut3s,n_cut1f,n_cut2f,n_cut3f,n_2,n_3,weight_PFS,weight_OS){
  alpha_1_0=alpha_1
  alpha_2_0=alpha_2
  n_2_0=n_2
  n_3_0=n_3
  power_ind=c(0,0,0,0)
  stop=rep(3,4)
  a1s=a1f=a2s=a2f=a3s=a3f=rep(0,4)
  
  if (sum(alpha_1)>0){
    
    ############ first layer hypothesis testing ############
    # IA1
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[1:2],w_data2=weight_OS[1:2],n_cut=n_cut1s,step=1)
    z1=out[1,]
    z1[which(alpha_1==0)]=0 # in order to let a=0 in ad_alpha_reallocation_iteration function
    n_1=out[2,]
    z_PFS_IA1=out[3,1]
    z_OS_IA1=out[3,2]
    # IA2
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[3:4],w_data2=weight_OS[3:4],n_cut=n_cut2s,step=2)
    z2=out[1,]
    z2[which(alpha_1==0)]=0 # in order to let a=0 in ad_alpha_reallocation_iteration function
    n_2=out[2,]
    z_PFS_IA2=out[3,1]
    z_OS_IA2=out[3,2]
    # FA
    out=ad_logrank(OS,PFS,stage1,select.pop,w_data1=weight_OS[5:6],w_data2=weight_PFS[5:6],n_cut=n_cut3s,step=3)
    z3=out[1,]
    z3[which(alpha_1==0)]=0 # in order to let a=0 in ad_alpha_reallocation_iteration function
    n_3=out[2,]
    z_PFS_FA=out[3,1]
    z_OS_FA=out[3,2]
    
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
    
    intersection=rep(0,4)
    alpha=alpha_1+alpha_2
    alpha_old=alpha-1
    # alpha_1+alpha_2: make sure to use minimum boundary(subgroup or full group) for intersection hypothesis when subgroup & full group are selected
    # z1[which(alpha_1==0)]=0: in order to let a=0 in ad_alpha_reallocation_iteration function and do not reallocate alpha_2
    # transition matrix guarantees alpha_1 won't be reallocated to the position of alpha_2>0
    while(sum(alpha-alpha_2)>0 & sum(abs(alpha_old-alpha))>0 & sum(as.numeric(a1s>0)+as.numeric(a2s>0)+as.numeric(a3s>0))<2){
      # IA1
      out=ad_alpha_reallocation_iteration(z=z1,z_PFS=z_PFS_IA1,z_OS=z_OS_IA1,n_1,n_2_0,n_2_0,n_3_0,n_3_0,alpha,transition,step=1,intersection)
      alpha=out[1,]
      a1s=a1s+out[2,]
      transition=out[3:nrow(out),]
      intersection[which(out[2,]!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
      for (ii in 1:4){
        if (a1s[ii]==1){stop[ii]=1}
      }
      alpha_old=alpha
      
      # IA2
      if (sum(alpha-alpha_2)>0){   
        out=ad_alpha_reallocation_iteration(z=z2,z_PFS=z_PFS_IA2,z_OS=z_OS_IA2,n_1,n_2,n_2_0,n_3_0,n_3_0,alpha,transition,step=2,intersection)
        alpha=out[1,]
        a2s=a2s+out[2,]
        transition=out[3:nrow(out),]
        for (ii in 1:4){
          if (a2s[ii]==1){stop[ii]=2}
        }
        
        if (sum(out[2,])>0){  # alpha reallocated in IA2
          next
        }
        
        # FA: only test OS
        if (sum((a1s+a2s)[c(3,4)])<1 & sum((alpha-alpha_2)[c(3,4)])>0){ 
          out=ad_alpha_reallocation_iteration(z=z3,z_PFS=z_PFS_FA,z_OS=z_OS_FA,n_1,n_2,n_2_0,n_3,n_3_0,alpha,transition,step=3,intersection)
          alpha=out[1,]
          a3s=a3s+out[2,]
          transition=out[3:nrow(out),]
        }
      }
    }
  }
  
  if (sum(alpha_1)>0 & sum(alpha_2)>0){
    n_cut1f=n_cut1s
    n_cut2f=n_cut2s
    n_cut3f=n_cut3s
  }
  
  
  if (sum(alpha_2)>0){
    
    ############ second layer hypothesis testing ############
    # IA1
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[1:2],w_data2=weight_OS[1:2],n_cut=n_cut1f,step=1)
    z1=out[1,]
    z1[which(alpha_2==0)]=0 # in order to let a=0 in ad_alpha_reallocation_iteration function
    n_1=out[2,]
    z_PFS_IA1=out[3,1]
    z_OS_IA1=out[3,2]
    # IA2
    out=ad_logrank(PFS,OS,stage1,select.pop,w_data1=weight_PFS[3:4],w_data2=weight_OS[3:4],n_cut=n_cut2f,step=2)
    z2=out[1,]
    z2[which(alpha_2==0)]=0 # in order to let a=0 in ad_alpha_reallocation_iteration function
    n_2=out[2,]
    z_PFS_IA2=out[3,1]
    z_OS_IA2=out[3,2]
    # FA
    out=ad_logrank(OS,PFS,stage1,select.pop,w_data1=weight_OS[5:6],w_data2=weight_PFS[5:6],n_cut=n_cut3f,step=3)
    z3=out[1,]
    z3[which(alpha_2==0)]=0 # in order to let a=0 in ad_alpha_reallocation_iteration function
    n_3=out[2,]
    z_PFS_FA=out[3,1]
    z_OS_FA=out[3,2]
    
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

    
    intersection=rep(0,4)
    alpha=alpha_1+alpha_2
    # alpha_1+alpha_2: make sure to use minimum boundary(subgroup or full group) for intersection hypothesis when subgroup & full group are selected
    alpha_old=alpha-1
    while(sum(alpha-alpha_1)>0 & sum(abs(alpha_old-alpha))>0 & sum(as.numeric(a1f>0)+as.numeric(a2f>0)+as.numeric(a3f>0))<2){
      # IA1
      out=ad_alpha_reallocation_iteration(z=z1,z_PFS=z_PFS_IA1,z_OS=z_OS_IA1,n_1,n_2_0,n_2_0,n_3_0,n_3_0,alpha,transition,step=1,intersection)
      alpha=out[1,]
      a1f=a1f+out[2,]
      transition=out[3:nrow(out),]
      intersection[which(out[2,]!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
      for (ii in 1:4){
        if (a1f[ii]==1){stop[ii]=1}
      }
      alpha_old=alpha
      
      # IA2
      if (sum(alpha-alpha_1)>0){
        out=ad_alpha_reallocation_iteration(z=z2,z_PFS=z_PFS_IA2,z_OS=z_OS_IA2,n_1,n_2,n_2_0,n_3_0,n_3_0,alpha,transition,step=2,intersection)
        alpha=out[1,]
        a2f=a2f+out[2,]
        transition=out[3:nrow(out),]
        for (ii in 1:4){
          if (a2f[ii]==1){stop[ii]=2}
        }
        if (sum(out[2,])>0){  # alpha reallocated in IA2
          next
        }
        
        # FA: only test OS
        if (sum((a1f+a2f)[c(3,4)])<1 & sum((alpha-alpha_1)[c(3,4)])>0){
          out=ad_alpha_reallocation_iteration(z=z3,z_PFS=z_PFS_FA,z_OS=z_OS_FA,n_1,n_2,n_2_0,n_3,n_3_0,alpha,transition,step=3,intersection)
          alpha=out[1,]
          a3f=a3f+out[2,]
          transition=out[3:nrow(out),]
        }
      }
    }
  }
  
  
  if (sum(alpha_1)>0 & sum(a1s)>0){
    
    # power for 4 hypotheses respectively
    power_ind=power_ind+as.numeric(a1s>0)+as.numeric(a2s>0)+as.numeric(a3s>0)+as.numeric(a1f>0)+as.numeric(a2f>0)+as.numeric(a3f>0)
  } else if (sum(alpha_1)>0 & sum(a1s)==0 & sum(a2s)>0){
    
    if (sum(alpha_2)>0){
      a2f=a3f=rep(0,4)
      ############ second layer IA2 & FA (IA2 and FA group sequential) ############
      intersection=rep(0,4)
      alpha=alpha_2+alpha_1
      alpha_old=alpha-1
      while(sum(alpha-alpha_1)>0 & sum(abs(alpha_old-alpha))>0 & sum(as.numeric(a2f>0)+as.numeric(a3f>0))<2){
        # IA2
        out=alpha_reallocation_iteration_ht(z=z2,z_PFS=z_PFS_IA2,z_OS=z_OS_IA2,n_2,n_3_0,n_3_0,alpha,transition,step=2,intersection)
        alpha=out[1,]
        a2f=a2f+out[2,]
        transition=out[3:nrow(out),]
        intersection[which(out[2,]!=0)]=1 # =1 if this hypothesis is significant: be useful for the intersection test when select.pop='both'
        for (ii in 1:4){
          if (a2f[ii]==1){stop[ii]=2}
        }
        alpha_old=alpha
        
        # FA
        if (sum(a2f[c(3,4)])<1 & sum((alpha-alpha_1)[c(3,4)])>0){
          out=alpha_reallocation_iteration_ht(z=z3,z_PFS=z_PFS_FA,z_OS=z_OS_FA,n_2,n_3,n_3_0,alpha,transition,step=3,intersection)
          alpha=out[1,]
          a3f=a3f+out[2,]
          transition=out[3:nrow(out),]
        }
      }
    }
    
    # power for 4 hypotheses respectively
    power_ind=power_ind+as.numeric(a2s>0)+as.numeric(a3s>0)+as.numeric(a2f>0)+as.numeric(a3f>0)
    
  } else if (sum(alpha_1)>0 & sum(a1s)==0 & sum(a2s)==0 & sum(a3s)>0){
    
    ############ second layer FA ############
    boundary=qnorm(1-alpha_2)
    a3f=as.numeric(z3>boundary)
    
    # power for 4 hypotheses respectively
    power_ind=power_ind+as.numeric(a3s>0)+as.numeric(a3f>0)
    
  } else if(sum(alpha_1)>0 & sum(a1s)==0 & sum(a2s)==0 & sum(a3s)==0){
    
    # power for 4 hypotheses respectively
    power_ind=power_ind+0
    
  } else if (sum(alpha_1)==0 & sum(alpha_2)>0){
    
    # power for 4 hypotheses respectively
    power_ind=power_ind+as.numeric(a1f>0)+as.numeric(a2f>0)+as.numeric(a3f>0)
  }
  out=rbind(power_ind,stop)
  #return(power_ind)
  return(out)
}