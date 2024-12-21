#return decimals
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


target.range<-function(H1,H0){
  mean.diff=(H1-H0)/2
  span=seq(-0.5,0.5,0.5)*mean.diff
  target=span+H1
  return(target)
}

target.range.lg<-function(H1,H0){
  mean.diff=(H1-H0)/2
  span=seq(-1,1,0.5)*mean.diff
  target=span+H1
  return(target)
}

##################################################################################
# Completely randomized data
###################################################################################

single_one_stage_size=function(tau,tau_e=0,theta ,rho, alpha, beta, direction, var_type){
  
  if(direction==1 & var_type==1) {
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    n1=((1+1/rho)*((tau^2)*(z_a+z_b)^2))/theta^2   
    n2=rho*n1
    n1=2*ceiling(n1/2)
    n2=2*ceiling(n2/2)
    n=n1+n2
    n_val=n1+n2
  }
  
  
  if(direction==2 & var_type==1) {
    
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    n1= ((1+1/rho)*((tau^2)*(z_a+z_b)^2))/theta^2   
    n2=rho*n1
    n1=2*ceiling(n1/2)
    n2=2*ceiling(n2/2)
    n=n1+n2
    n_val=n1+n2
  }
  
  
  if(direction==1 & var_type==2) {
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    a=tau_e/tau
    n1=((1+((a^2)/rho) )*((tau^2)*(z_a+z_b)^2))/theta^2   
    n2=rho*n1
    n1=2*ceiling(n1/2)
    n2=2*ceiling(n2/2)
    n=n1+n2
    n_val=n1+n2
  }
  
  
  if(direction==2 & var_type==2) {
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    a=tau_e/tau
    n1=((1+((a^2)/rho) )*((tau^2)*(z_a+z_b)^2))/theta^2   
    n2=rho*n1
    n1=2*ceiling(n1/2)
    n2=2*ceiling(n2/2)
    n=n1+n2
    n_val=n1+n2
  }
  
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val
  )
  
  return(final_result)
}


#------------------------------------------------------------------------------

single_one_stage_dat<-function(tau, tau_e=0, theta, rho, alpha, beta, direction, var_type){
  target=target.range.lg(H1=theta,H0=0)
  nround=max(decimalplaces(theta),decimalplaces(tau))+1
  dat<- data.frame(rbind(round(target,nround), sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist(single_one_stage_size(tau, tau_e, target[x], rho, alpha, beta,
                                   direction, var_type))
      
    ))))
  
  
  rownames(dat)=c("Treatment Effect (Offsite)", "Sample Size (Control Arm)",
                  "Sample Size (Experimental Arm)", "Total sample size")
  dat[c(2,3,4),dat[1,]<0]=NA  #difference in mean suppose to > 0
  return(dat)
}


################################################################################


hybrid_one_stage_size=function(tau1, tau1_e=0, tau2, tau2_e=0, theta, pi_IF, r, rho,
                               alpha, beta, direction, var_type){
  
  
  if(direction==1 & var_type==1) {
    #------------------------------------
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    #n1=round(4*(tau1^2)*(z_a+z_b)^2/((theta^2)*(1+ (r/(k^2) ) )))
    n1=((2+rho+1/rho)*(tau1^2)*(z_a+z_b)^2/((theta^2)))* (1/(1 + r*((1+pi_IF)^2)/(k))) 
    #n1
    n2=r*n1
    #n2
    n1_t = ceiling(n1 * (rho/(1+rho)))
    n1_c = ceiling(n1 * (1/(1+rho)))
    n2_t = ceiling(n2 * (rho/(1+rho)))
    n2_c = ceiling(n2 * (1/(1+rho)))
    
    
    n1=(n1_t+n1_c)
    n2=(n2_t+n2_c)
    n_val=n1+n2
  }
  
  
  if(direction==2 & var_type==1) {
    
    #------------------------------------
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    #n1=round(4*(tau1^2)*(z_a+z_b)^2/((theta^2)*(1+ (r/(k^2) ) )))
    n1=((2+rho+1/rho)*(tau1^2)*(z_a+z_b)^2/((theta^2)))* (1/(1 + r*((1+pi_IF)^2)/(k)))  
    #n1
    n2=r*n1
    #n2
    n1_t = ceiling(n1 * (rho/(1+rho)))
    n1_c = ceiling(n1 * (1/(1+rho)))
    n2_t = ceiling(n2 * (rho/(1+rho)))
    n2_c = ceiling(n2 * (1/(1+rho)))
    
    
    n1=(n1_t+n1_c)
    n2=(n2_t+n2_c)
    n_val=n1+n2
  }
  
  
  if(direction==1 & var_type==2) {
    #------------------------------------
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(rho+1)
    
    parta = 1/(1+rho + (a1^2)*((1+rho)/rho) )
    partb= 1/(1+rho + (a2^2)*((1+rho)/rho) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    
    n1=((tau1^2)*((z_a+z_b)^2)/(theta^2))* der
    #n1
    n2=r*n1
    #n2
    n1_t = ceiling(n1 * (rho/(1+rho)))
    n1_c = ceiling(n1 * (1/(1+rho)))
    n2_t = ceiling(n2 * (rho/(1+rho)))
    n2_c = ceiling(n2 * (1/(1+rho)))
    
    
    n1=(n1_t+n1_c)
    n2=(n2_t+n2_c)
    n_val=n1+n2
  }
  
  
  if(direction==2 & var_type==2) {
    #------------------------------------
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(rho+1)
    
    parta = 1/(1+rho + (a1^2)*((1+rho)/rho) )
    partb= 1/(1+rho + (a2^2)*((1+rho)/rho) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    
    n1=((tau1^2)*((z_a+z_b)^2)/(theta^2))* der 
    n2=r*n1
    
    n1_t = ceiling(n1 * (rho/(1+rho)))
    n1_c = ceiling(n1 * (1/(1+rho)))
    n2_t = ceiling(n2 * (rho/(1+rho)))
    n2_c = ceiling(n2 * (1/(1+rho)))
    
    n1=(n1_t+n1_c)
    n2=(n2_t+n2_c)
    
    n_val=n1+n2
  }
  
  
  #-------------------------------
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val,
                    n1_c = n1_c,
                    n1_t = n1_t,
                    n2_c = n2_c,
                    n2_t = n2_t)
  
  return(final_result)
  
  
}


# single_one_stage_size(tau=20, theta=10, rho=1, alpha=0.05, beta=0.2)
#hybrid_one_stage_size(tau1=20,tau2=25, theta=10, rho=1, r=3, alpha=0.05, beta=0.2)
#----------------------------------------------------------------------
hybrid_one_stage_dat<-function(tau1, tau1_e=0, tau2, tau2_e=0, theta, pi_IF, r, rho,
                               alpha, beta, direction, var_type){
  
  
  target=rep(target.range(H1=theta,H0=0), each=3)
  r_bias = rep(c(pi_IF-0.1, pi_IF, pi_IF+0.1),3)
  target_r = target*(1+r_bias)
  nround=max(decimalplaces(theta),decimalplaces(tau1))+1
  dat<- data.frame(rbind(round(target,nround),round(target_r,nround),r_bias, sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist(hybrid_one_stage_size(tau1, tau1_e, tau2, tau2_e, target[x],r_bias[x], r, rho,
                                   alpha, beta, direction, var_type))
      
    ))))
  
  dat = dat[1:6,]
  
  
  dat_standard<- data.frame(rbind(round(target,nround), sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist(single_one_stage_size(tau=tau1, tau_e=tau1_e, target[x], rho, alpha, beta,
                                   direction, var_type))
      
    ))))
  
  #dat = rbind(dat, dat_standard[4,])
  
  
  rownames(dat)=c("Treatment Effect \u03B4\u2081 (Onsite)", "Treatment Effect \u03B4\u2082 (Offsite)",
                  "Relative Bias (\u03B4\u2082/\u03B4\u2081-1)"  ,  "Sample Size (Onsite)",
                  "Sample Size (Offsite)", "Total Sample Size")
  dat[c(4,5,6),dat[1,]<0]=NA  #difference in mean suppose to > 0
  return(dat)
}

#--------------------------------------------------------------------------------------



##################################################################################
#Correlated data
###################################################################################


cd_single_one_stage_size=function( tau, tau_e=0,theta,lambda, rho,n_s,
                                   alpha, beta, direction, var_type){
  
  
  
  if(direction==1 & var_type==1) {
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    n_cluster =((lambda+1/lambda+2)*((tau^2)*(z_a+z_b)^2)*(1+ (n_s-1)*rho))/(n_s*(theta^2))
    n1_cluster=(1/(1+lambda))*n_cluster
    n2_cluster=(lambda/(1+lambda))*n_cluster

    n1=ceiling(n1_cluster)
    n2=ceiling(n2_cluster)
    
    n_val=n1+n2
  }
  
  
  if(direction==2 & var_type==1) {
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    n_cluster =((lambda+1/lambda+2)*((tau^2)*(z_a+z_b)^2)*(1+ (n_s-1)*rho))/(n_s*(theta^2))
    n1_cluster=(1/(1+lambda))*n_cluster
    n2_cluster=(lambda/(1+lambda))*n_cluster
    
    n1=ceiling(n1_cluster)
    n2=ceiling(n2_cluster)
    
    n_val=n1+n2
  }
  
  
  if(direction==1 & var_type==2) {
    
    #----------------------------------------
    #adjusted
    
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    a=tau_e/tau
    
    n_cluster =((1+lambda+(a^2)/lambda+a^2)*((tau^2)*(z_a+z_b)^2)*(1+ (n_s-1)*rho))/(n_s*(theta^2))
    n1_cluster=(1/(1+lambda))*n_cluster
    n2_cluster=(lambda/(1+lambda))*n_cluster
    
    n1=ceiling(n1_cluster)
    n2=ceiling(n2_cluster)
    
    n_val=n1+n2
    
  }
  
  
  if(direction==2 & var_type==2) {
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    a=tau_e/tau
    
    n_cluster =((1+lambda+(a^2)/lambda+a^2)*((tau^2)*(z_a+z_b)^2)*(1+ (n_s-1)*rho))/(n_s*(theta^2))
    n1_cluster=(1/(1+lambda))*n_cluster
    n2_cluster=(lambda/(1+lambda))*n_cluster
    
    n1=ceiling(n1_cluster)
    n2=ceiling(n2_cluster)
    
    n_val=n1+n2
  }
  
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val
  )
  
  return(final_result)
  
}


#--------------------------------------------------------------------------------------
cd_single_one_stage_dat=function( tau, tau_e=0,theta,lambda, rho,n_s,
                                   alpha, beta, direction, var_type){
  
  target=target.range.lg(H1=theta,H0=0) 
  nround=max(decimalplaces(theta),decimalplaces(tau))+1
  dat<- data.frame(rbind(round(target,nround), sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist( cd_single_one_stage_size( tau, tau_e,target[x], lambda, rho,n_s,
                                        alpha, beta, direction, var_type))
    ))))
  
  
  rownames(dat)=c("Treatment Effect (Offsite)", "Sample Size (Control Arm)",
                  "Sample Size (Experimental Arm)", "Total sample size")
  
  
  dat[c(2,3,4),dat[1,]<0]=NA  #difference in mean suppose to > 0
  return(dat)

}




#--------------------------------------------------------------------------------------
cd_hybrid_one_stage_size=function(tau1,tau1_e=0, tau2, tau2_e=0,theta, pi_IF, r, lambda, rho1, rho2, n_s1, n_s2,
                                 alpha, beta, direction, var_type){
  
  
  if(direction==1 & var_type==1) {
    
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    v=lambda+1/lambda+2
    nu1 = sqrt((n_s1-1)*rho1+1)
    nu2 = sqrt((n_s2-1)*rho2+1)
    del= ((1+pi_IF)^2)*r/(k)
    
    
    
    n1_cluster=(((tau1^2)*((z_a+z_b)^2)*v)/(theta^2)) * ( 1/(((n_s1/(nu1^2)) + (del*n_s2/(nu2^2)))) )
    n2_cluster=n1_cluster*r
    n1_cluster_t = ceiling(n1_cluster * (lambda/(1+lambda)))
    n1_cluster_c = ceiling(n1_cluster * (1/(1+lambda)))
    n2_cluster_t = ceiling(n2_cluster * (lambda/(1+lambda)))
    n2_cluster_c = ceiling(n2_cluster * (1/(1+lambda)))
    n1=(n1_cluster_t+n1_cluster_c)
    n2=(n2_cluster_t+n2_cluster_c)
    
    n_val=n1+n2
    
  }
  
  
  if(direction==2 & var_type==1) {
    
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    v=lambda+1/lambda+2
    nu1 = sqrt((n_s1-1)*rho1+1)
    nu2 = sqrt((n_s2-1)*rho2+1)
    del= ((1+pi_IF)^2)*r/(k)
    
    
    n1_cluster=(((tau1^2)*((z_a+z_b)^2)*v)/(theta^2)) * ( 1/(((n_s1/(nu1^2)) + (del*n_s2/(nu2^2)))) )
    n2_cluster=n1_cluster*r
    n1_cluster_t = ceiling(n1_cluster * (lambda/(1+lambda)))
    n1_cluster_c = ceiling(n1_cluster * (1/(1+lambda)))
    n2_cluster_t = ceiling(n2_cluster * (lambda/(1+lambda)))
    n2_cluster_c = ceiling(n2_cluster * (1/(1+lambda)))
    n1=(n1_cluster_t+n1_cluster_c)
    n2=(n2_cluster_t+n2_cluster_c)
    
    n_val=n1+n2
  }
  
  
  if(direction==1 & var_type==2) {
    
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    #adjusted
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(lambda+1)
    nu1 = sqrt((n_s1-1)*rho1+1)
    nu2 = sqrt((n_s2-1)*rho2+1)
    
    del= ((1+pi_IF)^2)*r/(k)
    
    parta = n_s1/(   (1+lambda)*(nu1^2) + ((1+lambda)/lambda)*(a1^2)*(nu1^2) )
    partb= n_s2/(   (1+lambda)*(nu2^2) + ((1+lambda)/lambda)*(a2^2)*(nu2^2) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    n1_cluster=(((tau1^2)*((z_a+z_b)^2))/(theta^2)) * der 
    
    
    
    n2_cluster=n1_cluster*r
    n1_cluster_t = ceiling(n1_cluster * (lambda/(1+lambda)))
    n1_cluster_c = ceiling(n1_cluster * (1/(1+lambda)))
    n2_cluster_t = ceiling(n2_cluster * (lambda/(1+lambda)))
    n2_cluster_c = ceiling(n2_cluster * (1/(1+lambda)))
    n1=(n1_cluster_t+n1_cluster_c)
    n2=(n2_cluster_t+n2_cluster_c)
    
    n_val=n1+n2
    
  }
  
  
  if(direction==2 & var_type==2) {
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    #adjusted
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(lambda+1)
    nu1 = sqrt((n_s1-1)*rho1+1)
    nu2 = sqrt((n_s2-1)*rho2+1)
    
    del= ((1+pi_IF)^2)*r/(k)
    
    parta = n_s1/(   (1+lambda)*(nu1^2) + ((1+lambda)/lambda)*(a1^2)*(nu1^2) )
    partb= n_s2/(   (1+lambda)*(nu2^2) + ((1+lambda)/lambda)*(a2^2)*(nu2^2) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    n1_cluster=(((tau1^2)*((z_a+z_b)^2))/(theta^2)) * der
    
    n2_cluster=n1_cluster*r
    n1_cluster_t = ceiling(n1_cluster * (lambda/(1+lambda)))
    n1_cluster_c = ceiling(n1_cluster * (1/(1+lambda)))
    n2_cluster_t = ceiling(n2_cluster * (lambda/(1+lambda)))
    n2_cluster_c = ceiling(n2_cluster * (1/(1+lambda)))
    n1=(n1_cluster_t+n1_cluster_c)
    n2=(n2_cluster_t+n2_cluster_c)
    
    n_val=n1+n2
    
  }
  
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val,
                    n1_c = n1_cluster_c,
                    n1_t = n1_cluster_t,
                    n2_c = n2_cluster_c,
                    n2_t = n2_cluster_t
  )
  
  return(final_result)

  
}



#--------------------------------------------------------------------------------------
cd_hybrid_one_stage_dat=function( tau1,tau1_e=0, tau2, tau2_e=0,theta,pi_IF,  r, lambda, rho1, rho2, n_s1, n_s2,
                                  alpha, beta, direction, var_type){
  
  
  target=rep(target.range(H1=theta,H0=0), each=3)
  r_bias = rep(c(pi_IF-0.1, pi_IF, pi_IF+0.1),3)
  target_r = target*(1+r_bias)
  nround=max(decimalplaces(theta),decimalplaces(tau1))+1
  dat<- data.frame(rbind(round(target,nround),round(target_r,nround), r_bias, sapply(1:length(target), function(x) 
    ceiling(
      
      unlist( cd_hybrid_one_stage_size(tau1,tau1_e, tau2, tau2_e,target[x], r_bias[x], r, lambda, rho1, rho2, n_s1, n_s2,
                                       alpha, beta, direction, var_type))
    ))))
  
  dat = dat[1:6,]
  
  dat_standard<- data.frame(rbind(round(target,nround), sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist( cd_single_one_stage_size( tau=tau1, tau_e=tau1_e,target[x], lambda, rho=rho1,n_s=n_s1,
                                        alpha, beta, direction, var_type))
    ))))
  
  #dat = rbind(dat, dat_standard[4,])
  
  
  rownames(dat)=c("Treatment Effect \u03B4\u2081 (Onsite)", "Treatment Effect \u03B4\u2082 (Offsite)",
                  "Relative Bias (\u03B4\u2082/\u03B4\u2081-1)"  ,  "Sample Size (Onsite)",
                  "Sample Size (Offsite)", "Total Sample Size")
  dat[c(4,5,6),dat[1,]<0]=NA  #difference in mean suppose to > 0
  
  return(dat)
  
}


##################################################################################
##################################################################################
##################################################################################

#Binary data:

##################################################################################


##################################################################################
# Completely randomized data
###################################################################################

BIN_single_one_stage_size=function(p,p_e, rho, alpha, beta, direction){
  
  theta = p_e - p

  tau = sqrt(p*(1-p))
  tau_e = sqrt(p_e*(1-p_e))
  
  
  if(direction==1 ) {
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    a=tau_e/tau
    n1=((1+((a^2)/rho) )*((tau^2)*(z_a+z_b)^2))/theta^2   
    n2=rho*n1
    n1=2*ceiling(n1/2)
    n2=2*ceiling(n2/2)
    n=n1+n2
    n_val=n1+n2
  }
  
  
  if(direction==2 ) {
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    a=tau_e/tau
    n1=((1+((a^2)/rho) )*((tau^2)*(z_a+z_b)^2))/theta^2   
    n2=rho*n1
    n1=2*ceiling(n1/2)
    n2=2*ceiling(n2/2)
    n=n1+n2
    n_val=n1+n2
  }
  
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val
  )
  
  return(final_result)
}


#------------------------------------------------------------------------------

BIN_single_one_stage_dat<-function(p, p_e, rho, alpha, beta, direction){
  
  theta = p_e - p
  tau = sqrt(p*(1-p))
  tau_e = sqrt(p_e*(1-p_e))
  
  target=target.range.lg(H1=theta,H0=0)
  nround=max(decimalplaces(theta),decimalplaces(tau))+1
  dat<- data.frame(rbind(round(target,nround), sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist(BIN_single_one_stage_size(p, p_e=p+target[x], rho, alpha, beta,
                                   direction))
      
    ))))
  
  
  rownames(dat)=c("Treatment Effect (Offsite)", "Sample Size (Control Arm)",
                  "Sample Size (Experimental Arm)", "Total sample size")
  dat[c(2,3,4),dat[1,]<0]=NA  #difference in mean suppose to > 0
  return(dat)
}


################################################################################


BIN_hybrid_one_stage_size=function(p1, p1_e, p2, p2_e, r, rho,
                               alpha, beta, direction){
  
  theta = p1_e - p1
  pi_IF = (p2_e-p2)/(p1_e - p1) - 1
  
  tau1 = sqrt(p1*(1-p1))
  tau1_e = sqrt(p1_e*(1-p1_e))
  tau2 = sqrt(p2*(1-p2))
  tau2_e = sqrt(p2_e*(1-p2_e))
  
  
  if(direction==1) {
    #------------------------------------
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(rho+1)
    
    parta = 1/(1+rho + (a1^2)*((1+rho)/rho) )
    partb= 1/(1+rho + (a2^2)*((1+rho)/rho) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    
    n1=((tau1^2)*((z_a+z_b)^2)/(theta^2))* der
    #n1
    n2=r*n1
    #n2
    n1_t = ceiling(n1 * (rho/(1+rho)))
    n1_c = ceiling(n1 * (1/(1+rho)))
    n2_t = ceiling(n2 * (rho/(1+rho)))
    n2_c = ceiling(n2 * (1/(1+rho)))
    
    
    n1=(n1_t+n1_c)
    n2=(n2_t+n2_c)
    n_val=n1+n2
  }
  
  
  if(direction==2) {
    #------------------------------------
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(rho+1)
    
    parta = 1/(1+rho + (a1^2)*((1+rho)/rho) )
    partb= 1/(1+rho + (a2^2)*((1+rho)/rho) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    
    n1=((tau1^2)*((z_a+z_b)^2)/(theta^2))* der 
    n2=r*n1
    
    n1_t = ceiling(n1 * (rho/(1+rho)))
    n1_c = ceiling(n1 * (1/(1+rho)))
    n2_t = ceiling(n2 * (rho/(1+rho)))
    n2_c = ceiling(n2 * (1/(1+rho)))
    
    n1=(n1_t+n1_c)
    n2=(n2_t+n2_c)
    
    n_val=n1+n2
  }
  
  
  #-------------------------------
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val,
                    n1_c = n1_c,
                    n1_t = n1_t,
                    n2_c = n2_c,
                    n2_t = n2_t)
  
  return(final_result)
  
  
}



#----------------------------------------------------------------------
BIN_hybrid_one_stage_dat<-function(p1, p1_e, p2, p2_e, r, rho,
                               alpha, beta, direction){
  
  theta = p1_e - p1
  pi_IF = (p2_e-p2)/(p1_e - p1) - 1
  
  tau1 = sqrt(p1*(1-p1))
  tau1_e = sqrt(p1_e*(1-p1_e))
  tau2 = sqrt(p2*(1-p2))
  tau2_e = sqrt(p2_e*(1-p2_e))
  
  target=rep(target.range(H1=theta,H0=0), each=3)
  r_bias = rep(c(pi_IF-0.1, pi_IF, pi_IF+0.1),3)
  target_r = target*(1+r_bias)
  nround=max(decimalplaces(theta),decimalplaces(tau1))+1
  dat<- data.frame(rbind(round(target,nround),round(target_r,nround),round(r_bias,nround), sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist(BIN_hybrid_one_stage_size(p1, p1_e=p1+target[x], p2, p2_e = p2+target_r[x], r, rho,
                                   alpha, beta, direction))
      
    ))))
  
  dat = dat[1:6,]
  
  
  rownames(dat)=c("Treatment Effect \u03B4\u2081 (Onsite)", "Treatment Effect \u03B4\u2082 (Offsite)",
                  "Relative Bias (\u03B4\u2082/\u03B4\u2081-1)"  ,  "Sample Size (Onsite)",
                  "Sample Size (Offsite)", "Total Sample Size")
  dat[c(4,5,6),dat[1,]<0]=NA  #difference in mean suppose to > 0
  return(dat)
}

#--------------------------------------------------------------------------------------



##################################################################################
#Correlated data
###################################################################################


BIN_cd_single_one_stage_size=function( p, p_e, lambda, rho,n_s,
                                   alpha, beta, direction){
  
  
  
  theta = p_e - p
  
  tau = sqrt(p*(1-p))
  tau_e = sqrt(p_e*(1-p_e))

  
  

  if(direction==1 ) {
    
    #----------------------------------------
    #adjusted
    
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    a=tau_e/tau
    
    n_cluster =((1+lambda+(a^2)/lambda+a^2)*((tau^2)*(z_a+z_b)^2)*(1+ (n_s-1)*rho))/(n_s*(theta^2))
    n1_cluster=(1/(1+lambda))*n_cluster
    n2_cluster=(lambda/(1+lambda))*n_cluster
    
    n1=ceiling(n1_cluster)
    n2=ceiling(n2_cluster)
    
    n_val=n1+n2
    
  }
  
  
  if(direction==2) {
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    a=tau_e/tau
    
    n_cluster =((1+lambda+(a^2)/lambda+a^2)*((tau^2)*(z_a+z_b)^2)*(1+ (n_s-1)*rho))/(n_s*(theta^2))
    n1_cluster=(1/(1+lambda))*n_cluster
    n2_cluster=(lambda/(1+lambda))*n_cluster
    
    n1=ceiling(n1_cluster)
    n2=ceiling(n2_cluster)
    
    n_val=n1+n2
  }
  
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val
  )
  
  return(final_result)
  
}


#--------------------------------------------------------------------------------------
BIN_cd_single_one_stage_dat=function( p, p_e, lambda, rho,n_s,
                                  alpha, beta, direction){
  
  target=target.range.lg(H1=theta,H0=0) 
  nround=max(decimalplaces(theta),decimalplaces(tau))+1
  dat<- data.frame(rbind(round(target,nround), sapply(1:length(target), function(x) 
    
    ceiling(
      
      unlist( BIN_cd_single_one_stage_size( p, p_e=p+target[x], lambda, rho,n_s,
                                        alpha, beta, direction))
    ))))
  
  
  rownames(dat)=c("Treatment Effect (Offsite)", "Sample Size (Control Arm)",
                  "Sample Size (Experimental Arm)", "Total sample size")
  
  
  dat[c(2,3,4),dat[1,]<0]=NA  #difference in mean suppose to > 0
  return(dat)
  
}




#--------------------------------------------------------------------------------------
BIN_cd_hybrid_one_stage_size=function(p1,p1_e, p2, p2_e, r, lambda, rho1, rho2, n_s1, n_s2,
                                  alpha, beta, direction){
  

  
  theta = p1_e - p1
  pi_IF = (p2_e-p2)/(p1_e - p1) - 1
  
  tau1 = sqrt(p1*(1-p1))
  tau1_e = sqrt(p1_e*(1-p1_e))
  tau2 = sqrt(p2*(1-p2))
  tau2_e = sqrt(p2_e*(1-p2_e))
  
  if(direction==1 ) {
    
    z_a=qnorm(1-alpha)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    #adjusted
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(lambda+1)
    nu1 = sqrt((n_s1-1)*rho1+1)
    nu2 = sqrt((n_s2-1)*rho2+1)
    
    del= ((1+pi_IF)^2)*r/(k)
    
    parta = n_s1/(   (1+lambda)*(nu1^2) + ((1+lambda)/lambda)*(a1^2)*(nu1^2) )
    partb= n_s2/(   (1+lambda)*(nu2^2) + ((1+lambda)/lambda)*(a2^2)*(nu2^2) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    n1_cluster=(((tau1^2)*((z_a+z_b)^2))/(theta^2)) * der 
    
    
    
    n2_cluster=n1_cluster*r
    n1_cluster_t = ceiling(n1_cluster * (lambda/(1+lambda)))
    n1_cluster_c = ceiling(n1_cluster * (1/(1+lambda)))
    n2_cluster_t = ceiling(n2_cluster * (lambda/(1+lambda)))
    n2_cluster_c = ceiling(n2_cluster * (1/(1+lambda)))
    n1=(n1_cluster_t+n1_cluster_c)
    n2=(n2_cluster_t+n2_cluster_c)
    
    n_val=n1+n2
    
  }
  
  
  if(direction==2 ) {
    z_a=qnorm(1-alpha/2)
    #z_a
    z_b=qnorm(1-beta)
    #z_b
    
    #adjusted
    k=(tau2/tau1)^2
    #k2=tau2_e/tau1_e
    a1=tau1_e/tau1
    a2=tau2_e/tau2
    b=1/(lambda+1)
    nu1 = sqrt((n_s1-1)*rho1+1)
    nu2 = sqrt((n_s2-1)*rho2+1)
    
    del= ((1+pi_IF)^2)*r/(k)
    
    parta = n_s1/(   (1+lambda)*(nu1^2) + ((1+lambda)/lambda)*(a1^2)*(nu1^2) )
    partb= n_s2/(   (1+lambda)*(nu2^2) + ((1+lambda)/lambda)*(a2^2)*(nu2^2) )
    partc = (((1+pi_IF)^2)*r/k)
    der=1/(parta + partc*partb)
    
    n1_cluster=(((tau1^2)*((z_a+z_b)^2))/(theta^2)) * der
    
    n2_cluster=n1_cluster*r
    n1_cluster_t = ceiling(n1_cluster * (lambda/(1+lambda)))
    n1_cluster_c = ceiling(n1_cluster * (1/(1+lambda)))
    n2_cluster_t = ceiling(n2_cluster * (lambda/(1+lambda)))
    n2_cluster_c = ceiling(n2_cluster * (1/(1+lambda)))
    n1=(n1_cluster_t+n1_cluster_c)
    n2=(n2_cluster_t+n2_cluster_c)
    
    n_val=n1+n2
    
  }
  
  
  final_result=list(n1=n1,
                    n2=n2,
                    n=n_val,
                    n1_c = n1_cluster_c,
                    n1_t = n1_cluster_t,
                    n2_c = n2_cluster_c,
                    n2_t = n2_cluster_t
  )
  
  return(final_result)
  
  
}



#--------------------------------------------------------------------------------------
BIN_cd_hybrid_one_stage_dat=function( p1,p1_e, p2, p2_e, r, lambda, rho1, rho2, n_s1, n_s2,
                                  alpha, beta, direction){
  
  
  target=rep(target.range(H1=theta,H0=0), each=3)
  r_bias = rep(c(pi_IF-0.1, pi_IF, pi_IF+0.1),3)
  target_r = target*(1+r_bias)
  nround=max(decimalplaces(theta),decimalplaces(tau1))+1
  dat<- data.frame(rbind(round(target,nround),round(target_r,nround), round(r_bias,nround), sapply(1:length(target), function(x) 
    ceiling(
      
      unlist( BIN_cd_hybrid_one_stage_size(p1,p1_e=p+target[x], p2, p2_e=p2+target_r[x], r, lambda, rho1, rho2, n_s1, n_s2,
                                       alpha, beta, direction))
    ))))
  
  dat = dat[1:6,]
  
  
  
  rownames(dat)=c("Treatment Effect \u03B4\u2081 (Onsite)", "Treatment Effect \u03B4\u2082 (Offsite)",
                  "Relative Bias (\u03B4\u2082/\u03B4\u2081-1)"  ,  "Sample Size (Onsite)",
                  "Sample Size (Offsite)", "Total Sample Size")
  dat[c(4,5,6),dat[1,]<0]=NA  #difference in mean suppose to > 0
  
  return(dat)
  
}





