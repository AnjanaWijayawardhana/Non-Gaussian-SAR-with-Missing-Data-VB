#### Require packages
library(LaplacesDemon)
library(extraDistr)
library(igraph)
library(spdep)
library(Matrix)
library(tictoc)
library(crayon)
library(spatialreg)
library(mvtnorm)
library(numDeriv)
library(Matrix)
# library(rstan)
library(MASS)
library(ggplot2)
library(reshape2)


########### Functions for simulation





simulateSEM<-function(para,p,weightmat){
  
  
  weightmat<-as(weightmat,"CsparseMatrix")
  sigma2=para[p+3]
  rho=para[p+2]
  n=ncol(weightmat)
  m=0
  std=1
  X<-matrix(rep(c(1),n))
  
  for(i in 1:p){
    X<-cbind(X,rnorm(n,m,std))
  }
  
  b<-para[0:p+1]
  errors_part<-rnorm(n,sd=sqrt(sigma2))
  
  Y<-X%*%b+solve(Diagonal(n)-rho*weightmat,errors_part)
  
  return(list("Independent"=X[,2:(p+1)],"Dependent"=Y))
  
} # Simulates data from a Gaussian Spatial Error Model (SEM-Gau).


simulateSEMterrors<-function(para,p, # Simulates data from a Student t Spatial Error Model (SEM-t).
                             weightmat){
  
  
  weightmat<-as(weightmat,"CsparseMatrix")
  rho=para[p+2]
  sigma2=para[p+3]
  df<-para[length(para)]
  n=ncol(weightmat)
  m=0
  std=1
  X<-matrix(rep(c(1),n))
  
  for(i in 1:p){
    X<-cbind(X,rnorm(n,m,std))
  }
  
  b<-para[0:p+1]
  # df=4
  # errors_part<-rt(n,df=dfree)*sqrt(sigma2)
  # errors_part<-rstudent_t(n,dfree,0,sqrt(sigma2))
  # errors_part=gen.t.errors(n,df=dfree,sigma2)
  
  # ?rmvt
  # errors_part=LaplacesDemon::rmvt(1,rep(0,n),sigma2*diag(n),df=dfree)
  # errors_part=as.numeric(errors_part)
  errors_part=rt(n, df = df) * sqrt(sigma2)
  # errors_part=r_fastmvt(1,dfree,rep(0,n),sigma2*Diagonal(n))
  # errors_part=as.numeric(errors_part)
  
  Y<-X%*%b+solve(Diagonal(n)-rho*weightmat,errors_part)
  
  return(list("Independent"=X[,2:(p+1)],"Dependent"=Y,"e"=errors_part))
  
} 



t_inverse<-function(z,gama){ 
  y<-c()
  # i=1
  for(i in 1:length(z)){
    if(z[i]<0){
      y[i]<-1-(1-z[i]*(2-gama))^(1/(2-gama))
    }else{
      y[i]<-(1+z[i]*gama)^(1/gama)-1
    }
  }
  return(y)
} # Computes the inverse of the Yeo-Johnson transformation.


row_standardize_to_sum_1<-function(w) {
  # Get the row sums
  row_sums <- rowSums(w)
  
  # Avoid division by zero (if there are any rows with sum 0, handle them separately)
  row_sums[row_sums == 0] <- 1  # Set zero sums to 1 to avoid division by zero
  
  # Divide each element in a row by the sum of that row
  standardized_w <- sweep(w, 1, row_sums, "/")
  
  return(standardized_w)
} # Makes  Row-standardized W 


################## SEM-Gau, SEM-t, YJ-SEM-Gau, and YJ-SEM-t algorithms

SEM_Gau<-function(x,y,N, # Fit SEM-Gau 
                  w,start.theta,p){ 
  # p=1
  x<-cbind(1,x)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=0     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-8)
  v=0.95
  
  n<-length(y)
  
  I<-Diagonal(n)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  # theta_u<-c(c(1,2),0.5,0.1,c(3,3,3))
  grad<-function(theta_u,wpluswt,wtw,w,x,y,I,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p mean dimension of X.
    
    beta<-theta_u[(1):(p)]
    gama<-theta_u[(p+1)] #transform of sigma2
    lamda<-theta_u[(p+2)] #transform of sigma rho
    
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    omega_psi=10000
    
    
    
    A=I-rho*w
    M=t(A)%*%A
    #betas
    div_h_beta<-exp(-gama)*t(y-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #gamma
    div_h_gamma<--0.5*n+(0.5*exp(-gama))*t(y-x%*%beta)%*%M%*%(y-x%*%beta)-gama/omega_gamma
    
    #Lamda
    drho_by_dlamda<-2*exp(lamda)/(1+exp(lamda))^2
    
    dM_by_dlamda<-(-wpluswt+2*rho*wtw)*drho_by_dlamda
    # dlogdetM_by_dlamda<-sum(diag(solve(M,dM_by_dlamda))) # use when real data w is used.
    dlogdetM_by_dlamda<-sum(diag(solve(Cholesky(M),dM_by_dlamda))) # use when simulated data (rook W) is used.
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  cal_Lbound<-function(theta_u,d,B,x,y,wpluswt,wtw){
    
    p<-ncol(x) #length of betas
    beta<-theta_u[1:p]
    gama<-theta_u[(p+1)]
    lamda<-theta_u[(p+2)]
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    
    n<-length(y)
    A=I-rho*w
    M=t(A)%*%A
    
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_psi=10000
    omega_lamda=10000
    D<-diag(as.vector(d))
    
    
    
    
    
    loglike1<-0 #p_m
    # loglike1<-sum(y*z-log(1+exp(z)))
    
    # sigma_<-exp(gama)*solve(M)
    # loglike2<-dmvnorm(x=as.vector(y),mean=x%*%beta,
    #                   sigma=sigma_, log=T,checkSymmetry = F) # direct p_y 
    
    log_det<-2*sum(log(diag(chol(M))))
    loglike2<--0.5*n*log(2*pi)-0.5*n*gama+
      0.5*log_det-0.5*exp(-gama)*t(y-x%*%beta)%*%M%*%(y-x%*%beta) # p_y efficient
    
    loglike3<-dmvnorm(x=as.vector(beta),mean = rep(0,ncol(x)),sigma=omega_beta*diag(ncol(x)),
                      log=T,checkSymmetry = F) #p_beta
    loglike4<-dnorm(x=gama,mean = 0,sd=omega_gamma) #p_gama (variance)
    
    loglike5<-dnorm(x=lamda,mean = 0,sd=omega_lamda) #p_lamda (variance)
    
    loglike6<-0
    
    logq<-dmvnorm(x=as.vector(theta_u),mean = as.vector(mu),sigma =B%*%t(B)+
                    D^2,log = T,checkSymmetry = F) # q
    
    
    return(as.numeric(loglike1)+as.numeric(loglike2)+
             as.numeric(loglike3)+as.numeric(loglike4)+
             as.numeric(loglike5)+as.numeric(loglike6)-as.numeric(logq))
    
    
  } #to calculate lower bound
  #initial values for lamdha
  #p=2
  mu<-start.theta
  d<-rep(0.1,totpara)
  B<-matrix(rep(0.1,(totpara)*p),ncol=p)
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara)
  E_g2_d<-rep(0,totpara)
  E_g2_B<-rep(0,totpara) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara)
  E_delta2_d<-rep(0,totpara)
  E_delta2_B<-rep(0,totpara) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara)
  rho_B<-rep(0,totpara)
  rho_d<-rep(0,totpara)
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N)
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    theta<-mu+B%*%z+d*as.vector(epsilon)
    # theta_t<-matrix(theta_t,ncol=1)
    D<-diag(d)
    
    # Calculate log gradient of theta
    
    
    gradh_theta<-grad(theta_u=theta,w=w,x=x,y=y,I=I,
                      wpluswt=wpluswt,wtw=wtw,n=n) # gratient of log h(theta)
    part2<-solve(forceSymmetric((B%*%t(B)+D^2)))%*%(B%*%z+d*as.vector(epsilon))# the second part inside the expectation
    
    
    #Construct UB estimates
    grad_mu<-gradh_theta+part2
    grad_B<-gradh_theta%*%t(z)+part2%*%t(z)
    gradd<-diag(gradh_theta%*%t(epsilon)+part2%*%t(epsilon))
    grad_d<-matrix(gradd,ncol = 1)
    
    #Set new learning rates
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B<-B+rho_B*grad_B
    d<-d+rho_d*grad_d
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    # calculate lowerbound
    
    # Lbound[i]<-cal_Lbound(theta_u=theta,d=d,B=B,
    # x=x,y=y,wpluswt=wpluswt,wtw=wtw)
    Lbound[i]<-0
    all_paras[i,]<-as.vector(mu)
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,Lbound=Lbound,all_paras=all_paras))
} 

YJ_SEM_Gau<-function(x,y,N,w, #Fit YJ-SEM-Gau
                     start.theta,p){ 
  # p=1
  x<-cbind(1,x)
  rc=ncol(x) #number of regression covariates. i.e. betas
  totpara<-rc+3
  adapt_epsilon=10^(-6)
  v=0.99
  
  n<-length(y)
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  # theta_u<-c(c(1,2),0.5,0.1,c(3,3,3))
  cal_z<-function(y,gama){
    y_negative<-function(y,gama){
      return(-((-y+1)^(2-gama)-1)/(2-gama))
    }
    y_else<-function(y,gama){
      return(((y+1)^(gama)-1)/gama)
    }
    z <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
  }
  
  cal_der_tYJ_wrt_y<-function(y,gama){
    y_negative<-function(y,gama){
      return((-y+1)^(1-gama))
    }
    y_else<-function(y,gama){
      return((y+1)^(gama-1))
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  } #grad of t_jy wrt y
  
  cal_der_tYJ_wrt_gama<-function(y,gama){
    y_negative<-function(y,gama){
      der<-((2-gama)*(-y+1)^(2-gama)*log(-y+1)-(-y+1)^(2-gama)+1)/(2-gama)^2
      return(der)
    }
    y_else<-function(y,gama){
      der<-(gama*(y+1)^gama*log(y+1)-(y+1)^gama+1)/gama^2
      return(der)
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  } #grad of t_jy wrt gama
  
  cal_second_der_wrt_y_gama<-function(y,gama){
    y_negative<-function(y,gama){
      der<--(-y+1)^(1-gama)*log(-y+1)
      return(der)
    }
    y_else<-function(y,gama){
      der<-(y+1)^(gama-1)*log(y+1)
      return(der)
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama)) # there will be warnings like  In log(-y + 1) : NaNs produced and In log(y + 1) : NaNs produced 
    return(der_tYJ)                                 # maybe all elements in y apply for both function but return values correclty
  } #second der of t_jy wrt y and gama
  
  grad<-function(theta_u,wpluswt,wtw,w,x,I,y,mu,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p mean dimension of X.
    
    beta<-theta_u[(1):(p)]
    omega.dash<-theta_u[(p+1)] #transform of sigma2
    rho.dash<-theta_u[(p+2)] #transform of sigma rho
    gama.dash<-theta_u[(p+3)]
    
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_omega.dash=10000
    omega_rho.dash=10000
    omega_gama.dash=10000
    
    A=I-rho*w
    M=t(A)%*%A
    z<-cal_z(y,gama)
    #betas
    div_h_beta<-exp(-omega.dash)*t(z-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #omega.dash
    div_h_omega.dash<--0.5*n+(0.5*exp(-omega.dash))*t(z-x%*%beta)%*%M%*%(z-x%*%beta)-omega.dash/omega_omega.dash
    
    #rho.dash
    drho_by_drho.dash<-2*exp(rho.dash)/(1+exp(rho.dash))^2
    
    dM_by_drho.dash<-(-wpluswt+2*rho*wtw)*drho_by_drho.dash
    # dlogdetM_by_drho.dash<-sum(diag(solve(M,dM_by_drho.dash))) # use when real data w is used.
    dlogdetM_by_drho.dash<-sum(diag(solve(Cholesky(M),dM_by_drho.dash)))  # use when simulated data (rook W) is used.
    
    div_h_rho.dash<-0.5*dlogdetM_by_drho.dash-0.5*exp(-omega.dash)*t(z-x%*%beta)%*%dM_by_drho.dash%*%(z-x%*%beta)-rho.dash/omega_rho.dash
    
    #gama.dash
    d_gama_by_gama.dash<-2*exp(gama.dash)/(1+exp(gama.dash))^2
    
    der_tYJ_wrt_gama<-cal_der_tYJ_wrt_gama(y,gama)
    part1<--exp(-omega.dash)*t(z-x%*%beta)%*%M%*%der_tYJ_wrt_gama*d_gama_by_gama.dash
    
    second_der_wrt_y_gama<-cal_second_der_wrt_y_gama(y,gama)
    der_tYJ_wrt_y<-cal_der_tYJ_wrt_y(y,gama)
    part2<-sum((1/der_tYJ_wrt_y)*second_der_wrt_y_gama*d_gama_by_gama.dash)
    
    div_h_gama.dash<-part1+part2-gama.dash/omega_gama.dash
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_omega.dash),as.vector(div_h_rho.dash),as.vector(div_h_gama.dash))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  cal_Lbound<-function(theta_u,mu,d,B,x,y,wpluswt,wtw){
    
    p<-ncol(x) #length of betas
    beta<-theta_u[1:p]
    omega.dash<-theta_u[(p+1)]
    rho.dash<-theta_u[(p+2)]
    gama.dash<-theta_u[(p+3)]
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    
    
    n<-length(y)
    A=I-rho*w
    M=t(A)%*%A
    
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_omega.dash=10000
    omega_rho.dash=10000
    omega_gama.dash=10000
    D<-diag(as.vector(d))
    
    
    log_det<-2*sum(log(diag(chol(M))))
    loglike2<--0.5*n*log(2*pi)-0.5*n*omega.dash+
      0.5*log_det-0.5*exp(-omega.dash)*t(y-x%*%beta)%*%M%*%(y-x%*%beta) # p_y efficient
    
    loglike3<-dmvnorm(x=as.vector(beta),mean = rep(0,ncol(x)),sigma=omega_beta*diag(ncol(x)),
                      log=T,checkSymmetry = F) #p_beta
    loglike4<-dnorm(x=omega.dash,mean = 0,sd=omega_omega.dash) #p_omega.dash (variance)
    
    loglike5<-dnorm(x=rho.dash,mean = 0,sd=omega_rho.dash) #p_rho.dash (variance)
    
    loglike6<-dnorm(x=gama.dash,mean = 0,sd=omega_gama.dash) #p_rho.dash (variance)
    
    
    logq<-dmvnorm(x=as.vector(theta_u),mean = as.vector(mu),sigma =B%*%t(B)+
                    D^2,log = T,checkSymmetry = F) # q
    
    
    return(+as.numeric(loglike2)+
             as.numeric(loglike3)+as.numeric(loglike4)+
             as.numeric(loglike5)+as.numeric(loglike6)-as.numeric(logq))
    
    
  } #to calculate lower bound
  #initial values for lamdha
  #p=2
  mu<-start.theta
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara)
  E_g2_d<-rep(0,totpara)
  E_g2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara)
  E_delta2_d<-rep(0,totpara)
  E_delta2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara)
  rho_B<-rep(0,totpara*p)
  rho_d<-rep(0,totpara)
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) #To store lamda in each iteration
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    eta<-matrix(zeta[1:p],ncol = 1)
    
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    theta<-mu+B%*%eta+d*as.vector(epsilon)
    # theta_t<-matrix(theta_t,ncol=1)
    D<-diag(d)
    
    
    # Calculate log gradient of theta
    gradiant<-grad(theta_u=theta,wpluswt,wtw,w,x,I,y,mu,n)
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%eta+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    #Construct UB estimates
    grad_mu<-gradiant-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradiant))
    grad_B<-t(t(eta)%x%I_theta)%*%(gradiant-dlogq_by_dtheta)
    gradd<-epsilon*(gradiant-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    
    #Set new learning rates
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    d<-d+rho_d*grad_d
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    # calculate lowerbound
    
    # Lbound[i]<-cal_Lbound(theta_u=theta,mu,d,B,x,y,wpluswt,wtw)
    Lbound[i]<-0
    all_paras[i,]<-as.vector(mu)
    # print(i)
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,Lbound=Lbound,all_paras=all_paras))
} 


SEM_terrors_noninfo<-function(x,y,w,p, #Fit SEM-t
                              start.theta,N,shape=4,rate=1/2){
  # p=10
  rho.a<--1
  rho.b<-1
  x<-cbind(1,x)
  numcov<-ncol(x)
  n<-length(y)
  totpara<-3+n+numcov
  adapt_epsilon=10^(-6)
  v=0.95
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  ################################### starting values
  
  start.df.dash<-start.theta[length(start.theta)]
  start.df<-exp(start.df.dash)+3
  start.lamda<-extraDistr::rinvgamma(n,start.df/2,start.df/2)
  start.lamda.dash<-log(start.lamda)
  mu<-c(start.lamda.dash,start.theta) # mu=(start_lamda,start_theta except rho) 
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  # B<-as(B, "dgCMatrix")
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara)
  E_g2_d<-rep(0,totpara)
  E_g2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara)
  E_delta2_d<-rep(0,totpara)
  E_delta2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara)
  rho_B<-rep(0,totpara*p)
  rho_d<-rep(0,totpara)
  
  
  ##########  ############## ##
  
  # ?digamma
  # f<-function(nu_dash){(log(gamma(0.5*exp(nu_dash)+1.5)))}
  # f<-function(x){(digamma(x))}
  f<-function(df.dash){(digamma(0.5*exp(df.dash)+1.5))}
  f2<-function(df.dash){(lgamma(0.5*exp(df.dash)+1.5))}
  # tic()
  # f(nu_dash)*0.5*exp(nu_dash)
  # toc()
  # tic()
  # grad(f2,nu_dash)
  # toc()
  calgrad<-function(lamda_theta,x,y,w,wpluswt,wtw,rho.a,rho.b,I,n){
    
    numcov<-ncol(x)
    lamda.dash<-lamda_theta[1:n]
    beta<-lamda_theta[(n+1):(n+numcov)]
    omega.dash<-lamda_theta[(n+numcov+1)] #transform of sigma2 
    rho.dash<-lamda_theta[(n+numcov+2)] #transform of rho 
    df.dash<-lamda_theta[length(lamda_theta)]
    
    lamda<-exp(lamda.dash)
    omega<-exp(omega.dash)
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    df<-exp(df.dash)+3
    # print(df)
    
    
    #Hyper parameters
    omega_gama=10000
    omega_df=10000
    omega_beta=10000
    omega_rho.dash=10000
    
    # Precalculations
    r<-y-x%*%beta
    A<-I-rho*w
    Ar<-A%*%r
    
    #lamda.dash
    part1=-(1/2)
    r<-y-x%*%beta
    part2=0.5*exp(-omega.dash)*(Ar)^2*exp(-lamda.dash)
    div_sem_h_lamda<-part1+part2
    
    div_plamda_dlamda<--(df/2)+(df/2)*exp(-lamda.dash)
    div_h_lamda<-div_sem_h_lamda+div_plamda_dlamda # 1/lamda is the jec adjusment
    
    #beta
    Sigma_lamda_inv=Diagonal(n,(1/lamda))
    M=t(A)%*%Sigma_lamda_inv%*%A
    div_h_beta<-exp(-omega.dash)*t(r)%*%M%*%x-t(beta)/omega_beta
    
    
    #omega.dash
    div_h_omega.dash<--0.5*n+(0.5*exp(-omega.dash))*t(Ar)%*%Sigma_lamda_inv%*%(Ar)-omega.dash/omega_gama
    
    #rho.dash
    ATA=I-rho*wpluswt+rho*rho*wtw
    drho_by_drho.dash<--((rho.a-rho.b)*exp(rho.dash))/(1+exp(rho.dash))^2
    ATA_by_rho.dash<-(-wpluswt+2*rho*wtw)*drho_by_drho.dash
    dlogdetATA_by_rho.dash<-sum(diag(solve((ATA),ATA_by_rho.dash)))
    
    # der.of.quadratic<--0.5*exp(-omega.dash)*t(r)%*%(-Sigma_lamda_inv%*%w-t(w)%*%Sigma_lamda_inv+2*rho*t(w)%*%Sigma_lamda_inv%*%w)%*%r*drho_by_drho.dash
    der.of.quadratic<-exp(-omega.dash)*t(Ar)%*%(Sigma_lamda_inv)%*%w%*%r*drho_by_drho.dash # both expressions are correct 
    
    div_h_drho.dash<-0.5*dlogdetATA_by_rho.dash+der.of.quadratic-rho.dash/omega_rho.dash
    
    
    
    ############################
    
    #df.dash
    #1 flat prior on df.dash
    df.dash.prior=-df.dash/omega_df
    #2 Informative prior on df~gamma(alp,bet)
    # alp<-shape
    # lamb<-rate
    # df.dash.prior<-(((alp-1)*exp(df.dash))/(exp(df.dash)+3))-lamb*exp(df.dash)+1
    # 
    div_h_ddf_dash<-0.5*n*exp(df.dash)+0.5*n*log(1.5+0.5*exp(df.dash))*exp(df.dash)-
      n*f(df.dash)*0.5*exp(df.dash)-0.5*exp(df.dash)*sum(lamda.dash+exp(-lamda.dash))+df.dash.prior
    
    fullgrad<-c(as.vector(div_h_lamda),as.vector(div_h_beta),
                as.vector(div_h_omega.dash),as.vector(div_h_drho.dash),as.vector(div_h_ddf_dash))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  
  
  # i=1 N=10
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N)
  trace.xi<-c()
  no.of.accepts<-0
  for(i in 1:N){ 
    # cat(red(paste("Iteration number is:", i)), "\n") #print iteration #
    
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    lamda_theta.dash<-as.vector(mu+B%*%z+d*as.vector(epsilon))
    
    
    # lamda.dash<-lamda_theta[1:n]
    # # lamda=exp(lamda.dash)
    # omega.dash=lamda_theta[n+1]
    # df_dash=lamda_theta[n+2]
    # 
    # lamda_theta=c(lamda.dash,omega.dash,df_dash)
    # tic() lamda_theta[-(1:n)]
    
    gradiant=calgrad(lamda_theta=lamda_theta.dash,x,y,w,wpluswt,wtw,rho.a,rho.b,I,n)
    # toc()
    # print(gradiant[-(1:(n-2))])
    # print("solve")
    # tic()
    
    ##
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    ###
    
    
    ###
    
    # toc()
    #Construct UB estimates
    grad_mu<-gradiant-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradiant))
    grad_B<-t(t(z)%x%I_theta)%*%(gradiant-dlogq_by_dtheta)
    gradd<-epsilon*(gradiant-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    
    #Set new learning rates
    
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    d<-d+rho_d*grad_d
    
    
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    
    
    all_paras[i,]<-as.vector(mu)
    
    
  }
  
  return(list(mu=mu,B=B,d=d,trace.xi=trace.xi,
              all_paras=all_paras,
              no.of.accepts=no.of.accepts))
  
}



YJ_SEM_terrors_noninfo<-function(x,y,w,p, #Fit YJ-SEM-Gau
                                 start.theta,N,shape=4,rate=1/2){ 
  # p=4
  rho.a<--1
  rho.b<-1
  x<-cbind(1,x)
  rc=ncol(x) #number of regression covariates. i.e. betas
  n<-length(y)
  totpara<-n+rc+4
  adapt_epsilon=10^(-6)
  v=0.99
  
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  # theta_u<-c(c(1,2),0.5,0.1,c(3,3,3))
  cal_z<-function(y,gama){
    y_negative<-function(y,gama){
      return(-((-y+1)^(2-gama)-1)/(2-gama))
    }
    y_else<-function(y,gama){
      return(((y+1)^(gama)-1)/gama)
    }
    z <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
  } #YJ transformation
  
  cal_der_tYJ_wrt_y<-function(y,gama){
    y_negative<-function(y,gama){
      return((-y+1)^(1-gama))
    }
    y_else<-function(y,gama){
      return((y+1)^(gama-1))
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  } #grad of t_jy wrt y
  
  cal_der_tYJ_wrt_gama<-function(y,gama){
    y_negative<-function(y,gama){
      der<-((2-gama)*(-y+1)^(2-gama)*log(-y+1)-(-y+1)^(2-gama)+1)/(2-gama)^2
      return(der)
    }
    y_else<-function(y,gama){
      der<-(gama*(y+1)^gama*log(y+1)-(y+1)^gama+1)/gama^2
      return(der)
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  } #grad of t_jy wrt gama
  
  cal_second_der_wrt_y_gama<-function(y,gama){
    y_negative<-function(y,gama){
      der<--(-y+1)^(1-gama)*log(-y+1)
      return(der)
    }
    y_else<-function(y,gama){
      der<-(y+1)^(gama-1)*log(y+1)
      return(der)
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama)) # there will be warnings like  In log(-y + 1) : NaNs produced and In log(y + 1) : NaNs produced 
    return(der_tYJ)                                 # maybe all elements in y apply for both function but return values correclty
  } #second der of t_jy wrt y and gama
  
  f<-function(df.dash){(digamma(0.5*exp(df.dash)+1.5))}
  f2<-function(df.dash){(lgamma(0.5*exp(df.dash)+1.5))}
  
  calgrad<-function(theta_u,wpluswt,wtw,w,x,I,y,mu,n){
    # m<-ncol(x)
    numcov<-ncol(x) # Here p mean dimension of X.
    
    lamda.dash<-theta_u[1:n]
    beta<-theta_u[(n+1):(n+numcov)]
    omega.dash<-theta_u[(n+numcov+1)] #transform of sigma2
    rho.dash<-theta_u[(n+numcov+2)] #transform of sigma rho
    df.dash<-theta_u[(n+numcov+3)]
    gama.dash<-theta_u[length(theta_u)]
    
    lamda<-exp(lamda.dash)
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    df<-exp(df.dash)+3
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_omega.dash=10000
    omega_rho.dash=10000
    omega_gama.dash=10000
    omega_rho.dash=10000
    omega_df=10000
    
    A=I-rho*w
    # M=t(A)%*%A
    z<-cal_z(y,gama)
    Sigma_lamda_inv=Diagonal(n,(1/lamda))
    
    r<-z-x%*%beta
    A<-I-rho*w
    Ar<-A%*%r
    
    #lamda.dash
    part1=-(1/2)
    # r<-y-x%*%beta
    part2=0.5*exp(-omega.dash)*(Ar)^2*exp(-lamda.dash)
    div_sem_h_lamda<-part1+part2
    
    div_plamda_dlamda<--(df/2)+(df/2)*exp(-lamda.dash)
    div_h_lamda<-div_sem_h_lamda+div_plamda_dlamda # 1/lamda is the jec adjusment
    
    
    #betas
    div_h_beta<-exp(-omega.dash)*t(z-x%*%beta)%*%t(A)%*%Sigma_lamda_inv%*%A%*%x-t(beta)/omega_beta
    
    
    #omega.dash
    div_h_omega.dash<--0.5*n+(0.5*exp(-omega.dash))*t(Ar)%*%Sigma_lamda_inv%*%(Ar)-omega.dash/omega_omega.dash
    
    #rho.dash
    ATA=I-rho*wpluswt+rho*rho*wtw
    drho_by_drho.dash<--((rho.a-rho.b)*exp(rho.dash))/(1+exp(rho.dash))^2
    ATA_by_rho.dash<-(-wpluswt+2*rho*wtw)*drho_by_drho.dash
    dlogdetATA_by_rho.dash<-sum(diag(solve((ATA),ATA_by_rho.dash)))
    
    der.of.quadratic<--0.5*exp(-omega.dash)*t(r)%*%(-Sigma_lamda_inv%*%w-t(w)%*%Sigma_lamda_inv+2*rho*t(w)%*%Sigma_lamda_inv%*%w)%*%r*drho_by_drho.dash
    
    # der.of.quadratic_<-exp(-omega.dash)*t(Ar)%*%Sigma_lamda_inv%*%w%*%r*drho_by_drho.dash
    
    div_h_drho.dash<-0.5*dlogdetATA_by_rho.dash+der.of.quadratic-rho.dash/omega_rho.dash
    
    
    #df.dash
    #1 flat prior on df.dash
    df.dash.prior=-df.dash/omega_df
    #2 Informative prior on df~gamma(alp,bet)
    # alp<-shape
    # lamb<-rate
    # df.dash.prior<-(((alp-1)*exp(df.dash))/(exp(df.dash)+3))-lamb*exp(df.dash)+1
    
    div_h_ddf_dash<-0.5*n*exp(df.dash)+0.5*n*log(1.5+0.5*exp(df.dash))*exp(df.dash)-
      n*f(df.dash)*0.5*exp(df.dash)-0.5*exp(df.dash)*sum(lamda.dash+exp(-lamda.dash))+df.dash.prior
    
    
    #gama.dash
    d_gama_by_gama.dash<-2*exp(gama.dash)/(1+exp(gama.dash))^2
    
    der_tYJ_wrt_gama<-cal_der_tYJ_wrt_gama(y,gama)
    part1<--exp(-omega.dash)*t(Ar)%*%Sigma_lamda_inv%*%A%*%der_tYJ_wrt_gama*d_gama_by_gama.dash
    
    second_der_wrt_y_gama<-cal_second_der_wrt_y_gama(y,gama)
    der_tYJ_wrt_y<-cal_der_tYJ_wrt_y(y,gama)
    part2<-sum((1/der_tYJ_wrt_y)*second_der_wrt_y_gama*d_gama_by_gama.dash)
    
    div_h_gama.dash<-part1+part2-gama.dash/omega_gama.dash
    
    
    fullgrad<-c(as.vector(div_h_lamda),as.vector(div_h_beta),
                as.vector(div_h_omega.dash),as.vector(div_h_drho.dash),as.vector(div_h_ddf_dash),as.vector(div_h_gama.dash))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  
  #initial values for lamdha
  #p=2
  start.df.dash<-start.theta[length(start.theta)-1]
  start.df<-exp(start.df.dash)+3
  start.lamda<-extraDistr::rinvgamma(n,start.df/2,start.df/2)
  start.lamda.dash<-log(start.lamda)
  mu<-c(start.lamda.dash,start.theta)
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara)
  E_g2_d<-rep(0,totpara)
  E_g2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara)
  E_delta2_d<-rep(0,totpara)
  E_delta2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara)
  rho_B<-rep(0,totpara*p)
  rho_d<-rep(0,totpara)
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) #To store lamda in each iteration
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    eta<-matrix(zeta[1:p],ncol = 1)
    
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    theta<-mu+B%*%eta+d*as.vector(epsilon)
    # theta_t<-matrix(theta_t,ncol=1)
    D<-diag(d)
    
    # Calculate log gradient of theta
    gradiant<-calgrad(theta_u=theta,wpluswt,wtw,w,x,I,y,mu,n)
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%eta+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    #Construct UB estimates
    grad_mu<-gradiant-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradiant))
    grad_B<-t(t(eta)%x%I_theta)%*%(gradiant-dlogq_by_dtheta)
    gradd<-epsilon*(gradiant-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    
    #Set new learning rates
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    d<-d+rho_d*grad_d
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    # calculate lowerbound
    
    # Lbound[i]<-cal_Lbound(theta_u=theta,mu,d,B,x,y,wpluswt,wtw)
    Lbound[i]<-0
    all_paras[i,]<-as.vector(mu)
    # print(i)
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,Lbound=Lbound,all_paras=all_paras))
} # New code


################# Function for draw from the posterior (variational) distribution

simulate_theta_SEM_Gau<-function(mu,B,x,d){
  
  
  p<-ncol(B)
  m<-length(d)
  
  epsilon<-rnorm(m)
  z<-rnorm(p)
  
  theta.dash<-mu+B%*%z+as.vector(d)*epsilon # parameters in transformed space
  # theta<-theta.dash
  #rho
  rho.dash<-theta.dash[ncol(x)+1+1+1]
  rho<-(exp(rho.dash) -
          1)/(exp(rho.dash) +1)
  
  #sigma2
  sigma2<-exp(theta.dash[ncol(x)+1+1])
  
  #
  
  
  theta<-theta.dash
  theta[ncol(x)+1+1]<-sigma2 # parameters in actual space transformed space
  theta[ncol(x)+1+1+1]<-rho
  return(list(theta.dash=as.vector(theta.dash),theta=as.vector(theta)))
  
  
  
  
} # SEM-Gau 


simulate_theta_YJ_SEM_Gau<-function(mu,B,x,d){
  
  
  p<-ncol(B)
  m<-length(d)
  n<-nrow(w)
  
  epsilon<-rnorm(m)
  z<-rnorm(p)
  
  theta.dash<-mu+B%*%z+as.vector(d)*epsilon # parameters in transformed space
  # theta<-theta.dash
  #rho
  rho.dash<-theta.dash[ncol(x)+1+1+1]
  rho<-(exp(rho.dash) -
          1)/(exp(rho.dash) +1)
  
  #sigma2
  sigma2<-exp(theta.dash[ncol(x)+1+1])
  
  #gamma.dash
  gamma.dash<-theta.dash[ncol(x)+1+1+1+1]
  
  
  theta<-theta.dash
  theta[ncol(x)+1+1]<-sigma2 # parameters in actual space transformed space
  theta[ncol(x)+1+1+1]<-rho
  theta[ncol(x)+1+1+1+1]<-2*exp(gamma.dash)/(exp(gamma.dash)+1)
  # theta[n+ncol(x)+1+1+1+1+1]<-2*exp(gamma.dash)/(exp(gamma.dash)+1)
  return(list(theta.dash=as.vector(theta.dash),theta=as.vector(theta)))
  
  
  
  
} # YJ-SEM-Gau


simulate_theta_SEM_t<-function(mu,B,x,d){
  
  
  p<-ncol(B)
  m<-length(d)
  n<-nrow(w)
  
  epsilon<-rnorm(m)
  z<-rnorm(p)
  
  theta.dash<-mu+B%*%z+as.vector(d)*epsilon # parameters in transformed space
  # theta<-theta.dash
  #rho
  rho.dash<-theta.dash[n+ncol(x)+1+1+1]
  rho<-(exp(rho.dash) -
          1)/(exp(rho.dash) +1)
  
  #sigma2
  sigma2<-exp(theta.dash[n+ncol(x)+1+1])
  
  #
  df.dash<-theta.dash[n+ncol(x)+1+1+1+1]
  
  theta<-theta.dash
  theta[n+ncol(x)+1+1]<-sigma2 # parameters in actual space transformed space
  theta[n+ncol(x)+1+1+1]<-rho
  theta[n+ncol(x)+1+1+1+1]<-exp(df.dash)+3
  return(list(theta.dash=as.vector(theta.dash),theta=as.vector(theta)))
  
  
  
  
} # SEM-t

simulate_theta_YJ_SEM_t<-function(mu,B,x,d){
  
  
  p<-ncol(B)
  m<-length(d)
  n<-nrow(w)
  
  epsilon<-rnorm(m)
  z<-rnorm(p)
  
  theta.dash<-mu+B%*%z+as.vector(d)*epsilon # parameters in transformed space
  # theta<-theta.dash
  #rho
  rho.dash<-theta.dash[n+ncol(x)+1+1+1]
  rho<-(exp(rho.dash) -
          1)/(exp(rho.dash) +1)
  
  #sigma2
  sigma2<-exp(theta.dash[n+ncol(x)+1+1])
  
  #gamma.dash
  df.dash<-theta.dash[n+ncol(x)+1+1+1+1]
  #
  gamma.dash<-theta.dash[n+ncol(x)+1+1+1+1+1]
  
  theta<-theta.dash
  theta[n+ncol(x)+1+1]<-sigma2 # parameters in actual space transformed space
  theta[n+ncol(x)+1+1+1]<-rho
  theta[n+ncol(x)+1+1+1+1]<-exp(df.dash)+3 
  theta[n+ncol(x)+1+1+1+1+1]<-2*exp(gamma.dash)/(exp(gamma.dash)+1)
  return(list(theta.dash=as.vector(theta.dash),theta=as.vector(theta)))
  
  
  
  
} # YJ-SEM-t


################ Functions to compare model fit using DIC1 and DIC2 





cal_DIC_SEMGau<-function(SEM.Gau.theta.sample,y,x,w){
  
  cal.log.density<-function(SEM.Gau.theta,y,x,w){
    betas<-SEM.Gau.theta[1:ncol(cbind(1,x))]
    sigma2<-SEM.Gau.theta[ncol(cbind(1,x))+1]
    rho<-SEM.Gau.theta[ncol(cbind(1,x))+2]
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(log.density))
  }
  
  # tic()
  log.densities<-apply(SEM.Gau.theta.sample, 1, cal.log.density,y,x,w)
  # toc()
  part1<--4*mean(log.densities)
  ########################
  #DIC1
  theta<-apply(SEM.Gau.theta.sample, 2,mean)
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  # psi<-theta[(p+3):length(theta)]
  n<-ncol(w)
  I<-Diagonal(n)
  
  A<-I-rho*w
  M<-t(A)%*%A
  Sigma_y<-sigma2y*solve(M,I)
  
  part2<-2*mvnfast::dmvn(y,as.vector(x%*%beta),as.matrix(Sigma_y),log = T)
  
  DIC1<-as.numeric(part1+part2)
  # print(DIC1)
  
  #DIC2
  x<-x[, -1] # in the post cal function as well as in remaining codes leading 1 will be added, so remove it
  cal.log.post<-function(SEM.Gau.theta,y,x,w){
    
    betas<-SEM.Gau.theta[1:ncol(cbind(1,x))]
    sigma2<-SEM.Gau.theta[ncol(cbind(1,x))+1]
    rho<-SEM.Gau.theta[ncol(cbind(1,x))+2]
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    # log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
    #   0.5*(1/sigma2)*t(r)%*%M%*%r
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(SEM.Gau.theta[ncol(cbind(1,x))+1])
    rho.dash<-log(1+SEM.Gau.theta[ncol(cbind(1,x))+2])-log(1-SEM.Gau.theta[ncol(cbind(1,x))+2])
    
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    ##################
    
    
    log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(rho.dash,0,omega_rho.dash,log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(log.density+log.prior))
  }
  # tic()
  log.posts<-apply(SEM.Gau.theta.sample, 1, cal.log.post,y,x,w)
  # toc()
  theta<-SEM.Gau.theta.sample[which.max(log.posts),]# find the theta that max post of theta. p(theta | y)
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  # psi<-theta[(p+3):length(theta)]
  n<-ncol(w)
  I<-Diagonal(n)
  
  A<-I-rho*w
  M<-t(A)%*%A
  Sigma_y<-sigma2y*solve(M,I)
  
  
  part2<-2*mvnfast::dmvn(y,as.vector(x%*%beta),as.matrix(Sigma_y),log = T)
  
  DIC2<-as.numeric(part1+part2)
  
  
  return(list(DIC1=DIC1,DIC2=DIC2))
  
} # SEM-Gau



cal_DIC_YJSEMGau<-function(YJ.SEM.Gau.theta.sample,y,x,w){
  
  
  
  cal_z<-function(y,gama){
    y_negative<-function(y,gama){
      return(-((-y+1)^(2-gama)-1)/(2-gama))
    }
    y_else<-function(y,gama){
      return(((y+1)^(gama)-1)/gama)
    }
    z <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
  }
  
  cal_der_tYJ_wrt_y<-function(y,gama){
    y_negative<-function(y,gama){
      return((-y+1)^(1-gama))
    }
    y_else<-function(y,gama){
      return((y+1)^(gama-1))
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  } #grad of t_jy wrt y
  
  
  cal.log.density<-function(YJ.SEM.Gau.theta,y,x,w){
    betas<-YJ.SEM.Gau.theta[1:ncol(cbind(1,x))]
    sigma2<-YJ.SEM.Gau.theta[ncol(cbind(1,x))+1]
    rho<-YJ.SEM.Gau.theta[ncol(cbind(1,x))+2]
    gama<-YJ.SEM.Gau.theta[ncol(cbind(1,x))+3]
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(log.density))
  }
  
  log.densities<-apply(YJ.SEM.Gau.theta.sample, 1, cal.log.density,y,x,w)
  
  part1<--4*mean(log.densities)
  # print(part1)
  ######################## cal part 2
  
  # DIC1
  theta<-apply(YJ.SEM.Gau.theta.sample, 2,mean)
  
  
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  gama<-theta[(p+3)]
  n<-ncol(w)
  I<-Diagonal(n)
  
  A<-I-rho*w
  M<-t(A)%*%A
  Sigma_y<-sigma2y*solve(M,I)
  
  part2<-2*mvnfast::dmvn(cal_z(y,gama),as.vector(x%*%beta),as.matrix(Sigma_y),log = T)+
    sum(log(cal_der_tYJ_wrt_y(y,gama)))*2
  
  DIC1<-as.numeric(part1+part2)
  ##DIC 2
  x<-x[, -1] # in the post cal function as well as in remaining codes leading 1 will be added, so remove it
  cal.log.post<-function(YJ.SEM.Gau.theta,y,x,w){
    
    betas<-YJ.SEM.Gau.theta[1:ncol(cbind(1,x))]
    sigma2<-YJ.SEM.Gau.theta[ncol(cbind(1,x))+1]
    rho<-YJ.SEM.Gau.theta[ncol(cbind(1,x))+2]
    gama<-YJ.SEM.Gau.theta[ncol(cbind(1,x))+3]
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    # log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
    #   0.5*(1/sigma2)*t(r)%*%M%*%r
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(YJ.SEM.Gau.theta[ncol(cbind(1,x))+1])
    rho.dash<-log(1+YJ.SEM.Gau.theta[ncol(cbind(1,x))+2])-log(1-YJ.SEM.Gau.theta[ncol(cbind(1,x))+2])
    gama.dash<-log(gama/(2-gama))
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    ##################
    
    
    log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(rho.dash,0,omega_rho.dash,log=T)+
      dnorm(gama.dash,0,omega_rho.dash,log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(log.density+log.prior))
  }
  log.posts<-apply(YJ.SEM.Gau.theta.sample, 1, cal.log.post,y,x,w)
  theta<-YJ.SEM.Gau.theta.sample[which.max(log.posts),]# find the theta that max post of theta. p(theta | y)
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  gama<-theta[(p+3)]
  n<-ncol(w)
  I<-Diagonal(n)
  A<-I-rho*w
  M<-t(A)%*%A
  
  ####
  x.beta<-(x)%*%beta
  r<-cal_z(y,gama)-x.beta
  
  part2<-2*(-0.5*n*log(2*pi)-0.5*n*log(sigma2y)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
              0.5*(1/sigma2y)*t(r)%*%M%*%r+sum(log(cal_der_tYJ_wrt_y(y,gama))))
  # print(part2)
  
  ####
  # Sigma_y<-sigma2y*solve(M,I)
  # 
  # part2<-2*mvnfast::dmvn(cal_z(y,gama),as.vector(x%*%beta),as.matrix(Sigma_y),log = T)+
  #   sum(log(cal_der_tYJ_wrt_y(y,gama)))*2
  # print(part2)
  ###
  
  DIC2<-as.numeric(part1+part2)
  
  return(list(DIC1=DIC1,DIC2=DIC2))
  
} # YJ-SEM-Gau




cal_DIC_SEMt<-function(SEM.t.lambda.dash.theta.sample,y,x,w){
  
  cal.log.density<-function(SEM.t.lambda.dash.theta,y,x,w){
    theta.para<-SEM.t.lambda.dash.theta[-(1:n)]
    betas<-theta.para[1:ncol(cbind(1,x))]
    sigma2<-theta.para[ncol(cbind(1,x))+1]
    rho<-theta.para[ncol(cbind(1,x))+2]
    df<-theta.para[ncol(cbind(1,x))+3]
    # lamda<-exp(SEM.t.lambda.dash.theta[1:n])
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    # log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
    #   0.5*(1/sigma2)*t(r)%*%M%*%r
    
    
    log.density<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(log.density))
  }
  
  log.densities<-apply(SEM.t.lambda.dash.theta.sample, 1, cal.log.density,y,x,w)
  
  part1<--4*mean(log.densities)
  # print(part1)
  ## part 2  ###
  
  # DIC1
  theta<-apply(SEM.t.lambda.dash.theta.sample, 2,mean)
  
  
  theta.para<-theta[-(1:n)]
  betas<-theta.para[1:ncol(cbind(1,x))]
  sigma2<-theta.para[ncol(cbind(1,x))+1]
  rho<-theta.para[ncol(cbind(1,x))+2]
  df<-theta.para[ncol(cbind(1,x))+3]
  # lamda<-exp(SEM.t.lambda.dash.theta[1:n])
  
  n<-nrow(w)
  A<-Diagonal(n)-rho*w
  M<-t(A)%*%A
  x.beta<-cbind(1,x)%*%betas
  r<-y-x.beta
  
  part2<-2*(lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
              lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r))
  
  DIC1<-as.numeric(part1+part2)
  
  # DIC 2, find the value that maximizes the posterior
  
  cal.log.post<-function(SEM.t.lambda.dash.theta,y,x,w){
    theta.para<-SEM.t.lambda.dash.theta[-(1:n)]
    betas<-theta.para[1:ncol(cbind(1,x))]
    sigma2<-theta.para[ncol(cbind(1,x))+1]
    rho<-theta.para[ncol(cbind(1,x))+2]
    df<-theta.para[ncol(cbind(1,x))+3]
    # lamda<-exp(SEM.t.lambda.dash.theta[1:n])
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    # log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
    #   0.5*(1/sigma2)*t(r)%*%M%*%r
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(theta.para[ncol(cbind(1,x))+1])
    rho.dash<-log(1+theta.para[ncol(cbind(1,x))+2])-log(1-theta.para[ncol(cbind(1,x))+2])
    df.dash<-log(theta.para[ncol(cbind(1,x))+3]-3)
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    ##################
    
    
    log.density<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(rho.dash,0,omega_rho.dash,log=T)+
      dnorm(df.dash,0,omega_omega.dash,log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(log.density+log.prior))
  }
  cal.log.posts<-apply(SEM.t.lambda.dash.theta.sample, 1, cal.log.post,y,x,w)
  
  
  # theta<-SEM.t.lambda.dash.theta.sample[which.max(log.densities),] 
  theta<-SEM.t.lambda.dash.theta.sample[which.max(cal.log.posts),]# find the theta that max post of theta. p(theta | y)
  
  theta.para<-theta[-(1:n)]
  betas<-theta.para[1:ncol(cbind(1,x))]
  sigma2<-theta.para[ncol(cbind(1,x))+1]
  rho<-theta.para[ncol(cbind(1,x))+2]
  df<-theta.para[ncol(cbind(1,x))+3]
  # lamda<-exp(SEM.t.lambda.dash.theta[1:n])
  
  n<-nrow(w)
  A<-Diagonal(n)-rho*w
  M<-t(A)%*%A
  x.beta<-cbind(1,x)%*%betas
  r<-y-x.beta
  
  
  
  part2<-2*(lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
              lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r))
  
  DIC2<-as.numeric(part1+part2)
  return(list(DIC1=DIC1,DIC2=DIC2))
  
} # SEM-t



cal_DIC_YJSEMt<-function(YJ.SEM.t.lambda.dash.theta.sample,y,x,w){
  
  
  
  cal_z<-function(y,gama){
    y_negative<-function(y,gama){
      return(-((-y+1)^(2-gama)-1)/(2-gama))
    }
    y_else<-function(y,gama){
      return(((y+1)^(gama)-1)/gama)
    }
    z <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
  }
  
  cal_der_tYJ_wrt_y<-function(y,gama){
    y_negative<-function(y,gama){
      return((-y+1)^(1-gama))
    }
    y_else<-function(y,gama){
      return((y+1)^(gama-1))
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  } #grad of t_jy wrt y
  
  
  cal.log.density<-function(YJ.SEM.t.lambda.dash.theta,y,x,w){
    theta.para<-YJ.SEM.t.lambda.dash.theta[-(1:n)]
    betas<-theta.para[1:ncol(cbind(1,x))]
    sigma2<-theta.para[ncol(cbind(1,x))+1]
    rho<-theta.para[ncol(cbind(1,x))+2]
    df<-theta.para[ncol(cbind(1,x))+3]
    gama<-theta.para[ncol(cbind(1,x))+4]
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    
    log.density<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(log.density))
  }
  
  log.densities<-apply(YJ.SEM.t.lambda.dash.theta.sample, 1, cal.log.density,y,x,w)
  
  part1<--4*mean(log.densities)
  # print(part1)
  ######################## cal part 2
  
  # DIC1
  theta<-apply(YJ.SEM.t.lambda.dash.theta.sample, 2,mean)
  theta.para<-theta[-(1:n)]
  betas<-theta.para[1:ncol(cbind(1,x))]
  sigma2<-theta.para[ncol(cbind(1,x))+1]
  rho<-theta.para[ncol(cbind(1,x))+2]
  df<-theta.para[ncol(cbind(1,x))+3]
  gama<-theta.para[ncol(cbind(1,x))+4]
  
  n<-nrow(w)
  A<-Diagonal(n)-rho*w
  M<-t(A)%*%A
  x.beta<-cbind(1,x)%*%betas
  r<-cal_z(y,gama)-x.beta
  
  log.density<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
    lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))
  part2<-2*log.density
  DIC1<-as.numeric(part1+part2)
  
  ##DIC 2
  
  cal.log.post<-function(YJ.SEM.t.lambda.dash.theta,y,x,w){
    theta.para<-YJ.SEM.t.lambda.dash.theta[-(1:n)]
    betas<-theta.para[1:ncol(cbind(1,x))]
    sigma2<-theta.para[ncol(cbind(1,x))+1]
    rho<-theta.para[ncol(cbind(1,x))+2]
    df<-theta.para[ncol(cbind(1,x))+3]
    gama<-theta.para[ncol(cbind(1,x))+4]
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    # log.density<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
    #   0.5*(1/sigma2)*t(r)%*%M%*%r
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(sigma2)
    rho.dash<-log(1+theta.para[ncol(cbind(1,x))+2])-log(1-theta.para[ncol(cbind(1,x))+2])
    df.dash<-log(theta.para[ncol(cbind(1,x))+3]-3)
    gama.dash<-log(gama/(2-gama))
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    ##################
    
    
    log.density<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(gama.dash,0,omega_omega.dash,log=T)+
      dnorm(rho.dash,0,omega_rho.dash,log=T)+dnorm(df.dash,0,omega_rho.dash,log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(log.density+log.prior))
  }
  cal.log.posts<-apply(YJ.SEM.t.lambda.dash.theta.sample, 1, cal.log.post,y,x,w)
  # theta<-SEM.t.lambda.dash.theta.sample[which.max(log.densities),] 
  theta<-YJ.SEM.t.lambda.dash.theta.sample[which.max(cal.log.posts),]# find the theta that max post of theta. p(theta | y)
  
  theta.para<-theta[-(1:n)]
  betas<-theta.para[1:ncol(cbind(1,x))]
  sigma2<-theta.para[ncol(cbind(1,x))+1]
  rho<-theta.para[ncol(cbind(1,x))+2]
  df<-theta.para[ncol(cbind(1,x))+3]
  gama<-theta.para[ncol(cbind(1,x))+4]
  
  n<-nrow(w)
  A<-Diagonal(n)-rho*w
  M<-t(A)%*%A
  x.beta<-cbind(1,x)%*%betas
  r<-cal_z(y,gama)-x.beta
  
  log.density<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
    lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))
  part2<-2*log.density
  
  DIC2<-as.numeric(part1+part2)
  
  return(list(DIC1=DIC1,DIC2=DIC2))
  
} # YJ-SEM-t




