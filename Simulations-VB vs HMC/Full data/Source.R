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
library(rstan)
library(MASS)
library(ggplot2)
library(reshape2)


########### Functions for simulation and VB approximation

make.SATN.data<-function(x,w,y){
  w<-as.matrix(w)
  n<-ncol(w)
  I<-diag(n)
  
  
  
  x_<-cbind(1,x)
  
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  p=ncol(x_)
  
  ############## priors
  sigma2_beta=1000*diag(p)
  mu_beta=rep(0,p)
  
  
  
  stan.data.list<-list(n=n,p=p,x=x_,y=y,
                       w=w,I=I,wpluwt=wpluswt,
                       wtw=wtw,mu_beta=mu_beta, sigma2_beta=sigma2_beta)
  
  return(stan.data.list=stan.data.list)
  
  
} # Prepares and returns a list of input data required for Bayesian modeling in STAN. 
                                     # Includes covariates, spatial weight matrix, priors, and other derived matrices.





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



YJ_SEM_Gau<-function(x,y,N,w,
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
} # Variational Bayes (VB) algorithm for estimating YJ-SEM-Gau.



################## Function for draw from the posterior distribution


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
  
  
  
  
} # Draws posterior samples of model parameters from the approximated distribution 


