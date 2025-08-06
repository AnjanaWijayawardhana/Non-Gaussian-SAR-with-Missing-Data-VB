#### Require packages

library(LaplacesDemon)
library(extraDistr)
library(igraph)
library(MASS)
library(spdep)
library(tictoc)
library(Matrix)
library(mvtnorm)
library(coda)
library(ggplot2)
library(mvnfast) # fast multivariate normal sampling
library(spatialreg)
library(moments)
library(gridExtra)


make.SATN.data<-function(x,xm,w,m,yo){
  w<-as.matrix(w)
  n<-ncol(w)
  nu<-sum(m)
  no<-n-nu
  I<-diag(n)
  
  
  
  x_<-cbind(1,x)
  xm<-cbind(1,xm)
  q<-ncol(xm)+1
  xo=x_[1:no,]
  xu=x_[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  p=ncol(x_)
  
  ############## priors
  sigma2_beta=1000*diag(p)
  mu_beta=rep(0,p)
  sigma2_psi=1000*diag(q)
  mu_psi=rep(0,q)
  
  stan.data.list<-list(n=n,p=p,q=q,no=no,nu=nu,xm=xm,x=x_,xu=xu,xo=xo,yo=yo,m=m,
                       w=w,I=I,wpluwt=wpluswt,
                       wtw=wtw,mu_beta=mu_beta, sigma2_beta=sigma2_beta,
                       mu_psi=mu_psi,sigma2_psi=sigma2_psi)
  
  return(stan.data.list=stan.data.list)
  
  
} # Prepares and returns a list of input data required for Bayesian modeling in STAN. 
# Includes covariates, spatial weight matrix, priors, and other derived matrices.





######################################################  YJ-SEM-Gau codes ######################################################



# In MCMC stpe proposals are from zu |zo, but before calculating accept.prob, 
#Generate zu |zo-> yu |yo. Then accept.prob is calculated only including missing val model. 
YJ_SEM_Gau_MNAR_NoB<-function(x,xm,yo,m,N,w, # This is the correct, working code
                              start.theta,start.yu,p,N1){ 
  # p=1
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rc=ncol(x) #number of regression covariates. i.e. betas
  rcm=ncol(xm)+1
  totpara<-rc+3+rcm
  adapt_epsilon=10^(-6)
  v=0.99
  
  n<-length(y)
  no<-length(yo)
  nu<-n-no
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  
  # theta_u<-c(c(1,2),0.5,0.1,c(3,3,3))
  cal_z<-function(y,gama){
    y_negative<-function(y,gama){
      return(-((-y+1)^(2-gama)-1)/(2-gama))
    }
    y_else<-function(y,gama){
      return(((y+1)^(gama)-1)/gama)
    }
    z <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(z)
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
  
  
  cal_der__log_der_tYJ_wrt_y__wrt_gamma<-function(y,gama){
    y_negative<-function(y,gama){
      return(-log(-y+1))
    }
    y_else<-function(y,gama){
      return(log(y+1))
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  }#der wrt gamma
  #of log of der wrt y.
  
  
  
  grad<-function(theta_u,wpluswt,wtw,w,x,xm,m,I,yo,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p mean dimension of X.
    
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    omega.dash<-theta_u[(nu+p+1)] #transform of sigma2
    rho.dash<-theta_u[(nu+p+2)] #transform of rho
    gama.dash<-theta_u[(nu+p+3)] #transform of gama
    psi<-theta_u[(nu+p+4):length(theta_u)]
    
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_omega.dash=10000
    omega_rho.dash=10000
    omega_gama.dash=10000
    omega_psi=10000
    
    A=I-rho*w
    M=t(A)%*%A
    y<-c(yo,yu)
    z<-cal_z(y,gama)
    
    #betas
    div_h_beta<-exp(-omega.dash)*t(z-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #omega.dash
    div_h_omega.dash<--0.5*n+(0.5*exp(-omega.dash))*t(z-x%*%beta)%*%M%*%(z-x%*%beta)-omega.dash/omega_omega.dash
    
    #rho.dash
    drho_by_drho.dash<-2*exp(rho.dash)/(1+exp(rho.dash))^2
    
    dM_by_drho.dash<-(-wpluswt+2*rho*wtw)*drho_by_drho.dash
    
    #1
    # dlogdetM_by_drho.dash<-sum(diag(solve(M)%*%dM_by_drho.dash))
    #2
    # dlogdetM_by_drho.dash<-sum(diag(solve(M,dM_by_drho.dash))) #For real
    #3
    dlogdetM_by_drho.dash<-sum(diag(solve(Cholesky(M),dM_by_drho.dash))) # For simulations
    #
    
    
    div_h_rho.dash<-0.5*dlogdetM_by_drho.dash-0.5*exp(-omega.dash)*t(z-x%*%beta)%*%dM_by_drho.dash%*%(z-x%*%beta)-rho.dash/omega_rho.dash
    
    #gama.dash
    d_gama_by_gama.dash<-2*exp(gama.dash)/(1+exp(gama.dash))^2
    
    der_tYJ_wrt_gama<-cal_der_tYJ_wrt_gama(y,gama)
    part1<--exp(-omega.dash)*t(z-x%*%beta)%*%M%*%der_tYJ_wrt_gama*d_gama_by_gama.dash
    
    # second_der_wrt_y_gama<-cal_second_der_wrt_y_gama(y,gama)
    # der_tYJ_wrt_y<-cal_der_tYJ_wrt_y(y,gama)
    # part2<-sum((1/der_tYJ_wrt_y)*second_der_wrt_y_gama*d_gama_by_gama.dash)
    ##
    part2<-sum(cal_der__log_der_tYJ_wrt_y__wrt_gamma(y,gama))*d_gama_by_gama.dash
    ##
    
    div_h_gama.dash<-part1+part2-gama.dash/omega_gama.dash
    
    
    #psi
    zmiss<-cbind(xm,y)
    div_h_psi<-(t(m-(exp(zmiss%*%psi)/((1+exp(zmiss%*%psi)))))%*%zmiss)-t(psi)/omega_psi
    
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_omega.dash),as.vector(div_h_rho.dash),as.vector(div_h_gama.dash),
                as.vector(div_h_psi))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  sim_yu.g.yom.theta<-function(x,xu,xo,xm,w,wpluswt,wtw,yo,yu.start,m,beta, # New yu | yo, using bp direct (not the function in "VB_YJ.SEM.Gaus_MNAR_BP_Direct.yu.gyo")
                               sigma2,rho,gama,psi,N){
    log.target<-function(m,yo,yu,x,xm,psi,gama){
      
      y=c(yo,yu)
      lin.comb=cbind(xm,y)%*%psi
      pxy<-1/(1+exp(-lin.comb))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    con.dis.ugo<-function(beta,rho,sigma2,zo,xo,xu,w,wpluswt,wtw,I){
      M=I-rho*wpluswt+rho*rho*wtw
      no<-nrow(xo)
      nu<-nrow(xu)
      
      M_oo=M[1:no,1:no]
      M_uu=M[(no+1):n,(no+1):n]
      M_ou=M[1:no,(no+1):n]
      M_uo=t(M_ou)
      In<-Diagonal(nu)
      cholMuu<-Cholesky(M_uu)
      mu_ugo_=xu%*%beta-solve(cholMuu,M_uo)%*%(zo-xo%*%beta)
      sigma_ugo_=forceSymmetric(sigma2*solve(cholMuu,In))
      # toc()
      
      return(list(as.vector(mu_ugo_),
                  as.matrix(sigma_ugo_)))
    }
    
    nu=length(yu.start)
    samples<-matrix(rep(0,nu*N),ncol = nu)
    
    accepts<-0
    zo<-cal_z(yo,gama)
    co.dis<-con.dis.ugo(beta,rho,sigma2,zo,xo,xu,w,wpluswt,wtw,I)
    con.mean<-co.dis[[1]]
    con.var<-co.dis[[2]]
    
    # tic()
    zu.star.mat=rmvn(N,as.vector(con.mean),as.matrix(con.var))
    yu.star.mat<-apply(zu.star.mat, 2, t_inverse,gama)
    
    # toc()
    yu.star.mat[1,]<-yu.start # replace first relement in yu.mat with yu.start
    yu.chain<-matrix(rep(0,ncol(yu.star.mat)*nrow(yu.star.mat)),ncol=ncol(yu.star.mat))
    yu.chain[1,]<- yu.star.mat[1,]
    current_state<- yu.star.mat[1,]
    
    
    for(i in 2:N){#i=2
      
      yu.star=yu.star.mat[i,] # start adding new yu from the second element
      ration=min(0,log.target(m,yo,yu=yu.star,x,xm,psi,gama)-log.target(m,yo,yu=current_state,x,xm,psi,gama))
      
      
      if ((!is.na(ration))){
        if((log(runif(1)) <= ration))
        {
          current_state <- yu.star
          accepts<-accepts+1
        }
        
      }
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    return(list("yu.chain"=yu.chain,"accepts"=accepts))
    
  }
  
  
  
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
  accepts<-c()
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
    
    
    # sampling step
    beta<-theta[(1:rc)]
    omega.dash.rho.dash<-theta[(rc+1):(rc+2)]
    gama.dash<-theta[rc+3]
    
    rho<-(exp(omega.dash.rho.dash[2])-1)/(exp(omega.dash.rho.dash[2])+1)
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    omega<-exp(omega.dash.rho.dash[1])
    # N1=10
    psi<-theta[-(1:(rc+3))]
    mysamples<-sim_yu.g.yom.theta(x,xu,xo,xm,w,wpluswt,wtw,yo,yu.start=start.yu,m,
                                  beta,sigma2=omega,rho,gama,psi,N=N1)
    # mysamples<-sim_yu_g_yomtheta(x,xu,xo,xm,w,wpluswt,wtw,yo,yu=start.yu,m,beta=c(2.0,3.0),sigma2y=1,rho=0.8,gama=0.5,psi=c(-1,rep(0.5,1),-0.1),N=N1)
    accepts[i]<-mysamples$accepts
    yunote<-(mysamples$yu.chain)[N1,]
    start.yu<-yunote 
    theta_u<-c(as.vector(yunote),as.vector(theta)) # combine thetas and yu
    
    gradiant<-grad(theta_u=theta_u,wpluswt,wtw,w,x,xm,m,I,yo,n)
    
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
    all_paras[i,]<-as.vector(mu)
    
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,
              all_paras=all_paras,accepts=accepts))
} # VB approximation, B is treated as a vector and then transfer into a matrix 


# In MCMC stpe proposals are from zu |zo, but before calculating accept.prob, 
#Generate zu |zo-> yu |yo. Then accept.prob is calculated only including missing val model. 



YJ_SEM_Gau_MNAR_Buall_new2<-function(x,xm,yo,m,N,w, # This is the correct, working code
                                     start.theta,start.yu,p,bsize,N1){ 
  # p=1
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rc=ncol(x) #number of regression covariates. i.e. betas
  rcm=ncol(xm)+1
  totpara<-rc+3+rcm
  adapt_epsilon=10^(-6)
  v=0.99
  
  # bsize=1000
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  
  n<-nrow(x)
  no<-length(yo)
  nu<-n-no
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  
  # theta_u<-c(c(1,2),0.5,0.1,c(3,3,3))
  cal_z<-function(y,gama){
    y_negative<-function(y,gama){
      return(-((-y+1)^(2-gama)-1)/(2-gama))
    }
    y_else<-function(y,gama){
      return(((y+1)^(gama)-1)/gama)
    }
    z <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(z)
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
  
  
  cal_der__log_der_tYJ_wrt_y__wrt_gamma<-function(y,gama){
    y_negative<-function(y,gama){
      return(-log(-y+1))
    }
    y_else<-function(y,gama){
      return(log(y+1))
    }
    der_tYJ <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(der_tYJ)
  }
  
  
  
  grad<-function(theta_u,wpluswt,wtw,w,x,xm,m,I,yo,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p mean dimension of X.
    
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    omega.dash<-theta_u[(nu+p+1)] #transform of sigma2
    rho.dash<-theta_u[(nu+p+2)] #transform of rho
    gama.dash<-theta_u[(nu+p+3)] #transform of gama
    psi<-theta_u[(nu+p+4):length(theta_u)]
    
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_omega.dash=10000
    omega_rho.dash=10000
    omega_gama.dash=10000
    omega_psi=10000
    
    A=I-rho*w
    M=t(A)%*%A
    y<-c(yo,yu)
    z<-cal_z(y,gama)
    
    #betas
    div_h_beta<-exp(-omega.dash)*t(z-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #omega.dash
    div_h_omega.dash<--0.5*n+(0.5*exp(-omega.dash))*t(z-x%*%beta)%*%M%*%(z-x%*%beta)-omega.dash/omega_omega.dash
    
    #rho.dash
    drho_by_drho.dash<-2*exp(rho.dash)/(1+exp(rho.dash))^2
    
    dM_by_drho.dash<-(-wpluswt+2*rho*wtw)*drho_by_drho.dash
    dlogdetM_by_drho.dash<-sum(diag(solve(Cholesky(M),dM_by_drho.dash))) # For simulations with rook W
    # dlogdetM_by_drho.dash<-sum(diag(solve(M)%*%dM_by_drho.dash)) # For real data
    
    div_h_rho.dash<-0.5*dlogdetM_by_drho.dash-0.5*exp(-omega.dash)*t(z-x%*%beta)%*%dM_by_drho.dash%*%(z-x%*%beta)-rho.dash/omega_rho.dash
    
    #gama.dash
    d_gama_by_gama.dash<-2*exp(gama.dash)/(1+exp(gama.dash))^2
    
    der_tYJ_wrt_gama<-cal_der_tYJ_wrt_gama(y,gama)
    part1<--exp(-omega.dash)*t(z-x%*%beta)%*%M%*%der_tYJ_wrt_gama*d_gama_by_gama.dash
    
    # second_der_wrt_y_gama<-cal_second_der_wrt_y_gama(y,gama)
    # der_tYJ_wrt_y<-cal_der_tYJ_wrt_y(y,gama)
    # part2<-sum((1/der_tYJ_wrt_y)*second_der_wrt_y_gama*d_gama_by_gama.dash)
    
    ##
    part2<-sum(cal_der__log_der_tYJ_wrt_y__wrt_gamma(y,gama))*d_gama_by_gama.dash
    ##
    
    div_h_gama.dash<-part1+part2-gama.dash/omega_gama.dash
    
    
    #psi
    zmiss<-cbind(xm,y)
    # div_h_psi<-(t(m-(exp(zmiss%*%psi)/((1+exp(zmiss%*%psi)))))%*%zmiss)-t(psi)/omega_psi
    div_h_psi<-t(zmiss)%*%(m-(1/(1+exp(-zmiss%*%psi))))-matrix(t(psi)/omega_psi,ncol=1) # The above one also correct.
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_omega.dash),as.vector(div_h_rho.dash),as.vector(div_h_gama.dash),
                as.vector(div_h_psi))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  
  make.w_x_<-function(w,x,xm,m,k,reminder){
    
    if(reminder==0){
      
      w.list<-list()
      x.list<-list()
      xm.list<-list()
      m.list<-list()
      
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
        # w[1:15,1:15]
        xui<-matrix(x[ui,],nrow=bsize,)
        xui_<-as.matrix(x[ui_,],nrow=nu-bsize)
        xo<-as.matrix(x[o,])
        x.list[[i]]<-rbind(xui,rbind(xui_,xo))
        
        xmui<-matrix(xm[ui,],nrow=bsize,)
        xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
        xmo<-as.matrix(xm[o,])
        xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
        
        m.list[[i]]<-rev(m)
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        w.list<-list()
        x.list<-list()
        xm.list<-list()
        m.list<-list()
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=bsize,)
            xui_<-matrix(x[ui_,],nrow=nu-bsize)
            xo<-x[o,]
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=bsize,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=reminder,)
            xui_<-as.matrix(x[ui_,],nrow=nu-reminder)
            xo<-as.matrix(x[o,])
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=reminder,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-reminder)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(w.list=w.list,x.list=x.list,xm.list=xm.list,m.list=m.list))
  }
  
  sim_yu.g.yom.theta.uallb<-function(zo,yo,listw_x_m,    # sim zu g zo not blocking, using bp 
                                     start.yu,beta,
                                     sigma2,rho,gama,psi,
                                     N1,k){
    nu=length(start.yu)
    
    
    cal.con.dis<-function(x,zo,zui_,w,I,beta,rho,sigma2,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      M<-t(I-rho*w)%*%(I-rho*w)
      if(reminder==0){
        
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(zui_,zo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the incomplete block
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo) # ui given ui_not and yo
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2*solve(cholMuj.uj,In))
          
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yo,yu)
      lin.comb=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-lin.comb))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    yu.chain[1,]<-start.yu
    current_state<-start.yu
    
    # tic()
    accepts<-0
    for(i in 2:N1){#i=3 N=5
      
      
      # j=1
      current<-0
      for (j in 1:(k+ifelse(reminder>0,1,0))) {
        # accepts<-0
        # yubj<-current_state[(current+1):(current+bsize)]
        if(j<=k){ # for full blocks j=2
          zubj_<-cal_z(current_state[-((current+1):(current+bsize))],gama)
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2=sigma2,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=t_inverse(rmvn(1,as.vector(con.mean),as.matrix((con.var))),gama)
          proposal_state<-current_state
          
          proposal_state[(current+1):(current+bsize)]<-yubj.star
          
          # current<-5
          ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                       log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
          
          
          
          if ((!is.na(ration))){
            if((log(runif(1)) <= ration))
            {
              current_state <- proposal_state
              accepts<-accepts+1
            }
            
          }
          
          current<-current+bsize
        }else{ # for the half block j=2 current=5
          zubj_<-cal_z(current_state[-((current+1):(current+bsize))],gama)
          # current=10
          # print("in half block")
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2=sigma2,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          # print("out half block")
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=t_inverse(rmvn(1,as.vector(con.mean),as.matrix((con.var))),gama)
          
          
          proposal_state<-current_state
          
          proposal_state[(current+1):(current+reminder)]<-yubj.star
          
          # current<-5
          ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                       log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
          
          
          
          if ((!is.na(ration))){
            if((log(runif(1)) <= ration))
            {
              current_state <- proposal_state
              accepts<-accepts+1
            }
            
          }
          
        }
      }
      
      
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(yugzom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  
  
  
  
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
  accepts<-c()
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
    
    
    # sampling step
    beta<-theta[(1:rc)]
    omega.dash.rho.dash<-theta[(rc+1):(rc+2)]
    gama.dash<-theta[rc+3]
    
    rho<-(exp(omega.dash.rho.dash[2])-1)/(exp(omega.dash.rho.dash[2])+1)
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    omega<-exp(omega.dash.rho.dash[1])
    
    
    
    psi<-theta[-(1:(rc+3))]
    # start.zu<-cal_z(start.yu,gama)
    zo<-cal_z(yo,gama)
    
    
    fullblockmysamples<-sim_yu.g.yom.theta.uallb(zo,yo,listw_x_m,    # After adding un-balanced # of blocks 
                                                 start.yu,beta,
                                                 sigma2=omega,rho,gama,psi,
                                                 N1,k)
    
    # mysamples<-sim_zu.g.zom.theta(x,xu,xo,xm,w,wpluswt,wtw,zo,
    #                               zu=start.zu,m,beta,sigma2=omega,rho,
    #                               psi,N=N1)
    accepts[i]<-fullblockmysamples$accepts
    
    yunote<-(fullblockmysamples$yugzom.chain)[N1,]
    # yunote<-t_inverse(zunote,gama)
    start.yu<-yunote
    theta_u<-c(as.vector(yunote),as.vector(theta)) # combine thetas and yu
    
    
    # yunote<-start.yu
    # theta_u<-c(as.vector(start.yu),as.vector(theta))
    # 
    # 
    # theta<-c(yu,theta)
    # Calculate log gradient of theta
    gradiant<-grad(theta_u=theta_u,wpluswt,wtw,w,x,xm,m,I,yo,n)
    
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
    
    # Lbound[i]<-0
    
    all_paras[i,]<-as.vector(mu)
    # print(i)
    # Lbound[i]<-0
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,
              all_paras=all_paras,accepts=accepts))
} # VB approximation, B is treated as a vector and then transfer into a matrix 




######################################################  Functions to simulate data and make missing values ################


MNAR_splitter_all_xm<-function(x,xm,y,w,psi){
  n=nrow(y)
  x<-matrix(x,nrow = n)
  xm<-matrix(xm,nrow = n)
  lead_1=rep(1,n)
  z=cbind(lead_1,xm,y)%*%psi 
  pr= 1/(1+exp(-z)) 
  #m=rbinom(nrow(y),1,pr) 
  
  u=runif(nrow(y))
  
  m=as.numeric(pr>u)
  
  y_us=rep(0,n)
  
  for(i in 1:n){
    if(m[i]==1){
      y_us[i]=NA
    }else{
      y_us[i]=y[i]
    }
  }
  
  ####
  
  nu<-sum(is.na(y_us))
  no<-n-nu
  
  
  missingindexes<-which(is.na(y_us))
  xu<-x[missingindexes,]
  xo<-x[setdiff(1:n,missingindexes),]
  yo<-y[setdiff(1:n,missingindexes)]
  yu<-y[missingindexes]
  xmu<-xm[missingindexes,]
  xmo<-xm[setdiff(1:n,missingindexes),]
  u<-missingindexes
  o<-setdiff(1:n,u)
  
  w<-w[c(o,u),c(o,u)]
  mo<-m[o]
  mu_<-m[u] #To overcome the confusion with mean of prior, which is also mu
  m<-c(mo,mu_)
  
  x=rbind(matrix(xo,nrow=no),matrix(xu,nrow=nu))
  xm=rbind(matrix(xmo,nrow=no),matrix(xmu,nrow=nu))
  return(list("m"=m,"X"=x,"xm"=xm,"yo"=yo,"yu"=yu,"W"=w))
  
}  #Generate missing values undeer MNAR, and 
                                                  #create partioned matrices, output, w=(woo,wou,
                                                  #wuo,wuu), x=(xo,xu),yo,yu, m=(mo,mu)



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
  
} # Generate SSEM, efficient generator



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
} # this calculates y from z~N using the inverse of YJ transformation

row_standardize_to_sum_1 <- function(w) {
  # Get the row sums
  row_sums <- rowSums(w)
  
  # Avoid division by zero (if there are any rows with sum 0, handle them separately)
  row_sums[row_sums == 0] <- 1  # Set zero sums to 1 to avoid division by zero
  
  # Divide each element in a row by the sum of that row
  standardized_w <- sweep(w, 1, row_sums, "/")
  
  return(standardized_w)
}



##### Functions for draw from the posterior distributions 


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
  
  
  
  
} # Generate draws of model parameters 




generate_posterior_yu_YJSEMGau<-function(theta,yo,start.yu,
                                x,xm,w,m,bsize,N1){
  
  
  
  
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  p=ncol(x) #number of regression covariates. i.e. betas
  
  
  
  ### parameters
  
  
  beta<-theta[(1):(p)]
  omega.dash<-theta[(p+1)] #transform of sigma2
  rho.dash<-theta[(p+2)] #transform of rho
  gama.dash<-theta[(p+3)] #transform of gama
  psi<-theta[(p+4):length(theta)]
  
  # Dont need to do tranfomation, as theta is in the original scale
  gama<-gama.dash
  rho<-rho.dash
  omega<-(omega.dash)
  
  
  
  # bsize=1000
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  
  n<-nrow(x)
  no<-length(yo)
  nu<-n-no
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  
  
  sim_yu.g.yom.theta.uallb<-function(zo,yo,listw_x_m,    # sim zu g zo not blocking, using bp 
                                     start.yu,beta,
                                     sigma2,rho,gama,psi,
                                     N1,k){
    nu=length(start.yu)
    
    
    cal.con.dis<-function(x,zo,zui_,w,I,beta,rho,sigma2,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      M<-t(I-rho*w)%*%(I-rho*w)
      if(reminder==0){
        
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(zui_,zo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the incomplete block
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo) # ui given ui_not and yo
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2*solve(cholMuj.uj,In))
          
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yo,yu)
      lin.comb=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-lin.comb))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    yu.chain[1,]<-start.yu
    current_state<-start.yu
    
    # tic()
    accepts<-0
    for(i in 2:N1){#i=3 N=5
      
      
      # j=1
      current<-0
      for (j in 1:(k+ifelse(reminder>0,1,0))) {
        # accepts<-0
        # yubj<-current_state[(current+1):(current+bsize)]
        if(j<=k){ # for full blocks j=2
          zubj_<-cal_z(current_state[-((current+1):(current+bsize))],gama)
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2=sigma2,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=t_inverse(rmvn(1,as.vector(con.mean),as.matrix((con.var))),gama)
          proposal_state<-current_state
          
          proposal_state[(current+1):(current+bsize)]<-yubj.star
          
          # current<-5
          ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                       log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
          
          
          
          if ((!is.na(ration))){
            if((log(runif(1)) <= ration))
            {
              current_state <- proposal_state
              accepts<-accepts+1
            }
            
          }
          
          current<-current+bsize
        }else{ # for the half block j=2 current=5
          zubj_<-cal_z(current_state[-((current+1):(current+bsize))],gama)
          # current=10
          # print("in half block")
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2=sigma2,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          # print("out half block")
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=t_inverse(rmvn(1,as.vector(con.mean),as.matrix((con.var))),gama)
          
          
          proposal_state<-current_state
          
          proposal_state[(current+1):(current+reminder)]<-yubj.star
          
          # current<-5
          ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                       log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
          
          
          
          if ((!is.na(ration))){
            if((log(runif(1)) <= ration))
            {
              current_state <- proposal_state
              accepts<-accepts+1
            }
            
          }
          
        }
      }
      
      
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(yugzom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  
  
  make.w_x_<-function(w,x,xm,m,k,reminder){
    
    if(reminder==0){
      
      w.list<-list()
      x.list<-list()
      xm.list<-list()
      m.list<-list()
      
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
        # w[1:15,1:15]
        xui<-matrix(x[ui,],nrow=bsize,)
        xui_<-as.matrix(x[ui_,],nrow=nu-bsize)
        xo<-as.matrix(x[o,])
        x.list[[i]]<-rbind(xui,rbind(xui_,xo))
        
        xmui<-matrix(xm[ui,],nrow=bsize,)
        xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
        xmo<-as.matrix(xm[o,])
        xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
        
        m.list[[i]]<-rev(m)
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        w.list<-list()
        x.list<-list()
        xm.list<-list()
        m.list<-list()
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=bsize,)
            xui_<-matrix(x[ui_,],nrow=nu-bsize)
            xo<-x[o,]
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=bsize,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=reminder,)
            xui_<-as.matrix(x[ui_,],nrow=nu-reminder)
            xo<-as.matrix(x[o,])
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=reminder,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-reminder)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(w.list=w.list,x.list=x.list,xm.list=xm.list,m.list=m.list))
  }
  
  cal_z<-function(y,gama){
    y_negative<-function(y,gama){
      return(-((-y+1)^(2-gama)-1)/(2-gama))
    }
    y_else<-function(y,gama){
      return(((y+1)^(gama)-1)/gama)
    }
    z <- ifelse(y<0, y_negative(y,gama), y_else(y,gama))
    return(z)
  }
  
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  
  zo<-cal_z(yo,gama)
  
  
  fullblockmysamples<-sim_yu.g.yom.theta.uallb(zo,yo,listw_x_m,    # After adding un-balanced # of blocks 
                                               start.yu,beta,
                                               sigma2=omega,rho,gama,psi,
                                               N1,k)
  
  return(apply((fullblockmysamples$yugzom.chain)[((N1/2):N1),],2,mean))
  
} # Generate draws of missing values





