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


########### Functions to generate data and make missing values





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
  
}  #spliting and re arranging, output, w=(woo,wou,
#wuo,wuu), x=(xo,xu),yo,yu, m=(mo,mu)



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


################## SEM-Gau, SEM-t, YJ-SEM-Gau, and YJ-SEM-t algorithms (HVB-AllB)

SEM_Gau_MNAR_Buall<-function(x,xm,m,w,yo,startyu,p, #SEM-Gau
                             start_theta,N,N1,bsize){ 
  # p=10
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=ncol(xm)+1     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  
  n<-length(m)
  no<-length(yo)
  nu<-n-no
  
  
  # bsize=10
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  
  grad_2<-function(theta_u,wpluswt,wtw,w,x,I,yo,xo,xu,xm,m,mu_,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    psi<-theta_u[(nu+p+3):length(theta_u)]
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    omega_psi=10000
    
    
    y=c(yo,yu)
    # A=I-rho*w
    # M=t(A)%*%A
    M=I-rho*wpluswt+rho*rho*wtw
    #betas
    div_h_beta<-exp(-gama)*t(y-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #gamma
    div_h_gamma<--0.5*n+(0.5*exp(-gama))*t(y-x%*%beta)%*%M%*%(y-x%*%beta)-gama/omega_gamma
    
    #Lamda
    drho_by_dlamda<-2*exp(lamda)/(1+exp(lamda))^2
    
    dM_by_dlamda<-(-wpluswt+2*rho*wtw)*drho_by_dlamda
    
    dlogdetM_by_dlamda<-sum(diag(solve((M),dM_by_dlamda))) #For real
    #
    # M.in.dMbdlamda<-solve(Cholesky(M),dM_by_dlamda) #2, For simulations
    # dlogdetM_by_dlamda<-sum(diag(M.in.dMbdlamda))
    
    #
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    #Psi
    z<-cbind(xm,y)
    div_h_psi<-(t(m-(exp(z%*%psi)/((1+exp(z%*%psi)))))%*%z)-t(psi)/omega_psi
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda),as.vector(div_h_psi))
    
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
  
  sim_yugyotheta_<-function(yo,listw_x_m,    # After adding un-balanced # of blocks 
                            yu.start,beta,
                            sigma2y,rho,psi,
                            N1,k){
    nu=length(yu.start)
    
    
    cal.con.dis<-function(x,yo,yui_,w,I,beta,rho,sigma2y,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      M<-t(I-rho*w)%*%(I-rho*w)
      # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
      if(reminder==0){
        # Sigma_ui.ui_o<-Sigma[1:bsize,(bsize+1):n] #between cov(yui and yui_,yo)
        # Sigma_ui_o.ui_o<-Sigma[(bsize+1):n,(bsize+1):n] #cov(yui_,yo and yui_,yo)
        # Sigma_uiui<-Sigma[1:bsize,1:bsize]
        # 
        # 
        # mu_uigui_o<-x[1:bsize,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
        #                                                      c(yui_,yo)-(x[(bsize+1):n,])%*%beta)
        # sigma_uigui_o<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
        
        
        # sigma_uj.v<-Sigma[1:bsize,(bsize+1):n]
        # sigma_v.v<-Sigma[(bsize+1):n,(bsize+1):n] 
        # sigma_uj.uj<-Sigma[1:bsize,1:bsize] 
        # mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        # v<-c(yui_,yo)
        # 
        # mu_ujguv<-x[1:bsize,]%*%beta+sigma_uj.v%*%solve(sigma_v.v,(v-mu_v))
        # sigma_ujgv<-sigma_uj.uj-sigma_uj.v%*%solve(sigma_v.v,t(sigma_uj.v))
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
          # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
          # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
          # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
          # 
          # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
          #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
          # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
          
          ################################
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yu,yo)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    yu.chain[1,]<-yu.start
    current_state<-yu.start
    
    # tic()
    accepts<-0
    for(i in 2:N1){#i=3 N=5
      
      
      # j=1
      current<-0
      for (j in 1:(k+ifelse(reminder>0,1,0))) {
        # accepts<-0
        # yubj<-current_state[(current+1):(current+bsize)]
        if(j<=k){ # for full blocks j=2
          yubj_<-current_state[-((current+1):(current+bsize))]
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
          yubj_<-current_state[-((current+1):(current+reminder))]
          # current=10
          # print("in half block")
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          # print("out half block")
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  #to calculate lower bound
  #initial values for lamdha
  #p=2
  #mu<-rep(0.1,totpara)
  
  mu<-c(start_theta)
  startyu<-yu
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
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) #To store lamda in each iteration
  tlowerbound<-c()
  accepts<-c()
  
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  
  #i=1
  for(i in 1:N){
    
    # thetas
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    theta<-mu+B%*%z+d*as.vector(epsilon)
    D<-diag(d)
    
    # generate yus from p(yu | yo, theta)
    beta<-theta[(1:rc)]
    gamma.lamda<-theta[(rc+1):(rc+2)]
    rho<-(exp(gamma.lamda[2])-1)/(exp(gamma.lamda[2])+1)
    sigma2y<-exp(gamma.lamda[1])
    psi<-theta[-(1:(rc+2))]
    # N1=100
    
    mysamples<-sim_yugyotheta_(yo=yo,listw_x_m=listw_x_m,
                               yu.start=startyu,beta=beta,
                               sigma2y=sigma2y,rho=rho,
                               psi=psi,N1=N1,k=k)
    # plot(ts((mysamples$samples)[,1]))
    accepts[i]<-mysamples$accepts
    # dim(mysamples$samples)
    # yunote<-apply((mysamples$samples)[-(1:N1/2),], 2, mean)
    yunote<-(mysamples$yugyom.chain)[N1,]
    startyu<-yunote # set starting yu of the MCMC samples of the next iteration
    
    theta_u<-c(as.vector(yunote),as.vector(theta)) # combine thetas and yu
    
    # tic()
    gradg_theta<-grad_2(theta_u=theta_u,w=w,x,I=I,yo=yo,xo=xo,xu=xu,xm=xm,
                        wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n,mu_=mu_,m=m) # gradient of log h(theta)
    
    # toc()
    
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    #Construct UB estimates
    grad_mu<-gradg_theta-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradg_theta))
    grad_B<-t(t(z)%x%I_theta)%*%(gradg_theta-dlogq_by_dtheta)
    gradd<-epsilon*(gradg_theta-dlogq_by_dtheta)
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
    
    # calculate lowerbound
    
    
    # tlowerbound[i]=system.time(Lbound[i]<-cal_Lbound(theta_u=theta,nu=nu,mu=mu,d=d,B=B,m=m,
    # x=x,yo=yo,wpluswt=wpluswt,wtw=wtw))[[3]]
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    
    
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
    # print(i)
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,
              all_paras=all_paras,accepts=accepts,
              dis_yug_yo=mysamples[[3]]))
} 



YJ_SEM_Gau_MNAR_Buall_new2<-function(x,xm,yo,m,N,w, # YJ-SEM-Gau
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
    # dlogdetM_by_drho.dash<-sum(diag(solve(Cholesky(M),dM_by_drho.dash))) # For simulations with rook W
    dlogdetM_by_drho.dash<-sum(diag(solve(M)%*%dM_by_drho.dash)) # For real data
    
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
} 

SEM_t_MNAR_Buall_noninfo<-function(x,xm,m,w,yo,yu,startyu,p, 
                                   start_theta,N,bsize,N1,shape=4,rate=1/2){
  
  
  rho.a<--1
  rho.b<-1
  # p=10
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  numcov<-ncol(x)
  xm<-cbind(1,xm)
  n<-length(m)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=ncol(xm)+1     #number of logistic covariates. i.e. betas
  totpara<-n+rc+lc+3
  adapt_epsilon=10^(-6)
  v=0.95
  
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  
  no<-length(yo)
  nu<-n-no
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  f<-function(df.dash){(digamma(0.5*exp(df.dash)+1.5))}
  f2<-function(df.dash){(lgamma(0.5*exp(df.dash)+1.5))}
  
  
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
  
  make.lamda.partion<-function(lamda,k,reminder){
    
    if(reminder==0){
      
      lamda.list<-list()
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        lamda.list[[i]]<-lamda[c(ui,ui_,o)]
        
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        lamda.list<-list()
        
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(lamda.list=lamda.list))
  }
  
  grad_2_new<-function(yu_lamda_theta,wpluswt,wtw,w,x,I,yo,xo,xu,xm, # SEM-Gau
                       m,nu,no,n){
    
    
    
    
    #################
    numcov<-ncol(x)
    yu<-yu_lamda_theta[1:nu]
    lamda.dash<-yu_lamda_theta[(nu+1):(n+nu)]
    beta<-yu_lamda_theta[(nu+n+1):(n+nu+numcov)]
    omega.dash<-yu_lamda_theta[(nu+n+numcov+1)] #transform of sigma2 
    rho.dash<-yu_lamda_theta[(nu+n+numcov+2)] #transform of rho 
    df.dash<-yu_lamda_theta[(nu+n+numcov+3)]
    psi<-yu_lamda_theta[(nu+n+numcov+4):length(yu_lamda_theta)]
    
    lamda<-exp(lamda.dash)
    omega<-exp(omega.dash)
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    df<-exp(df.dash)+3
    
    #################
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_omega.dash=10000
    omega_rho.dash=10000
    omega_psi=10000
    omega_df=10000
    
    y=c(yo,yu)
    # A=I-rho*w
    # M=t(A)%*%A
    
    # Precalculations
    r<-y-x%*%beta
    A<-I-rho*w
    Ar<-A%*%r
    
    #lamda.dash
    # length(lamda.dash)
    part1=-(1/2)
    # r<-y-x%*%beta
    part2=0.5*exp(-omega.dash)*(Ar)^2*exp(-lamda.dash)
    div_sem_h_lamda<-part1+part2
    
    div_plamda_dlamda<--(df/2)+(df/2)*exp(-lamda.dash)
    div_h_lamda<-div_sem_h_lamda+div_plamda_dlamda # 1/lamda is the jec adjusment
    
    #beta
    Sigma_lamda_inv=Diagonal(n,(1/lamda))
    M=t(A)%*%Sigma_lamda_inv%*%A
    div_h_beta<-exp(-omega.dash)*t(r)%*%M%*%x-t(beta)/omega_beta
    
    
    #omega.dash
    div_h_omega.dash<--0.5*n+(0.5*exp(-omega.dash))*t(Ar)%*%Sigma_lamda_inv%*%(Ar)-omega.dash/omega_omega.dash
    
    #rho.dash
    ATA=I-rho*wpluswt+rho*rho*wtw
    drho_by_drho.dash<--((rho.a-rho.b)*exp(rho.dash))/(1+exp(rho.dash))^2
    ATA_by_rho.dash<-(-wpluswt+2*rho*wtw)*drho_by_drho.dash
    
    #
    # dlogdetATA_by_rho.dash<-sum(diag(solve(Cholesky(ATA),ATA_by_rho.dash))) #For sim
    dlogdetATA_by_rho.dash<-sum(diag(solve((ATA),ATA_by_rho.dash))) # For real
    
    der.of.quadratic<--0.5*exp(-omega.dash)*t(r)%*%(-Sigma_lamda_inv%*%w-t(w)%*%Sigma_lamda_inv+2*rho*t(w)%*%Sigma_lamda_inv%*%w)%*%r*drho_by_drho.dash
    
    div_h_drho.dash<-0.5*dlogdetATA_by_rho.dash+der.of.quadratic-rho.dash/omega_rho.dash
    
    
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
    
    #Psi
    lincomb<-cbind(xm,y)
    div_h_psi<-(t(m-(exp(lincomb%*%psi)/((1+exp(lincomb%*%psi)))))%*%lincomb)-t(psi)/omega_psi
    
    fullgrad<-c(as.vector(div_h_lamda),as.vector(div_h_beta),
                as.vector(div_h_omega.dash),as.vector(div_h_drho.dash),as.vector(div_h_ddf_dash),as.vector(div_h_psi))
    
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  
  sim_yugyotheta_<-function(x,xu,xo,xm,w,wpluswt,wtw,yo,# efficiently calculate con mean and var
                            yu.start,lamda.dash,m,beta,
                            omega.dash,rho.dash,psi,N){
    nu=length(yu.start)
    lamda<-exp(lamda.dash)
    omega<-exp(omega.dash)
    sigma2y<-omega
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    
    lamda.list<-make.lamda.partion(lamda,k,reminder)
    
    cal.con.dis<-function(x,yo,yui_,w,I,lamda,beta,rho,sigma2y,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      Sigma_lamda_inv=Diagonal(n,(1/lamda))
      A<-I-rho*w
      M<-t(A)%*%Sigma_lamda_inv%*%A
      # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
      if(reminder==0){
        # Sigma_ui.ui_o<-Sigma[1:bsize,(bsize+1):n] #between cov(yui and yui_,yo)
        # Sigma_ui_o.ui_o<-Sigma[(bsize+1):n,(bsize+1):n] #cov(yui_,yo and yui_,yo)
        # Sigma_uiui<-Sigma[1:bsize,1:bsize]
        # 
        # 
        # mu_uigui_o<-x[1:bsize,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
        #                                                      c(yui_,yo)-(x[(bsize+1):n,])%*%beta)
        # sigma_uigui_o<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
        
        
        # sigma_uj.v<-Sigma[1:bsize,(bsize+1):n]
        # sigma_v.v<-Sigma[(bsize+1):n,(bsize+1):n] 
        # sigma_uj.uj<-Sigma[1:bsize,1:bsize] 
        # mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        # v<-c(yui_,yo)
        # 
        # mu_ujguv<-x[1:bsize,]%*%beta+sigma_uj.v%*%solve(sigma_v.v,(v-mu_v))
        # sigma_ujgv<-sigma_uj.uj-sigma_uj.v%*%solve(sigma_v.v,t(sigma_uj.v))
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
          # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
          # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
          # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
          # 
          # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
          #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
          # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
          
          ################################
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
          
        }
        
      }
      # print(lamda[1])
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yu,yo)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    yu.chain[1,]<-yu.start
    current_state<-yu.start
    
    # tic()
    accepts<-0
    for(i in 2:N1){#i=3 N=5
      
      
      # j=1
      current<-0
      for (j in 1:(k+ifelse(reminder>0,1,0))) {
        # accepts<-0
        # yubj<-current_state[(current+1):(current+bsize)]
        if(j<=k){ # for full blocks j=2
          yubj_<-current_state[-((current+1):(current+bsize))]
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
          yubj_<-current_state[-((current+1):(current+reminder))]
          # current=10
          # print("in half block")
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          # print("out half block")
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  start.df.dash<-start_theta[length(start_theta)-ncol(xm)-1] # additional -1 is for psi_y
  start.df<-exp(start.df.dash)+3
  mu<-c(log(extraDistr::rinvgamma(n,start.df/2,start.df/2)),start_theta) # mu=(start_lamda,start_theta except rho) 
  # startyu<-yu
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
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) #To store lamda in each iteration
  tlowerbound<-c()
  accepts<-c()
  #i=1
  for(i in 1:N){
    
    # thetas
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    lamda_theta<-mu+B%*%z+d*as.vector(epsilon)
    D<-diag(d)
    
    # generate yus from p(yu | yo, theta)
    lamda.dash<-lamda_theta[1:n]
    # lamda_theta[-(1:n)]
    beta<-lamda_theta[(n+1):(n+numcov)]
    omega.dash<-lamda_theta[(n+numcov+1)] 
    rho.dash<-lamda_theta[(n+numcov+2)] #transform of rho 
    df.dash<-lamda_theta[(n+numcov+3)]
    psi<-lamda_theta[(n+numcov+4):length(lamda_theta)]
    
    # N1=100
    
    mysamples<-sim_yugyotheta_(x=x,xu=xu,xo=xo,xm=xm,w=w,wpluswt = wpluswt,wtw =wtw ,yo=yo,
                               yu.start=startyu,lamda.dash,m=m,beta=beta,
                               omega.dash=omega.dash,rho.dash=rho.dash,psi=psi,N=N1)
    accepts[i]<-mysamples$accepts
    yunote<-(mysamples$yugyom.chain)[N1,]
    startyu<-yunote
    
    #############
    # accepts[i]<-1
    # yunote<-yu
    # startyu<-yu 
    ###########
    yu_lamda_theta<-c(as.vector(yunote),as.vector(lamda_theta)) # combine thetas and yu
    
    # tic()
    gradg_theta<-grad_2_new(yu_lamda_theta=yu_lamda_theta,w=w,x,I=I,yo=yo,xo=xo,xu=xu,xm=xm,
                            wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n,m=m) # gradient of log h(theta)
    
    # toc()
    
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    #Construct UB estimates
    grad_mu<-gradg_theta-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradg_theta))
    grad_B<-t(t(z)%x%I_theta)%*%(gradg_theta-dlogq_by_dtheta)
    gradd<-epsilon*(gradg_theta-dlogq_by_dtheta)
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
    
    # calculate lowerbound
    
    
    # tlowerbound[i]=system.time(Lbound[i]<-cal_Lbound(theta_u=theta,nu=nu,mu=mu,d=d,B=B,m=m,
    # x=x,yo=yo,wpluswt=wpluswt,wtw=wtw))[[3]]
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    
    
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
    # print(i)
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,
              all_paras=all_paras,accepts=accepts,
              dis_yug_yo=mysamples[[3]]))
  
  
}



YJ_SEM_t_MNAR_Buall_noninfo<-function(x,xm,m,w,yo,yu.start,p, #  Code is correct
                                      start.theta,N,N1,bsize,shape=4,rate=1/2){ 
  # p=4
  rho.a<--1
  rho.b<-1
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rcm=ncol(xm)+1
  rc=ncol(x) #number of regression covariates. i.e. betas
  numcov<-rc
  n<-nrow(x)
  totpara<-n+rc+4+rcm
  adapt_epsilon=10^(-6)
  v=0.99
  
  
  # bsize=1000
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  
  
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
  
  f<-function(df.dash){(digamma(0.5*exp(df.dash)+1.5))}
  f2<-function(df.dash){(lgamma(0.5*exp(df.dash)+1.5))}
  
  calgrad<-function(yu.lamda.theta,wpluswt,wtw,w,x,xm,I,yo,mu,n){
    # m<-ncol(x)
    numcov<-ncol(x) # Here p mean dimension of X.
    # print(sum(xm[,2]))
    yu<-yu.lamda.theta[(1:nu)]
    lamda_theta<-yu.lamda.theta[-(1:nu)]
    lamda.dash<-lamda_theta[1:n]
    beta<-lamda_theta[(n+1):(n+numcov)]
    omega.dash<-lamda_theta[(n+numcov+1)] #transform of sigma2
    rho.dash<-lamda_theta[(n+numcov+2)] #transform of sigma rho
    df.dash<-lamda_theta[(n+numcov+3)]
    gama.dash<-lamda_theta[(n+numcov+4)]
    psi<-lamda_theta[(n+numcov+5):length(lamda_theta)]
    
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
    omega_psi<-10000
    
    y<-c(yo,yu)
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
    dlogdetATA_by_rho.dash<-sum(diag(solve(Cholesky(ATA),ATA_by_rho.dash)))
    
    der.of.quadratic<--0.5*exp(-omega.dash)*t(r)%*%(-Sigma_lamda_inv%*%w-t(w)%*%Sigma_lamda_inv+2*rho*t(w)%*%Sigma_lamda_inv%*%w)%*%r*drho_by_drho.dash
    
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
    
    #psi
    zmiss<-cbind(xm,y)
    # div_h_psi_<-(t(m-(exp(zmiss%*%psi)/((1+exp(zmiss%*%psi)))))%*%zmiss)-t(psi)/omega_psi
    div_h_psi<-t(zmiss)%*%(m-(1/(1+exp(-zmiss%*%psi))))-matrix(t(psi)/omega_psi,ncol=1)
    # dim(t(zmiss))
    # message("new",paste(div_h_psi))
    # message("old",paste(div_h_psi_))
    
    fullgrad<-c(as.vector(div_h_lamda),as.vector(div_h_beta),
                as.vector(div_h_omega.dash),as.vector(div_h_drho.dash),as.vector(div_h_ddf_dash),as.vector(div_h_gama.dash),
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
  
  
  make.lamda.partion<-function(lamda,k,reminder){
    
    if(reminder==0){
      
      lamda.list<-list()
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        lamda.list[[i]]<-lamda[c(ui,ui_,o)]
        
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        lamda.list<-list()
        
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(lamda.list=lamda.list))
  }
  
  
  #YJ-SEM-t-NoB
  
  # 
  sim_yu.g.yom.theta.uallb<-function(yo,listw_x_m,    # After adding un-balanced # of blocks 
                                     start.yu,lamda.dash,beta,
                                     omega.dash,gama.dash,rho.dash,psi,
                                     N1,k){
    
    lamda<-exp(lamda.dash)
    omega<-exp(omega.dash)
    sigma2y<-omega
    rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    gama<-2*exp(gama.dash)/(exp(gama.dash)+1)
    lamda.list<-make.lamda.partion(lamda,k,reminder)
    # df<-exp(df.dash)+3
    
    nu=length(start.yu)
    zo<-cal_z(yo,gama)
    start.zu<-cal_z(start.yu,gama)
    
    cal.con.dis<-function(x,zo,zui_,w,I,lamda,beta,rho,sigma2y,bsize,
                          no,nu,reminder,j,k){ #j is the iteration number
      
      Sigma_lamda_inv=Diagonal(n,(1/lamda))
      A<-I-rho*w
      M<-t(A)%*%Sigma_lamda_inv%*%A
      # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
      if(reminder==0){
        
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(zui_,zo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        # print(isSymmetric(M_uj.uj)) #
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          # print(isSymmetric(M_uj.uj)) #
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
          # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
          # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
          # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
          # 
          # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
          #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
          # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
          
          ################################
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yo,yu)
      lincomb=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-lincomb))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    
    # zu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    # zu.chain[1,]<-start.zu
    # current_state<-start.zu
    
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
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
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
          zubj_<-cal_z(current_state[-((current+1):(current+reminder))],gama)
          # current=10
          # print("in half block")
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
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
        # print(ration)
      }
      
      
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder=reminder)
  
  #initial values for lamdha
  #p=2
  start.df.dash<-start.theta[length(start.theta)-1]
  start.df<-exp(start.df.dash)+3
  start.lamda<-rinvgamma(n,start.df/2,2/start.df)
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
  accepts<-c()
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    eta<-matrix(zeta[1:p],ncol = 1)
    
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    lamda_theta<-mu+B%*%eta+d*as.vector(epsilon)
    # theta_t<-matrix(theta_t,ncol=1)
    D<-diag(d)
    
    ############# Sampling of yu | yo,m, lambda
    lamda.dash<-lamda_theta[1:n]
    beta<-lamda_theta[(n+1):(n+numcov)]
    omega.dash<-lamda_theta[(n+numcov+1)] #transform of sigma2
    rho.dash<-lamda_theta[(n+numcov+2)] #transform of sigma rho
    df.dash<-lamda_theta[(n+numcov+3)]
    gama.dash<-lamda_theta[(n+numcov+4)]
    psi<-lamda_theta[(n+numcov+5):length(lamda_theta)]
    
    # message(paste("psi:", psi))
    #######################
    # N1=1000
    # N=100
    mysamples<-sim_yu.g.yom.theta.uallb(yo,listw_x_m,    # After adding un-balanced # of blocks
                                        start.yu=yu.start,lamda.dash,beta,
                                        omega.dash,gama.dash,rho.dash,psi,
                                        N1,k)
    
    accepts[i]<-mysamples$accepts
    yunote<-(mysamples$yugyom.chain)[N1,]
    yu.start<-yunote
    yu.lamda.theta<-c(as.vector(yunote),as.vector(lamda_theta)) # combine thetas and yu
    con.mean<-mysamples$co.dis[[1]]
    # 
    #######################
    # yunote<-yu.start
    # yu.start<-yu.start
    # yu.lamda.theta<-c(as.vector(yunote),as.vector(lamda_theta)) # combine thetas and yu
    # con.mean<-yu.start
    
    ############################## # Calculate log gradient of theta
    # sum((density(yu.start))$x)
    # print(gradiant[(length(gradiant)-2):length(gradiant)])
    
    gradiant<-calgrad(yu.lamda.theta=yu.lamda.theta,wpluswt,wtw,w,x,xm,I,yo,mu,n)
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
    # print(psi)
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,accepts=accepts,con.mean=con.mean,all_paras=all_paras))
} 



################## Function to draw samples from the posterior (variational) distribution of parameters




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
  
  
  
  
} #SEM-Gau 

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
  
  
  
  
} #  YJ-SEM-Gau

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


################## Function to draw samples from the posterior distribution of missing values




generate_posterior_yu_SEMGau<-function(theta,yo,start.yu,x,xm,w,m,bsize,N1){
  
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  
  
  n<-length(m)
  no<-length(yo)
  nu<-n-no
  
  
  # bsize=10
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  psi<-theta[(p+3):length(theta)]
  
  # rho<-(exp(lamda)-1)/(exp(lamda)+1)
  # sigma2y<-exp(gama)
  
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
  
  sim_yugyotheta_<-function(yo,listw_x_m,    # After adding un-balanced # of blocks 
                            yu.start,beta,
                            sigma2y,rho,psi,
                            N1,k){
    nu=length(yu.start)
    
    
    cal.con.dis<-function(x,yo,yui_,w,I,beta,rho,sigma2y,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      M<-t(I-rho*w)%*%(I-rho*w)
      # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
      if(reminder==0){
        # Sigma_ui.ui_o<-Sigma[1:bsize,(bsize+1):n] #between cov(yui and yui_,yo)
        # Sigma_ui_o.ui_o<-Sigma[(bsize+1):n,(bsize+1):n] #cov(yui_,yo and yui_,yo)
        # Sigma_uiui<-Sigma[1:bsize,1:bsize]
        # 
        # 
        # mu_uigui_o<-x[1:bsize,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
        #                                                      c(yui_,yo)-(x[(bsize+1):n,])%*%beta)
        # sigma_uigui_o<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
        
        
        # sigma_uj.v<-Sigma[1:bsize,(bsize+1):n]
        # sigma_v.v<-Sigma[(bsize+1):n,(bsize+1):n] 
        # sigma_uj.uj<-Sigma[1:bsize,1:bsize] 
        # mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        # v<-c(yui_,yo)
        # 
        # mu_ujguv<-x[1:bsize,]%*%beta+sigma_uj.v%*%solve(sigma_v.v,(v-mu_v))
        # sigma_ujgv<-sigma_uj.uj-sigma_uj.v%*%solve(sigma_v.v,t(sigma_uj.v))
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
          # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
          # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
          # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
          # 
          # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
          #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
          # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
          
          ################################
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yu,yo)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    yu.chain[1,]<-yu.start
    current_state<-yu.start
    
    # tic()
    accepts<-0
    for(i in 2:N1){#i=3 N=5
      
      
      # j=1
      current<-0
      for (j in 1:(k+ifelse(reminder>0,1,0))) {
        # accepts<-0
        # yubj<-current_state[(current+1):(current+bsize)]
        if(j<=k){ # for full blocks j=2
          yubj_<-current_state[-((current+1):(current+bsize))]
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
          yubj_<-current_state[-((current+1):(current+reminder))]
          # current=10
          # print("in half block")
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          # print("out half block")
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  
  mysamples<-sim_yugyotheta_(yo=yo,listw_x_m=listw_x_m,
                             yu.start=start.yu,beta=beta,
                             sigma2y=sigma2y,rho=rho,
                             psi=psi,N1=N1,k=k)
  
  return(apply((mysamples$yugyom.chain)[((N1/2):N1),],2,mean))
  
} # SEM-Gau 

generate_posterior_yu_YJSEMGau<-function(theta,yo,start.yu,x,xm,w,m,bsize,N1){
  
  
  
  
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
  
} # YJ-SEM-Gau

generate_posterior_yu_SEMt<-function(theta,yo,start.yu,x,xm,w,m,bsize,N1){
  
  # theta=SEM.t.theta.sample[1,]
  rho.a<--1
  rho.b<-1
  # p=10
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  numcov<-ncol(x)
  xm<-cbind(1,xm)
  n<-length(m)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=ncol(xm)+1     #number of logistic covariates. i.e. betas
  
  # start.yu
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  
  no<-length(yo)
  nu<-n-no
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  f<-function(df.dash){(digamma(0.5*exp(df.dash)+1.5))}
  f2<-function(df.dash){(lgamma(0.5*exp(df.dash)+1.5))}
  
  
  
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
  
  make.lamda.partion<-function(lamda,k,reminder){
    
    if(reminder==0){
      
      lamda.list<-list()
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        lamda.list[[i]]<-lamda[c(ui,ui_,o)]
        
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        lamda.list<-list()
        
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(lamda.list=lamda.list))
  }
  
  
  sim_yugyotheta_<-function(x,xu,xo,xm,w,wpluswt,wtw,yo,# efficiently calculate con mean and var
                            yu.start,lamda.dash,m,beta,
                            sigma2y,rho,psi,N){
    nu=length(yu.start)
    lamda<-exp(lamda.dash)
    # omega<-exp(omega.dash)
    # sigma2y<-omega
    # rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    
    lamda.list<-make.lamda.partion(lamda,k,reminder)
    
    cal.con.dis<-function(x,yo,yui_,w,I,lamda,beta,rho,sigma2y,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      Sigma_lamda_inv=Diagonal(n,(1/lamda))
      A<-I-rho*w
      M<-t(A)%*%Sigma_lamda_inv%*%A
      # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
      if(reminder==0){
        # Sigma_ui.ui_o<-Sigma[1:bsize,(bsize+1):n] #between cov(yui and yui_,yo)
        # Sigma_ui_o.ui_o<-Sigma[(bsize+1):n,(bsize+1):n] #cov(yui_,yo and yui_,yo)
        # Sigma_uiui<-Sigma[1:bsize,1:bsize]
        # 
        # 
        # mu_uigui_o<-x[1:bsize,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
        #                                                      c(yui_,yo)-(x[(bsize+1):n,])%*%beta)
        # sigma_uigui_o<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
        
        
        # sigma_uj.v<-Sigma[1:bsize,(bsize+1):n]
        # sigma_v.v<-Sigma[(bsize+1):n,(bsize+1):n] 
        # sigma_uj.uj<-Sigma[1:bsize,1:bsize] 
        # mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        # v<-c(yui_,yo)
        # 
        # mu_ujguv<-x[1:bsize,]%*%beta+sigma_uj.v%*%solve(sigma_v.v,(v-mu_v))
        # sigma_ujgv<-sigma_uj.uj-sigma_uj.v%*%solve(sigma_v.v,t(sigma_uj.v))
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
          # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
          # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
          # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
          # 
          # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
          #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
          # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
          
          ################################
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
          
        }
        
      }
      # print(lamda[1])
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yu,yo)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    yu.chain[1,]<-yu.start
    current_state<-yu.start
    
    # tic()
    accepts<-0
    for(i in 2:N1){#i=3 N=5
      
      
      # j=1
      current<-0
      for (j in 1:(k+ifelse(reminder>0,1,0))) {
        # accepts<-0
        # yubj<-current_state[(current+1):(current+bsize)]
        if(j<=k){ # for full blocks j=2
          yubj_<-current_state[-((current+1):(current+bsize))]
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
          yubj_<-current_state[-((current+1):(current+reminder))]
          # current=10
          # print("in half block")
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          # print("out half block")
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
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
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  ##
  numcov<-ncol(x)
  
  # lamda.dash<-yu_lamda_theta[(nu+1):(n+nu)]
  beta<-theta[(1):(numcov)]
  omega.dash<-theta[(numcov+1)] #transform of sigma2 
  rho.dash<-theta[(numcov+2)] #transform of rho 
  df.dash<-theta[(numcov+3)]
  psi<-theta[(numcov+4):length(theta)]
  
  # Dont need to do tranfomation, as theta is in the original scale
  omega<-(omega.dash)
  rho<-rho.dash
  
  df<-df.dash
  lamda<-rinvgamma(n,df/2,df/2)
  lamda.dash<-log(lamda)
  ##
  
  
  mysamples<-sim_yugyotheta_(x=x,xu=xu,xo=xo,xm=xm,w=w,wpluswt = wpluswt,wtw =wtw ,yo=yo,
                             yu.start=start.yu,lamda.dash,m=m,beta=beta,
                             sigma2y=omega,rho=rho,psi=psi,N=N1)
  # accepts[i]<-mysamples$accepts
  # yunote<-(mysamples$yugyom.chain)[N1,]
  # startyu<-yunote
  
  return(apply((mysamples$yugyom.chain)[((N1/2):N1),],2,mean))
  
} # SEM-t

generate_posterior_yu_YJSEMt<-function(theta,yo,start.yu,x,xm,w,m,bsize,N1){
  
  # p=4
  rho.a<--1
  rho.b<-1
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  
  # bsize=1000
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  
  n<-length(m)
  I<-Diagonal(n)
  # I_p<-Diagonal(p)
  
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
  
  
  make.lamda.partion<-function(lamda,k,reminder){
    
    if(reminder==0){
      
      lamda.list<-list()
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        lamda.list[[i]]<-lamda[c(ui,ui_,o)]
        
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        lamda.list<-list()
        
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            lamda.list[[i]]<-lamda[c(ui,ui_,o)]
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(lamda.list=lamda.list))
  }
  
  
  #YJ-SEM-t-NoB
  
  # 
  sim_yu.g.yom.theta.uallb<-function(yo,listw_x_m,    # After adding un-balanced # of blocks 
                                     start.yu,lamda,beta,
                                     omega,gama,rho,psi,
                                     N1,k){
    
    # lamda<-exp(lamda.dash)
    # omega<-exp(omega.dash)
    sigma2y<-omega
    # rho<-(exp(rho.dash)-1)/(exp(rho.dash)+1)
    
    lamda.list<-make.lamda.partion(lamda,k,reminder)
    # df<-exp(df.dash)+3
    
    nu=length(start.yu)
    zo<-cal_z(yo,gama)
    start.zu<-cal_z(start.yu,gama)
    
    cal.con.dis<-function(x,zo,zui_,w,I,lamda,beta,rho,sigma2y,bsize,
                          no,nu,reminder,j,k){ #j is the iteration number
      
      Sigma_lamda_inv=Diagonal(n,(1/lamda))
      A<-I-rho*w
      M<-t(A)%*%Sigma_lamda_inv%*%A
      # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
      if(reminder==0){
        
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(zui_,zo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        # print(isSymmetric(M_uj.uj)) #
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          # print(isSymmetric(M_uj.uj)) #
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
          # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
          # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
          # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
          # 
          # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
          #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
          # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
          
          ################################
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(zui_,zo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yo,yu)
      lincomb=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-lincomb))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    
    # zu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    # zu.chain[1,]<-start.zu
    # current_state<-start.zu
    
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
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
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
          zubj_<-cal_z(current_state[-((current+1):(current+reminder))],gama)
          # current=10
          # print("in half block")
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],zo=zo,zui_=zubj_,w=listw_x_m$w.list[[j]],I=I,lamda.list$lamda.list[[j]],beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
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
        # print(ration)
      }
      
      
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder=reminder)
  
  
  numcov<-ncol(x)
  # lamda.dash<-lamda_theta[1:n]
  beta<-theta[(1):(numcov)]
  omega<-theta[(numcov+1)] #transform of sigma2
  rho<-theta[(numcov+2)] #transform of sigma rho
  df<-theta[(numcov+3)]
  gama<-theta[(numcov+4)]
  psi<-theta[(numcov+5):length(theta)]
  lamda<-rinvgamma(n,df/2,df/2)
  
  
  
  mysamples<-sim_yu.g.yom.theta.uallb(yo,listw_x_m,    # After adding un-balanced # of blocks
                                      start.yu=start.yu,lamda,beta,
                                      omega,gama,rho,psi,
                                      N1,k)
  
  
  
  
  
  
  
  
  
  
  return(apply((mysamples$yugyom.chain)[((N1/2):N1),],2,mean))
  
} # YJ-SEM-t



################ Functions to compare model fit using DIC5 and DIC7 




cal_DIC_5_7_MNAR_SEMGau<-function(SEM.Gau.theta.sample,SEM.Gau.post.sample,
                                  yo,x,xm,w,m){
  SEM.Gau.post.sample<-t(SEM.Gau.post.sample)
  yu.theta.sample<-cbind(SEM.Gau.post.sample,SEM.Gau.theta.sample)
  xm<-cbind(1,xm)
  
  
  cal.complete.log.density1<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    # gama<-theta[ncol(cbind(1,x))+3]
    psi<-theta[(ncol(cbind(1,x))+3):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    
    loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(complete.log.density))
  }
  
  # tic()
  complete.log.densities1<-apply(yu.theta.sample, 1, cal.complete.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(complete.log.densities1)
  
  
  
  #DIC5
  cal.comp.log.post<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    # gama<-theta[ncol(cbind(1,x))+3]
    psi<-theta[(ncol(cbind(1,x))+3):length(theta)]
    
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(sigma2)
    rho.dash<-log(1+rho)-log(1-rho)
    # gama.dash<-log(gama/(2-gama))
    
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    omega_psi=10000
    ##################
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    
    loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    
    complete.log.density<-loglike1+loglike2
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(rho.dash,0,omega_rho.dash,log=T)+
      mvnfast::dmvn(psi,rep(0,length(psi)),omega_beta*diag(length(psi)),log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(complete.log.density+log.prior))
  }
  
  log.comp.posts<-apply(yu.theta.sample, 1, cal.comp.log.post,yo,x,xm,w,m)
  # toc()
  theta<-SEM.Gau.theta.sample[which.max(log.comp.posts),]# find the theta that max post of theta. p(theta | y)
  yu<-SEM.Gau.post.sample[which.max(log.comp.posts),]
  y<-c(yo,yu)
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  # gama<-theta[p+3]
  psi<-theta[(p+3):length(theta)]
  n<-ncol(w)
  I<-Diagonal(n)
  
  A<-I-rho*w
  M<-t(A)%*%A
  Sigma_y<-sigma2y*solve(M,I)
  
  lin.comb=cbind(xm,y)%*%psi
  pxy<-1/(1+exp(-lin.comb))
  loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
  
  part2<-2*(mvnfast::dmvn(y,as.vector(x%*%beta),as.matrix(Sigma_y),log = T)+
              loglike1)
  
  
  
  DIC5<-as.numeric(part1+part2)
  print("DIC 5 done")
  
  ### DIC 7
  cal.condi.log.density1<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(x)]
    sigma2<-theta[ncol(x)+1]
    rho<-theta[ncol(x)+2]
    # gama<-theta[ncol(x)+3]
    psi<-theta[(ncol(x)+3):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-x%*%betas
    r<-y-x.beta
    
    loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    no<-length(yo)
    M_oo=M[1:no,1:no]
    M_uu=M[(no+1):n,(no+1):n]
    M_ou=M[1:no,(no+1):n]
    M_uo=t(M_ou)
    B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
    xu=x[-(1:no),]
    nu<-n-no
    ru<-r[-(1:no)]
    
    loglike3<--0.5*nu*log(2*pi)-0.5*nu*log(sigma2)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(ru)%*%B.uu%*%ru
    
    condi.log.density<-loglike1+loglike2-loglike3
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(condi.log.density))
    
  }
  
  condi.log.densities1<-apply(yu.theta.sample, 1, cal.condi.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(condi.log.densities1)
  # for part 2, the part 2 from DIC5, - 2*log(yu|theta).
  
  no<-length(yo)
  M_oo=M[1:no,1:no]
  M_uu=M[(no+1):n,(no+1):n]
  M_ou=M[1:no,(no+1):n]
  M_uo=t(M_ou)
  B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
  xu=x[-(1:no),]
  nu<-n-no
  x.beta<-x%*%beta
  r<-y-x.beta
  ru<-r[-(1:no)]
  
  loglike3<--0.5*nu*log(2*pi)-0.5*nu*log(sigma2y)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-
    0.5*(1/sigma2y)*t(ru)%*%B.uu%*%ru
  
  part2<-part2-2*loglike3 # for part 2 for DIC7=e part 2 from DIC5, - 2*log(yu|theta).
  DIC7<-as.numeric(part1+part2)
  
  return(list(DIC5=DIC5,DIC7=DIC7))
  
}


cal_DIC_5_7_MNAR_YJSEMGau<-function(YJ.SEM.Gau.theta.sample,YJ.SEM.Gau.post.sample,
                                    yo,x,xm,w,m){
  YJ.SEM.Gau.post.sample<-t(YJ.SEM.Gau.post.sample)
  yu.theta.sample<-cbind(YJ.SEM.Gau.post.sample,YJ.SEM.Gau.theta.sample)
  xm<-cbind(1,xm)
  
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
  
  
  cal.complete.log.density1<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    gama<-theta[ncol(cbind(1,x))+3]
    psi<-theta[(ncol(cbind(1,x))+4):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    
    loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(complete.log.density))
  }
  
  # tic()
  complete.log.densities1<-apply(yu.theta.sample, 1, cal.complete.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(complete.log.densities1)
  
  
  
  #DIC5
  cal.comp.log.post<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    gama<-theta[ncol(cbind(1,x))+3]
    psi<-theta[(ncol(cbind(1,x))+4):length(theta)]
    
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(sigma2)
    rho.dash<-log(1+rho)-log(1-rho)
    gama.dash<-log(gama/(2-gama))
    
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    omega_psi=10000
    ##################
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    
    loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    
    complete.log.density<-loglike1+loglike2
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(rho.dash,0,omega_rho.dash,log=T)+
      dnorm(gama.dash,0,omega_rho.dash,log=T)+
      mvnfast::dmvn(psi,rep(0,length(psi)),omega_beta*diag(length(psi)),log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(complete.log.density+log.prior))
  }
  
  log.comp.posts<-apply(yu.theta.sample, 1, cal.comp.log.post,yo,x,xm,w,m)
  # toc()
  theta<-YJ.SEM.Gau.theta.sample[which.max(log.comp.posts),]# find the theta that max post of theta. p(theta | y)
  yu<-YJ.SEM.Gau.post.sample[which.max(log.comp.posts),]
  y<-c(yo,yu)
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  gama<-theta[p+3]
  psi<-theta[(p+4):length(theta)]
  n<-ncol(w)
  I<-Diagonal(n)
  
  A<-I-rho*w
  M<-t(A)%*%A
  Sigma_y<-sigma2y*solve(M,I)
  
  lin.comb=cbind(xm,y)%*%psi
  pxy<-1/(1+exp(-lin.comb))
  loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
  
  part2<-2*(mvnfast::dmvn(cal_z(y,gama),as.vector(x%*%beta),as.matrix(Sigma_y),log = T)+
              sum(log(cal_der_tYJ_wrt_y(y,gama)))+loglike1)
  
  
  
  DIC5<-as.numeric(part1+part2)
  
  ### DIC 7
  cal.condi.log.density1<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(x)]
    sigma2<-theta[ncol(x)+1]
    rho<-theta[ncol(x)+2]
    gama<-theta[ncol(x)+3]
    psi<-theta[(ncol(x)+4):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-x%*%betas
    r<-cal_z(y,gama)-x.beta
    
    loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(r)%*%M%*%r+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    no<-length(yo)
    M_oo=M[1:no,1:no]
    M_uu=M[(no+1):n,(no+1):n]
    M_ou=M[1:no,(no+1):n]
    M_uo=t(M_ou)
    B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
    xu=x[-(1:no),]
    nu<-n-no
    ru<-r[-(1:no)]
    
    loglike3<--0.5*nu*log(2*pi)-0.5*nu*log(sigma2)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-
      0.5*(1/sigma2)*t(ru)%*%B.uu%*%ru+sum(log(cal_der_tYJ_wrt_y(yu,gama)))
    
    condi.log.density<-loglike1+loglike2-loglike3
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(condi.log.density))
    
  }
  
  condi.log.densities1<-apply(yu.theta.sample, 1, cal.condi.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(condi.log.densities1)
  # for part 2, the part 2 from DIC5, - 2*log(yu|theta).
  
  no<-length(yo)
  M_oo=M[1:no,1:no]
  M_uu=M[(no+1):n,(no+1):n]
  M_ou=M[1:no,(no+1):n]
  M_uo=t(M_ou)
  B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
  xu=x[-(1:no),]
  nu<-n-no
  x.beta<-x%*%beta
  r<-cal_z(y,gama)-x.beta
  ru<-r[-(1:no)]
  
  loglike3<--0.5*nu*log(2*pi)-0.5*nu*log(sigma2y)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-
    0.5*(1/sigma2y)*t(ru)%*%B.uu%*%ru+sum(log(cal_der_tYJ_wrt_y(yu,gama)))
  
  part2<-part2-2*loglike3 # for part 2 for DIC7=e part 2 from DIC5, - 2*log(yu|theta).
  DIC7<-as.numeric(part1+part2)
  
  return(list(DIC5=DIC5,DIC7=DIC7))
  
}


cal_DIC_5_7_MNAR_SEMt<-function(SEM.t.theta.sample,SEM.t.post.sample,
                                yo,x,xm,w,m){
  SEM.t.post.sample<-t(SEM.t.post.sample)
  yu.theta.sample<-cbind(SEM.t.post.sample,SEM.t.theta.sample)
  xm<-cbind(1,xm)
  
  
  cal.complete.log.density1<-function(yu.theta,yo,x,xm,w,m){
    # yu.theta<-yu.theta.sample[1,]
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    df<-theta[ncol(cbind(1,x))+3]
    psi<-theta[(ncol(cbind(1,x))+4):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    
    # loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
    # 0.5*(1/sigma2)*t(r)%*%M%*%r
    
    ######################
    # loglike2<-log(gamma(0.5*(df+n)))-log(gamma(0.5*(df)))-0.5*n*log(sigma2)
    
    
    loglike2<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)
    
    # dmvt(y, x.beta, sigma2*solve(M), df=df, log=T)
    
    ###########################
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(complete.log.density))
  }
  
  # tic()
  complete.log.densities1<-apply(yu.theta.sample, 1, cal.complete.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(complete.log.densities1)
  
  
  
  #DIC5
  cal.comp.log.post<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    df<-theta[ncol(cbind(1,x))+3]
    psi<-theta[(ncol(cbind(1,x))+4):length(theta)]
    
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(sigma2)
    rho.dash<-log(1+rho)-log(1-rho)
    # gama.dash<-log(gama/(2-gama))
    df.dash<-log(df-3)
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    omega_psi=10000
    ##################
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-y-x.beta
    
    loglike2<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    
    complete.log.density<-loglike1+loglike2
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(rho.dash,0,omega_rho.dash,log=T)+
      dnorm(df.dash,0,omega_omega.dash,log=T)+
      mvnfast::dmvn(psi,rep(0,length(psi)),omega_beta*diag(length(psi)),log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(complete.log.density+log.prior))
  }
  
  log.comp.posts<-apply(yu.theta.sample, 1, cal.comp.log.post,yo,x,xm,w,m)
  # toc()
  theta<-SEM.t.theta.sample[which.max(log.comp.posts),]# find the theta that max post of theta. p(theta | y)
  yu<-SEM.t.post.sample[which.max(log.comp.posts),]
  y<-c(yo,yu)
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  df<-theta[p+3]
  psi<-theta[(p+4):length(theta)]
  n<-ncol(w)
  I<-Diagonal(n)
  
  A<-I-rho*w
  M<-t(A)%*%A
  x.beta<-x%*%beta
  r<-y-x.beta
  
  lin.comb=cbind(xm,y)%*%psi
  pxy<-1/(1+exp(-lin.comb))
  loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
  
  part2<-2*(lgamma((df+n)*0.5)-0.5*n*log(sigma2y)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
              lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2y))*t(r)%*%M%*%r)+
              loglike1)
  
  
  
  DIC5<-as.numeric(part1+part2)
  print("DIC 5 done")
  
  ### DIC 7
  cal.condi.log.density1<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(x)]
    sigma2<-theta[ncol(x)+1]
    rho<-theta[ncol(x)+2]
    df<-theta[ncol(x)+3]
    psi<-theta[(ncol(x)+4):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-x%*%betas
    r<-y-x.beta
    
    loglike2<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    no<-length(yo)
    M_oo=M[1:no,1:no]
    M_uu=M[(no+1):n,(no+1):n]
    M_ou=M[1:no,(no+1):n]
    M_uo=t(M_ou)
    B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
    xu=x[-(1:no),]
    nu<-n-no
    ru<-r[-(1:no)]
    
    loglike3<-lgamma((df+nu)*0.5)-0.5*nu*log(sigma2)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-0.5*nu*log(pi*df)-
      lgamma(df/2)-0.5*(df+nu)*log(1+(1/(df*sigma2))*t(ru)%*%B.uu%*%ru)
    
    condi.log.density<-loglike1+loglike2-loglike3
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(condi.log.density))
    
  }
  
  condi.log.densities1<-apply(yu.theta.sample, 1, cal.condi.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(condi.log.densities1)
  # for part 2, the part 2 from DIC5, - 2*log(yu|theta).
  
  no<-length(yo)
  M_oo=M[1:no,1:no]
  M_uu=M[(no+1):n,(no+1):n]
  M_ou=M[1:no,(no+1):n]
  M_uo=t(M_ou)
  B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
  xu=x[-(1:no),]
  nu<-n-no
  x.beta<-x%*%beta
  r<-y-x.beta
  ru<-r[-(1:no)]
  
  loglike3<-lgamma((df+nu)*0.5)-0.5*nu*log(sigma2y)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-0.5*nu*log(pi*df)-
    lgamma(df/2)-0.5*(df+nu)*log(1+(1/(df*sigma2y))*t(ru)%*%B.uu%*%ru)
  
  part2<-part2-2*loglike3 # for part 2 for DIC7=e part 2 from DIC5, - 2*log(yu|theta).
  DIC7<-as.numeric(part1+part2)
  
  return(list(DIC5=DIC5,DIC7=DIC7))
  
}


cal_DIC_5_7_MNAR_YJSEMt<-function(YJ.SEM.t.theta.sample,YJ.SEM.t.post.sample,
                                  yo,x,xm,w,m){
  YJ.SEM.t.post.sample<-t(YJ.SEM.t.post.sample)
  yu.theta.sample<-cbind(YJ.SEM.t.post.sample,YJ.SEM.t.theta.sample)
  xm<-cbind(1,xm)
  
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
  
  cal.complete.log.density1<-function(yu.theta,yo,x,xm,w,m){
    # yu.theta<-yu.theta.sample[1,]
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    df<-theta[ncol(cbind(1,x))+3]
    gama<-theta[ncol(cbind(1,x))+4]
    psi<-theta[(ncol(cbind(1,x))+5):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    
    # loglike2<--0.5*n*log(2*pi)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-
    # 0.5*(1/sigma2)*t(r)%*%M%*%r
    
    ######################
    # loglike2<-log(gamma(0.5*(df+n)))-log(gamma(0.5*(df)))-0.5*n*log(sigma2)
    
    
    loglike2<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    # dmvt(y, x.beta, sigma2*solve(M), df=df, log=T)
    
    ###########################
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(complete.log.density))
  }
  
  # tic()
  complete.log.densities1<-apply(yu.theta.sample, 1, cal.complete.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(complete.log.densities1)
  
  
  
  #DIC5
  cal.comp.log.post<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(cbind(1,x))]
    sigma2<-theta[ncol(cbind(1,x))+1]
    rho<-theta[ncol(cbind(1,x))+2]
    df<-theta[ncol(cbind(1,x))+3]
    gama<-theta[ncol(cbind(1,x))+4]
    psi<-theta[(ncol(cbind(1,x))+5):length(theta)]
    
    
    ########## transforms
    # theta.para.dash<-theta.para
    omega.dash<-log(sigma2)
    rho.dash<-log(1+rho)-log(1-rho)
    gama.dash<-log(gama/(2-gama))
    df.dash<-log(df-3)
    
    # omega_gama=10000
    omega_omega.dash=10000
    omega_beta=10000
    omega_rho.dash=10000
    omega_psi=10000
    ##################
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-cbind(1,x)%*%betas
    r<-cal_z(y,gama)-x.beta
    
    loglike2<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    
    complete.log.density<-loglike1+loglike2
    
    
    
    log.prior<-mvnfast::dmvn(betas,rep(0,length(betas)),omega_beta*diag(length(betas)),log=T)+
      dnorm(omega.dash,0,omega_omega.dash,log=T)+dnorm(rho.dash,0,omega_rho.dash,log=T)+
      dnorm(df.dash,0,omega_omega.dash,log=T)+dnorm(gama.dash,0,omega_omega.dash,log=T)
    mvnfast::dmvn(psi,rep(0,length(psi)),omega_beta*diag(length(psi)),log=T)
    
    # ?dmvn()
    # ?lgamma
    # isSymmetric(solve(M))
    # mvtnorm::dmvt(matrix(y,ncol=n,nrow=1), delta  = as.vector(x.beta), sigma = as.matrix(sigma2*solve(M)), df = df, log = TRUE)
    # mvnfast::dmvt(matrix(y,ncol=n,nrow=1), mu = as.vector(x.beta), sigma = sigma2*solve(M), df = df, log = TRUE)
    # ?dmvt
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.numeric(complete.log.density+log.prior))
  }
  
  log.comp.posts<-apply(yu.theta.sample, 1, cal.comp.log.post,yo,x,xm,w,m)
  # toc()
  theta<-YJ.SEM.t.theta.sample[which.max(log.comp.posts),]# find the theta that max post of theta. p(theta | y)
  yu<-YJ.SEM.t.post.sample[which.max(log.comp.posts),]
  y<-c(yo,yu)
  x<-cbind(1,x)
  p<-ncol(x) # Here p means the dimension of X.
  # yu<-theta[1:nu]
  beta<-theta[(1):(p)]
  sigma2y<-theta[(p+1)] #transform of sigma2
  rho<-theta[(p+2)] #transform of sigma rho
  df<-theta[p+3]
  gama<-theta[p+4]
  psi<-theta[(p+5):length(theta)]
  n<-ncol(w)
  I<-Diagonal(n)
  
  A<-I-rho*w
  M<-t(A)%*%A
  x.beta<-x%*%beta
  r<-cal_z(y,gama)-x.beta
  
  lin.comb=cbind(xm,y)%*%psi
  pxy<-1/(1+exp(-lin.comb))
  loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
  
  part2<-2*(lgamma((df+n)*0.5)-0.5*n*log(sigma2y)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
              lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2y))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))+
              loglike1)
  
  
  
  DIC5<-as.numeric(part1+part2)
  print("DIC 5 done")
  
  ### DIC 7
  cal.condi.log.density1<-function(yu.theta,yo,x,xm,w,m){
    
    theta<-yu.theta[-(1:nu)]
    yu<-yu.theta[(1:nu)]
    y<-c(yo,yu)
    
    betas<-theta[1:ncol(x)]
    sigma2<-theta[ncol(x)+1]
    rho<-theta[ncol(x)+2]
    df<-theta[ncol(x)+3]
    gama<-theta[ncol(x)+4]
    psi<-theta[(ncol(x)+5):length(theta)]
    
    
    n<-nrow(w)
    A<-Diagonal(n)-rho*w
    M<-t(A)%*%A
    x.beta<-x%*%betas
    r<-cal_z(y,gama)-x.beta
    
    loglike2<-lgamma((df+n)*0.5)-0.5*n*log(sigma2)+as.numeric(determinant(Cholesky(M),logarithm = T)$modulus)-0.5*n*log(pi*df)-
      lgamma(df/2)-0.5*(df+n)*log(1+(1/(df*sigma2))*t(r)%*%M%*%r)+sum(log(cal_der_tYJ_wrt_y(y,gama)))
    
    lin.comb=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-lin.comb))
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    complete.log.density<-loglike1+loglike2
    
    no<-length(yo)
    M_oo=M[1:no,1:no]
    M_uu=M[(no+1):n,(no+1):n]
    M_ou=M[1:no,(no+1):n]
    M_uo=t(M_ou)
    B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
    xu=x[-(1:no),]
    nu<-n-no
    ru<-r[-(1:no)]
    
    loglike3<-lgamma((df+nu)*0.5)-0.5*nu*log(sigma2)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-0.5*nu*log(pi*df)-
      lgamma(df/2)-0.5*(df+nu)*log(1+(1/(df*sigma2))*t(ru)%*%B.uu%*%ru)+sum(log(cal_der_tYJ_wrt_y(yu,gama)))
    
    condi.log.density<-loglike1+loglike2-loglike3
    # dmvn(y,as.vector(x.beta),sigma2*solve(M),log=T)
    return(as.vector(condi.log.density))
    
  }
  
  condi.log.densities1<-apply(yu.theta.sample, 1, cal.condi.log.density1,yo,x,xm,w,m)
  # toc()
  part1<--4*mean(condi.log.densities1)
  # for part 2, the part 2 from DIC5, - 2*log(yu|theta).
  
  no<-length(yo)
  M_oo=M[1:no,1:no]
  M_uu=M[(no+1):n,(no+1):n]
  M_ou=M[1:no,(no+1):n]
  M_uo=t(M_ou)
  B.uu<-M_uu-M_uo%*%solve(M_oo,M_ou)
  xu=x[-(1:no),]
  nu<-n-no
  x.beta<-x%*%beta
  r<-cal_z(y,gama)-x.beta
  ru<-r[-(1:no)]
  
  loglike3<-lgamma((df+nu)*0.5)-0.5*nu*log(sigma2y)+as.numeric(determinant(Cholesky(B.uu),logarithm = T)$modulus)-0.5*nu*log(pi*df)-
    lgamma(df/2)-0.5*(df+nu)*log(1+(1/(df*sigma2y))*t(ru)%*%B.uu%*%ru)+sum(log(cal_der_tYJ_wrt_y(yu,gama)))
  
  part2<-part2-2*loglike3 # for part 2 for DIC7=e part 2 from DIC5, - 2*log(yu|theta).
  DIC7<-as.numeric(part1+part2)
  
  return(list(DIC5=DIC5,DIC7=DIC7))
  
}


