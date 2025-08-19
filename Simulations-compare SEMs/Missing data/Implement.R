source("Source.R")

# This script evaluates the performance of different SEMs 
# across datasets with varying characteristics in the presence of missing data.



## Extract spatial weight matrix for houses sold in 1998, from Lucas County, OH house price dataset.

# Define regression formula
hform <- formula(log(price) ~ age + I(age^2) + I(age^3) +
                   log(lotsize) + rooms + log(TLA) + beds + syear)
hmm0 <- model.matrix(hform, data = house)
X<-hmm0

# Create full spatial weight matrix (all years)
hlw <- nb2listw(LO_nb)
W<- as(hlw, "CsparseMatrix") 

# ---- Extract 1998 observations ----
# Logical index: TRUE for homes sold in 1998
indeces.1998 <- data.frame(hmm0)$syear1998 == 1


# Subset W to 1998 homes only
W_1998<-W[indeces.1998,indeces.1998]
w<-W_1998

# Convert W to binary adjacency (1 if neighbors, 0 otherwise)
w[w != 0] <- 1


# Row-standardize so rows sum to 1
w <- row_standardize_to_sum_1(w)



# -----------------------------------------------------------
# 1. Set true parameter values
# -----------------------------------------------------------

# Dataset 1: right-skewed with heavy tails
#   - Yeo–Johnson transformation parameter gamma = 0.5
#   - Student-t errors with df = 4
# Dataset 2: symmetric with light tails
#   - gamma = 1
#   - Student-t errors with df = 30

p=5 # Number of covariates

# Regression coefficients: draw (p+1) values from -3,…,+3 excluding 0
beta<-sample((-3:3)[-3:3 != 0],(p+1),replace = T)
sigma2=0.5 # error variance
df=4 # degrees of freedom (Student-t errors)
rho    <- 0.8      # spatial autocorrelation
gama   <- 0.5        # Yeo–Johnson transformation parameter (gamma)

para=c(beta,c(rho,sigma2,df))

psi=c(1,rep(-1,0,1),-0.1) # True parameter values of missing data model

# -----------------------------------------------------------
# 2. Simulate data from YJ-SEM-t
# -----------------------------------------------------------


# First simulate SEM-t data
sim=simulateSEMterrors(para=para,p=p,weightmat=w)

# Independent variables (X) — includes a leading column of 1's
x=sim$Independent
x=as.matrix(x,nrow=nrow(w))

# Dependent variable before transformation
z=sim$Dependent
z=as.matrix(z)


# Apply Yeo–Johnson inverse transformation 
# to obtain YJ-SEM-t data.
y=t_inverse(z,gama) 


# -----------------------------------------------------------
# 3. Generate missing values according to MNAR
# -----------------------------------------------------------


splitting=MNAR_splitter_all_xm(x=as.matrix(x),xm=xm,y=as.matrix(y),w=w,psi=psi)
m=splitting$m
x=splitting$X
yo=splitting$yo
yu=splitting$yu
w=splitting$W
xm=splitting$xm
# x==xm
mean(m) # This is the missing value percentage
nu=sum(m)
nu
no<-nrow(w)-nu
no



# draw missing vs observed density
density1 <- density(yu)
density2 <- density(yo)
density3<-density(y)
plot(density1, col = "red", lwd = 2,xlim=c(min(y)-1,max(y)+1), main = "yu vs yo", xlab = "Value", ylab = "Density")
lines(density3, col = "blue", lwd = 2,lty=1)
legend("topright", legend = c("yu", "y"), col = c("red", "blue"), lwd = 2)

density3<-density(y)
plot(density2, col = "green", lwd = 2,xlim=c(min(y)-1,max(y)+1), main = "yu vs yo", xlab = "Value", ylab = "Density")
lines(density3, col = "blue", lwd = 2,lty = "dotted")
legend("topright", legend = c("yo", "y"), col = c("green", "blue"), lwd = 2)


plot(density2, col = "green", lwd = 2,xlim=c(min(y)-1,205), main = "yu vs yo", xlab = "Value", ylab = "Density")
lines(density3, col = "blue", lwd = 2,lty = "dotted")
legend("topright", legend = c("yo", "y"), col = c("green", "blue"), lwd = 2)


# -----------------------------------------------------------
# 4. Set starting values
# -----------------------------------------------------------


##4.1 For parameters

woo<-w[1:no,1:no]
xo_<-x[1:no,]
obs.SEM.Gau<-errorsarlm(yo~xo_,listw = mat2listw(woo),zero.policy = T,method = "LU")


rho.dash.start=log(((1+obs.SEM.Gau$lambda)/(1-obs.SEM.Gau$lambda)))
beta.start<-obs.SEM.Gau$coefficients
sigma2.dash.start<-log(obs.SEM.Gau$s2)
df.start<-4 # start df is from fittet df for residuals of SEM
df.dash.start<-log(df.start-3)+0.00001
gama.start<-1
gama.dash.start<-log((gama.start/(2-gama.start))+0.000001) # Add 0.000001, to overcome numerical issues.


## 4.2 For missing values
start.theta=c(beta.start,obs.SEM.Gau$s2,obs.SEM.Gau$lambda,
              df.start,gama.start,psi) 
nu
bsize.set.start<-floor(nu*0.1)
N1.set.start<-100
start.yu.set.start<-yu

start.yu<-generate_posterior_yu_YJSEMt(theta=start.theta,yo,start.yu=start.yu.set.start,x,xm,w,m,bsize.set.start,N1.set.start)
plot(yu,start.yu)
abline(a=0,b=1,col="red") # Plot stsating yu vs true yu




# -----------------------------------------------------------
# 5. Fit different SEM models
# -----------------------------------------------------------



N  <- 300  # Number of HVB iterations
N1 <- 5    # Number of MCMC iterations within each HVB iteration
bsize<-floor(nu*0.1) # Block size (10% of nu, rounded down)
bsize.upall<-bsize

## 5.1 SEM with Gaussian errors


start.theta<-c(beta.start,sigma2.dash.start,rho.dash.start,psi) # use true psi



tic()
fit.SEM.Gau.Ball<-SEM_Gau_MNAR_Buall(x=x,xm=xm,m=m,w=w,yo=yo,startyu=start.yu,p=4,
                                     start_theta=start.theta,N=N,N1 = N1,bsize=bsize.upall)
toc()



## 5.2 YJ-SEM with Gaussian errors


start.theta<-c(beta.start,sigma2.dash.start,rho.dash.start,
               gama.dash.start,psi) 



tic()
fit.YJ.SEM.Gau.Ball<-YJ_SEM_Gau_MNAR_Buall_new2(x,xm,yo,m,N,w,
                                                start.theta,
                                                start.yu=start.yu,p=4,bsize,N1)
toc()



## 5.3 SEM with Student-t errors
#non info
start.theta<-c(beta.start,sigma2.dash.start,rho.dash.start,
               df.dash.start,psi) 


tic()
fit.SEM.t.Ball.n.info<-SEM_t_MNAR_Buall_noninfo(x,xm,m,w,yo,yu,startyu=start.yu,p=4, #
                                                start_theta=start.theta,N,bsize,N1)
toc()



## 5.4 YJ-SEM with Student-t errors
# Non info

start.theta<-c(beta.start,sigma2.dash.start,rho.dash.start,
               df.dash.start,gama.dash.start,psi) 



tic()
fit.YJ.SEM.t.Buall.n.info<-YJ_SEM_t_MNAR_Buall_noninfo(x,xm,m,w,yo,yu.start=start.yu,p=4, #
                                                       start.theta,N,N1,bsize=bsize)
toc()







# -----------------------------------------------------------
# 6. Analyse the results
# -----------------------------------------------------------



## 6.1 Posterior sampling of model parameters using the variational distribution



nparas=length(fit.SEM.Gau.Ball$mu)
sample.size=10000
SEM.Gau.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_SEM_Gau(mu=fit.SEM.Gau.Ball$mu,
                                        B=fit.SEM.Gau.Ball$B,
                                        d=fit.SEM.Gau.Ball$d,
                                        x=x)
  SEM.Gau.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  SEM.Gau.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

SEM.Gau.theta.sample<-SEM.Gau.theta.sample
colnames(SEM.Gau.theta.sample)<-c(rep("bi",ncol(x)+1),"s2","rho",
                                  "psi1","psi2","psi3")
head(SEM.Gau.theta.sample)
apply(SEM.Gau.theta.sample, 2, mean)

#YJ-SEM-Gau

nparas=length(fit.YJ.SEM.Gau.Ball$mu)
# sample.size=10000
YJ.SEM.Gau.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_YJ_SEM_Gau(mu=fit.YJ.SEM.Gau.Ball$mu,
                                           B=fit.YJ.SEM.Gau.Ball$B,
                                           d=fit.YJ.SEM.Gau.Ball$d,
                                           x=x)
  YJ.SEM.Gau.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  YJ.SEM.Gau.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

YJ.SEM.Gau.theta.sample<-YJ.SEM.Gau.theta.sample
colnames(YJ.SEM.Gau.theta.sample)<-c(rep("beta",ncol(x)+1),"s2","rho","gamma",
                                     "psi1","psi2","psi3")
head(YJ.SEM.Gau.theta.sample)
apply(YJ.SEM.Gau.theta.sample, 2, mean)



#SEM-t

nparas=length(fit.SEM.t.Ball.n.info$mu)
# sample.size=10000
SEM.t.n.info.lambda.dash.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_SEM_t(mu=fit.SEM.t.Ball.n.info$mu,
                                      B=fit.SEM.t.Ball.n.info$B,
                                      d=fit.SEM.t.Ball.n.info$d,
                                      x=x)
  SEM.t.n.info.lambda.dash.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  SEM.t.n.info.lambda.dash.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

n<-ncol(w)
SEM.t.n.info.theta.sample<-SEM.t.n.info.lambda.dash.theta.sample[,-(1:n)]
colnames(SEM.t.n.info.theta.sample)<-c(rep("bi",ncol(x)+1),"s2","rho","df",
                                       "psi1","psi2","psi3")

head(SEM.t.n.info.theta.sample)
round(apply(SEM.t.n.info.theta.sample, 2, mean),4)

#YJ-SEM-t


nparas=length(fit.YJ.SEM.t.Buall.n.info$mu)
# sample.size=10000
YJ.SEM.t.n.info.lambda.dash.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_YJ_SEM_t(mu=fit.YJ.SEM.t.Buall.n.info$mu,
                                         B=fit.YJ.SEM.t.Buall.n.info$B,
                                         d=fit.YJ.SEM.t.Buall.n.info$d,
                                         x=x)
  YJ.SEM.t.n.info.lambda.dash.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  YJ.SEM.t.n.info.lambda.dash.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

YJ.SEM.t.n.info.theta.sample<-YJ.SEM.t.n.info.lambda.dash.theta.sample[,-(1:n)]
colnames(YJ.SEM.t.n.info.theta.sample)<-c(rep("betai",p),"sigma2","rho","nu","gama",rep("psi",3))
head(YJ.SEM.t.n.info.theta.sample)
round(apply(YJ.SEM.t.n.info.theta.sample, 2, mean),4)


## 6.2 Posterior densities of model parameters 


# B0


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,1]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,1]
SEM.t.theta.i.sample<-SEM.t.n.info.theta.sample[,1]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,1]

data <- data.frame(
  value = c(SEM.Gau.theta.i.sample, YJ.SEM.Gau.theta.i.sample,
            SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c("Group 1", "Group 2", "Group 3", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c("Group 1" = "blue", "Group 2" = "red", "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c("SEM-Gau", "YJ-SEM-Gau", "SEM-t", "YJ-SEM-t")

# Create the plot
p.B0<-ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  beta[1], linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(beta[0]), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "bottom",  # Place legend at the bottom
    # legend.position = "none",
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
    
  )


p.B0

# B1


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,2]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,2]
SEM.t.theta.i.sample<-SEM.t.n.info.theta.sample[,2]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,2]

data <- data.frame(
  value = c(SEM.Gau.theta.i.sample, YJ.SEM.Gau.theta.i.sample,
            SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c("Group 1", "Group 2", "Group 3", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c("Group 1" = "blue", "Group 2" = "red", "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c("SEM-Gau", "YJ-SEM-Gau", "SEM-t", "YJ-SEM-t")

# Create the plot
ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  beta[2], linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(beta[1]), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Place legend at the bottom
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
    
  )


#s2

ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,ncov+1]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,ncov+1]
SEM.t.theta.i.sample<-SEM.t.n.info.theta.sample[,ncov+1]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,ncov+1]


data <- data.frame(
  value = c(SEM.Gau.theta.i.sample, YJ.SEM.Gau.theta.i.sample,
            SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c("Group 1", "Group 2", "Group 3", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c("Group 1" = "blue", "Group 2" = "red", "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c("SEM-Gau", "YJ-SEM-Gau", "SEM-t", "YJ-SEM-t")

# Create the plot
p.s2<-ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  sigma2, linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(sigma[e]^2), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Place legend at the bottom
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
  )


p.s2

# rho


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,ncov+2]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,ncov+2]
SEM.t.theta.i.sample<-SEM.t.n.info.theta.sample[,ncov+2]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,ncov+2]


data <- data.frame(
  value = c(SEM.Gau.theta.i.sample, YJ.SEM.Gau.theta.i.sample,
            SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c("Group 1", "Group 2", "Group 3", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c("Group 1" = "blue", "Group 2" = "red", "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c("SEM-Gau", "YJ-SEM-Gau", "SEM-t", "YJ-SEM-t")

# Create the plot
p.rho<-ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  rho, linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(rho), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Place legend at the bottom
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
  )


p.rho

#df

SEM.t.theta.i.sample<-SEM.t.n.info.theta.sample[,ncov+3]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,ncov+3]




data <- data.frame(
  value = c(
    SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c( "Group 3", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c( "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c( "SEM-t", "YJ-SEM-t")

# Create the plot
p.df<-ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  df.start, linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", nrow = 1, ncol = 2) +  # 2x2 panel layout
  labs(title = expression(nu), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Place legend at the bottom
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
  )


p.df

#gamma



ncov<-ncol(x)+1

YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,ncov+3]

YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,ncov+4]


data <- data.frame(
  value = c( YJ.SEM.Gau.theta.i.sample,
             YJ.SEM.t.theta.i.sample),
  group = rep(c("Group 2", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c( "Group 2" = "red",  "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c("YJ-SEM-Gau", "YJ-SEM-t")

# Create the plot
p.gamma<-ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  gama, linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(gamma), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Place legend at the bottom
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
  )


p.gamma

#psi0


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,ncov+3]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,ncov+4]
SEM.t.theta.i.sample<-SEM.t.n.info.theta.sample[,ncov+4]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,ncov+5]


data <- data.frame(
  value = c(SEM.Gau.theta.i.sample, YJ.SEM.Gau.theta.i.sample,
            SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c("Group 1", "Group 2", "Group 3", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c("Group 1" = "blue", "Group 2" = "red", "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c("SEM-Gau", "YJ-SEM-Gau", "SEM-t", "YJ-SEM-t")

# Create the plot
p.si0<-ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  psi[1], linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(psi[0]), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Place legend at the bottom
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
  )
p.si0


#psiy


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,ncov+5]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,ncov+6]
SEM.t.theta.i.sample<-SEM.t.n.info.theta.sample[,ncov+6]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.n.info.theta.sample[,ncov+7]


data <- data.frame(
  value = c(SEM.Gau.theta.i.sample, YJ.SEM.Gau.theta.i.sample,
            SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c("Group 1", "Group 2", "Group 3", "Group 4"), each = sample.size)
)

# Create the plot
custom_colors <- c("Group 1" = "blue", "Group 2" = "red", "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c("SEM-Gau", "YJ-SEM-Gau", "SEM-t", "YJ-SEM-t")

# Create the plot
p.siy<-ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  psi[3], linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(psi[y]), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "none",  # Place legend at the bottom
    panel.grid = element_blank(),  # Remove grid lines
    strip.text = element_blank(),  # Remove individual plot titles
    
    plot.title = element_text(size = 25,hjust = 0.5),  # Center the plot title
    legend.title = element_blank(),  # Adjust font size of legend title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "lines")
  )
p.siy


## 6.3 Credible intervals/ means


credible.intervals.SEMGau <- apply(SEM.Gau.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.SEMGau),4)
round(apply(SEM.Gau.theta.sample,2,mean),4)

credible.intervals.YJSEMGau <- apply(YJ.SEM.Gau.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.YJSEMGau),4)
round(apply(YJ.SEM.Gau.theta.sample,2,mean),4)


credible.intervals.SEMt <- apply(SEM.t.n.info.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.SEMt),4)
round(apply(SEM.t.n.info.theta.sample,2,mean),4)


credible.intervals.YJSEMt <- apply(YJ.SEM.t.n.info.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.YJSEMt),4)
round(apply(YJ.SEM.t.n.info.theta.sample,2,mean),4)

## 6.4 Sampling from the posterior of y_u using parallel processing (multiple cores)
# -------------------------------------------------------------------
# Missing values (y_u) are sampled from p(y_u | y_o, m).
# Since this distribution is not available in closed form, 
# MCMC algorithms are used.
# -------------------------------------------------------------------

# Load parallelization libraries
library(parallel)
library(doParallel)

# Set simulation parameters
N1 <- 50        # Number of MCMC iterations per draw
n.thetas <- 100 # Number of posterior parameter samples (theta_i) 
# used to generate y_u samples

# Initialize parallel backend
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Run posterior sampling of missing values (y_u) in parallel
tic()
results <- foreach(i = 1:4, 
                   .packages = c("LaplacesDemon", "extraDistr", "igraph", "MASS", "spdep", 
                                 "tictoc", "Matrix", "mvtnorm", "coda", "ggplot2", 
                                 "mvnfast", "spatialreg", "bestNormalize", "moments"), 
                   .export = c("generate_posterior_yu_SEMGau", "generate_posterior_yu_YJSEMGau", 
                               "generate_posterior_yu_SEMt", "generate_posterior_yu_YJSEMt", 
                               "rmvn")) %dopar% {  
                                 if (i == 1) {
                                   apply(SEM.Gau.theta.sample[1:n.thetas,], 1, generate_posterior_yu_SEMGau, yo, start.yu, x, xm, w, m, bsize, N1)
                                 } else if (i == 2) {
                                   apply(YJ.SEM.Gau.theta.sample[1:n.thetas,], 1, generate_posterior_yu_YJSEMGau, yo, start.yu, x, xm, w, m, bsize, N1)
                                 } else if (i == 3){
                                   apply(YJ.SEM.t.n.info.theta.sample[1:n.thetas,], 1, generate_posterior_yu_YJSEMt, yo, start.yu, x, xm, w, m, bsize, N1)
                                 }
                                 else {
                                   apply(SEM.t.n.info.theta.sample[1:n.thetas,], 1, generate_posterior_yu_SEMt, yo, start.yu, x, xm, w, m, bsize, N1)
                                 }
                               }

toc()
stopCluster(cl)


# Collect posterior samples for missing values
SEM.Gau.post.sample <- results[[1]]
YJ.SEM.Gau.post.sample <- results[[2]]
YJ.SEM.t.n.info.post.sample <- results[[3]]
SEM.t.n.info.post.sample <- results[[4]]


# Posterior estimates of missing values (means vs. true values)
plot(yu,apply(SEM.Gau.post.sample,1,mean)) # SEM-Gau
plot(yu,apply(YJ.SEM.Gau.post.sample,1,mean)) # YJ-SEM-Gau
plot(yu,apply(SEM.t.n.info.post.sample,1,mean)) #SEM-t
plot(yu,apply(YJ.SEM.t.n.info.post.sample,1,mean)) #YJ-SEM-t


# Mean Squared Error (MSE) of missing value predictions
mean((apply(SEM.Gau.post.sample,1,mean)-yu)^2) # SEM-Gau
mean((apply(YJ.SEM.Gau.post.sample,1,mean)-yu)^2) # YJ-SEM-Gau
mean((apply(SEM.t.n.info.post.sample,1,mean)-yu)^2) #SEM-t
mean((apply(YJ.SEM.t.n.info.post.sample,1,mean)-yu)^2) #YJ-SEM-t


## 6.5 Model fit evaluation using DIC (DIC-5 and DIC-7)




tic()
DIC.5.7.SEMGau<-cal_DIC_5_7_MNAR_SEMGau(SEM.Gau.theta.sample=SEM.Gau.theta.sample[1:n.thetas,],SEM.Gau.post.sample,
                                        yo,x,xm,w,m)
toc()


tic()
DIC.5.7.YJSEMGau<-cal_DIC_5_7_MNAR_YJSEMGau(YJ.SEM.Gau.theta.sample=YJ.SEM.Gau.theta.sample[1:n.thetas,],
                                            YJ.SEM.Gau.post.sample,yo,x,xm,w,m)
toc()


tic()
DIC.5.7.noninfo.SEMt<-cal_DIC_5_7_MNAR_SEMt(SEM.t.theta.sample=SEM.t.n.info.theta.sample[1:n.thetas,],SEM.t.n.info.post.sample,yo,x,xm,w,m)
toc()

tic()
DIC.5.7.YJ.SEM.t.noninfo.MNAR<-cal_DIC_5_7_MNAR_YJSEMt(YJ.SEM.t.theta.sample=YJ.SEM.t.n.info.theta.sample[1:n.thetas,],YJ.SEM.t.n.info.post.sample,yo,x,xm,w,m)
toc()

DIC.5.7.SEMGau$DIC5
DIC.5.7.YJSEMGau$DIC5
DIC.5.7.noninfo.SEMt$DIC5
DIC.5.7.YJ.SEM.t.noninfo.MNAR$DIC5








