source("Source.R")


# This script compares the performance of different SEMs 
# across datasets with varying characteristics, using full data.


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
# 3. Set starting values
# -----------------------------------------------------------

lw <- mat2listw(w, style="W", zero.policy=TRUE)
error_model <- errorsarlm(y ~ x, listw = lw, method = "LU", zero.policy=TRUE)

beta.start<-error_model$coefficients
sigma2.start<-error_model$s2
rho.dash.start<-log(((error_model$lambda)/(1-error_model$lambda)))
df.start<-4
gama.start<-1
gama.dash.start<-log((gama.start/(2-gama.start))+0.000001) # Add 0.000001, to overcome numerical issues.

# -----------------------------------------------------------
# 4. Fit different SEM models
# -----------------------------------------------------------


N=15000 # The number of VB iterations.


## 4.1 SEM with Gaussian errors

start.theta=c(beta.start,log(sigma2.start),rho.dash.start)
tic()
fit.SEM.Gau<-SEM_Gau(x=x,y=y,N=N,w=w,start.theta,p=4)
toc()


## 4.2 YJ-SEM with Gaussian errors
start.theta=c(beta.start,log(sigma2.start),rho.dash.start,gama.dash.start)
tic()
fit.YJ.SEM.Gau<-YJ_SEM_Gau(x=x,y=y,N=N
                           ,w=w,start.theta,p=4)
toc()


## 4.3 SEM with Student-t errors
start.theta=c(beta.start,log(sigma2.start),rho.dash.start,log(df.start-3)) 


tic()
fit.SEM.terrors.noninfo<-SEM_terrors_noninfo(x,y,w,p=4,start.theta=start.theta,N)
toc()



## 4.4 YJ-SEM with Student-t errors
start.theta=c(beta.start,log(sigma2.start),rho.dash.start,
              log(df.start-3),gama.dash.start)


tic()
fit.YJ.SEM.terrors.noninfo=YJ_SEM_terrors_noninfo(x,y,w,p=4,start.theta=start.theta,N)
toc()



# -----------------------------------------------------------
# 5. Analyse the results
# -----------------------------------------------------------



## 5.1 Posterior sampling of model parameters using the variational distribution


#SEM-Gau

nparas=length(fit.SEM.Gau$mu)
sample.size=100
SEM.Gau.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_SEM_Gau(mu=fit.SEM.Gau$mu,
                                        B=fit.SEM.Gau$B,
                                        d=fit.SEM.Gau$d,
                                        x=x)
  SEM.Gau.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  SEM.Gau.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

SEM.Gau.theta.sample<-SEM.Gau.theta.sample
p<-ncol(as.matrix(x))+1
colnames(SEM.Gau.theta.sample)<-c(rep("betai",p),"sigma2","rho")
head(SEM.Gau.theta.sample)
apply(SEM.Gau.theta.sample, 2, mean)

#YJ-SEM-Gau

nparas=length(fit.YJ.SEM.Gau$mu)
YJ.SEM.Gau.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_YJ_SEM_Gau(mu=fit.YJ.SEM.Gau$mu,
                                           B=fit.YJ.SEM.Gau$B,
                                           d=fit.YJ.SEM.Gau$d,
                                           x=x)
  YJ.SEM.Gau.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  YJ.SEM.Gau.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

YJ.SEM.Gau.theta.sample<-YJ.SEM.Gau.theta.sample
colnames(YJ.SEM.Gau.theta.sample)<-c(rep("betai",p),"sigma2","rho","gamma")
head(YJ.SEM.Gau.theta.sample)
apply(YJ.SEM.Gau.theta.sample, 2, mean)



#SEM-t
nparas=length(fit.SEM.terrors.noninfo$mu)
# sample.size=10000
SEM.t.lambda.dash.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_SEM_t(mu=fit.SEM.terrors.noninfo$mu,
                                      B=fit.SEM.terrors.noninfo$B,
                                      d=fit.SEM.terrors.noninfo$d,
                                      x=x)
  SEM.t.lambda.dash.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  SEM.t.lambda.dash.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

n<-ncol(w)
SEM.t.theta.sample<-SEM.t.lambda.dash.theta.sample[,-(1:n)]
colnames(SEM.t.theta.sample)<-c(rep("betai",p),"sigma2","rho","nu")
head(SEM.t.theta.sample)
apply(SEM.t.theta.sample, 2, mean)

#YJ-SEM-t

nparas=length(fit.YJ.SEM.terrors.noninfo$mu)
# sample.size=10000
YJ.SEM.t.lambda.dash.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_YJ_SEM_t(mu=fit.YJ.SEM.terrors.noninfo$mu,
                                         B=fit.YJ.SEM.terrors.noninfo$B,
                                         d=fit.YJ.SEM.terrors.noninfo$d,
                                         x=x)
  YJ.SEM.t.lambda.dash.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  YJ.SEM.t.lambda.dash.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

YJ.SEM.t.theta.sample<-YJ.SEM.t.lambda.dash.theta.sample[,-(1:n)]
colnames(YJ.SEM.t.theta.sample)<-c(rep("betai",p),"sigma2","rho","nu","gama")
head(YJ.SEM.t.theta.sample)
apply(YJ.SEM.t.theta.sample, 2, mean)



## 5.2 Posterior densities of model parameters 

# B0


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,1]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,1]
SEM.t.theta.i.sample<-SEM.t.theta.sample[,1]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.theta.sample[,1]

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
  geom_vline(xintercept =  beta[1], linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
  facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(beta[0]), x =  "", y = "Density") +  # Single heading
  scale_color_manual(values = custom_colors, labels = custom_labels) +  # Set custom colors and labels for each group
  theme_minimal() +  # Clean theme
  theme(
    legend.position = "bottom",  # Place legend at the bottom
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


#

# B1


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,2]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,2]
SEM.t.theta.i.sample<-SEM.t.theta.sample[,2]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.theta.sample[,2]

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
SEM.t.theta.i.sample<-SEM.t.theta.sample[,ncov+1]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.theta.sample[,ncov+1]


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



# rho


ncov<-ncol(x)+1
SEM.Gau.theta.i.sample<-SEM.Gau.theta.sample[,ncov+2]
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,ncov+2]
SEM.t.theta.i.sample<-SEM.t.theta.sample[,ncov+2]
YJ.SEM.t.theta.i.sample<-YJ.SEM.t.theta.sample[,ncov+2]


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


#df

SEM.t.theta.i.sample<-SEM.t.theta.sample[,ncov+3]

YJ.SEM.t.theta.i.sample<-YJ.SEM.t.theta.sample[,ncov+3]

data <- data.frame(
  value = c(
    SEM.t.theta.i.sample,YJ.SEM.t.theta.i.sample),
  group = rep(c( "Group 3", "Group 4"), each = length(YJ.SEM.t.theta.i.sample))
)

# Create the plot
custom_colors <- c( "Group 3" = "green", "Group 4" = "purple")

# Create custom legend labels
custom_labels <- c( "SEM-t", "YJ-SEM-t")

# Create the plot
ggplot(data, aes(x = value, color = group)) +
  geom_density(alpha = 0,size=2) +  # Remove fill, keep line color
  geom_vline(xintercept =  df, linetype = "dashed", color = "orange", size = 2) +  # Add vertical line at x = 5
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




#gamma



ncov<-ncol(x)+1

YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,ncov+3]

YJ.SEM.t.theta.i.sample<-YJ.SEM.t.theta.sample[,ncov+4]


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
ggplot(data, aes(x = value, color = group)) +
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



 
## 5.3 Credible intervals



credible.intervals.SEMGau <- apply(SEM.Gau.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.SEMGau),4)
round(apply(SEM.Gau.theta.sample,2,mean),4)


credible.intervals.YJSEMGau <- apply(YJ.SEM.Gau.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.YJSEMGau),4)
round(apply(YJ.SEM.Gau.theta.sample,2,mean),4)

credible.intervals.SEMt <- apply(SEM.t.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.SEMt),4)
round(apply(SEM.t.theta.sample,2,mean),4)

credible.intervals.YJSEMt <- apply(YJ.SEM.t.theta.sample, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
round(t(credible.intervals.YJSEMt),4)
round(apply(YJ.SEM.t.theta.sample,2,mean),4)



## 5.4 Cal model fit

#### DIC
no.of.samples.consider<-sample.size
tic()
DIC.SEM.Gau<-cal_DIC_SEMGau(SEM.Gau.theta.sample[1:no.of.samples.consider,],y,x,w)
toc()
tic()
DIC.SEM.Gau.YJ<-cal_DIC_YJSEMGau(YJ.SEM.Gau.theta.sample[1:no.of.samples.consider,],y,x,w)
toc()
tic()
DIC.SEM.Gau.t<-cal_DIC_SEMt(SEM.t.lambda.dash.theta.sample[1:no.of.samples.consider,],y,x,w)
toc()

tic()
DIC.YJ.SEM.t<-cal_DIC_YJSEMt(YJ.SEM.t.lambda.dash.theta.sample[1:no.of.samples.consider,],y,x,w)
toc()


DIC.SEM.Gau$DIC2
DIC.SEM.Gau.t$DIC2
DIC.SEM.Gau.YJ$DIC2
DIC.YJ.SEM.t$DIC2


DIC.SEM.Gau$DIC1
DIC.SEM.Gau.t$DIC1
DIC.SEM.Gau.YJ$DIC1
DIC.YJ.SEM.t$DIC1


