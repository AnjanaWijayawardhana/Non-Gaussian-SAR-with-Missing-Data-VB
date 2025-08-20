# This script compares parameter estimates of the YJ-SEM-Gau using VB and HMC methods-full data case

source("Source.R")



###################### Simulate from YJ-SEM-Gau

sqrtn=25 #(i.e.n=sqrtn^2=625)

g<-graph.lattice(c(sqrtn,sqrtn)) #5*4 pixels.

g<-connect.neighborhood(g,1) #connect all vertices by one edge.
#plot(g)
wmat<-as_adjacency_matrix(g) # wmat is a matrix.
w<-mat2listw(wmat,style="W") # w is a listw object.
w=listw2mat(w) # w is a matrix object, but a standardized weight matrix
w=Matrix(w) #convert w into a sparse matrix
# matrix(w,ncol=ncol(w))



p=5

beta<-sample((-3:3)[-3:3 != 0],(p+1),replace = T)
sigma2=1
rho<-0.8

gama<-1.25

para=c(beta,c(rho,sqrt(sigma2)))

sim=simulateSEM(para=para,p=p,weightmat=w)

x=sim$Independent 
x=as.matrix(x,nrow=nrow(w))
z=sim$Dependent
z=as.matrix(z)
y=t_inverse(z,gama)
# lm(y~x)
plot(density(y))
# Fit MLE on t errors.

stan.data.list<-make.SATN.data(x,w,y) # Data to be included in the SATN function to implement HMC algorithm


###################### Fit YJ-SEM-Gau using HMC

filepath<-"YJ_SEM_Gau_full.stan"
N<-100 # number of HMC samples
tic()
fit.YJSEMGau.HMC<-stan(file = filepath,
                     data = stan.data.list,chains = 1,iter = N) # default warm upsize=iter/2.
toc()


p<-ncol(x)
rstan::traceplot(fit.YJSEMGau.HMC, main="Traceplot of MCMC Chains")
mcmc.YJSEMGau.chains <- rstan::extract(fit.YJSEMGau.HMC) # Extract posterior sample from the HMC output.
HMC.YJSEMGau.theta.chains<-cbind(mcmc.YJSEMGau.chains$beta,exp(mcmc.YJSEMGau.chains$omegadash),
                               (exp(mcmc.YJSEMGau.chains$rhodash) -1)/(exp(mcmc.YJSEMGau.chains$rhodash) +1),
                              2*exp(mcmc.YJSEMGau.chains$gamadash)/(exp(mcmc.YJSEMGau.chains$gamadash)+1))
colnames(HMC.YJSEMGau.theta.chains)<-c(rep("betai",p+1),"sigma2","rho","gamma")


round(apply(HMC.YJSEMGau.theta.chains, 2, mean),4) # Posterior mean.
round(apply(HMC.YJSEMGau.theta.chains, 2, sd),4) # Posterior standard deviations.





###################### Fit YJ-SEM-Gau using VB


error_model <- errorsarlm(y ~ x, listw=mat2listw(w, style="W"),method = "LU", zero.policy=TRUE) # finding starting values


start.beta<-error_model$coefficients
start.sigma2.dash<-log((error_model$s2))
start.rho.dash=log(((error_model$lambda)/(1-error_model$lambda)))
start.gama<-1
start.gama.dash<-log((start.gama/(2-start.gama))+0.000001) # Add 0.000001, to overcome numerical issues.


start.theta=c(start.beta,start.sigma2.dash,
              start.rho.dash,start.gama.dash) #Starting values in the transformed space


### Fit YJ-SEM-Gau

N<-10000 # Number of VB iterations.
tic()
fit.YJ.SEM.Gau<-YJ_SEM_Gau(x=x,y=y,N=N
                           ,w=w,start.theta,p=4)
toc()



n=ncol(w)
p<-ncol(as.matrix(x))+1
YJ.SEM.Gau.estimates<-(fit.YJ.SEM.Gau$mu)
YJ.SEM.Gau.estimates[p+1]<-exp(YJ.SEM.Gau.estimates[p+1])
YJ.SEM.Gau.estimates[p+2]<-(exp(YJ.SEM.Gau.estimates[p+2])-1)/(exp(YJ.SEM.Gau.estimates[p+2])+1)
YJ.SEM.Gau.estimates[p+3]<-2*exp(YJ.SEM.Gau.estimates[p+3])/(exp(YJ.SEM.Gau.estimates[p+3])+1)
YJ.SEM.Gau.estimates<-as.vector(YJ.SEM.Gau.estimates)

names(YJ.SEM.Gau.estimates)<-c(rep("betai",p),"sigma2","rho","gama") 
round(YJ.SEM.Gau.estimates,4) # Estimated variational means

plot(ts(fit.YJ.SEM.Gau$all_paras[,-(1:(p))]),
     main="Trojectories of sigma2,rho,nu,gamma YJ-SEM-Gau") # Trajectories of variational means of some of the parameters
                                                             #over VB iterations.


### sample from the posterior of model parameters.

nparas=length(fit.YJ.SEM.Gau$mu)
sample.size=length(HMC.YJSEMGau.theta.chains[,1])

nparas=length(fit.YJ.SEM.Gau$mu)
# sample.size=10000
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
colnames(YJ.SEM.Gau.theta.sample)<-names(YJ.SEM.Gau.estimates)
head(YJ.SEM.Gau.theta.sample) # Posterior sample generated from the estimated variational distribution 


### Comparison of posterior densities of model parameters estimated via VB and HMC.


#beta_0


YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,1]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,1]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,1])))
)

# Create the density plot with ggplot2
ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("dotted","solid")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = beta[1], color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(beta[0]), x =  "", y = "Density") +  # Single heading
  theme_minimal() +
  theme(
    legend.position = c(0.2, 0.8),  # Place legend at the bottom
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
# ggsave("sim_full_HMCvsVB_B0.png", dpi = 150, width = 8, height = 7, units = "in")

# beta_1



YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,2]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,2]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,2])))
)

# Create the density plot with ggplot2
ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("dotted","solid")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = beta[2], color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(beta[1]), x =  "", y = "Density") +  # Single heading
  theme_minimal() +
  theme(
    # legend.position = "bottom",  # Place legend at the bottom
    legend.position = "none",
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

# ggsave("sim_full_HMCvsVB_B1.png", dpi = 150, width = 8, height = 7, units = "in")


# beta_2


YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,3]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,3]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,2])))
)

# Create the density plot with ggplot2
ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("dotted","solid")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = beta[3], color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(beta[2]), x =  "", y = "Density") +  # Single heading
  theme_minimal() +
  theme(
    # legend.position = "bottom",  # Place legend at the bottom
    legend.position = "none",
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

# ggsave("sim_full_HMCvsVB_B2.png", dpi = 150, width = 8, height = 7, units = "in")

# s2

YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,"sigma2"]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,"sigma2"]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,"sigma2"])))
)

# Create the density plot with ggplot2
ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("solid", "dotted")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = sigma2, color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(sigma[e]^2), x =  "", y = "Density") +  # Single heading
  theme_minimal() +
  theme(
    # legend.position = "bottom",  # Place legend at the bottom
    legend.position = "none",
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

# ggsave("sim_full_HMCvsVB_s2.png", dpi = 150, width = 8, height = 7, units = "in")

# rho

YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,"rho"]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,"rho"]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,"rho"])))
)

# Create the density plot with ggplot2
ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("dotted","solid")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = rho, color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(rho), x =  "", y = "Density") +  # Single heading
  theme_minimal() +
  theme(
    # legend.position = "bottom",  # Place legend at the bottom
    legend.position = "none",
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

# ggsave("sim_full_HMCvsVB_rho.png", dpi = 150, width = 8, height = 7, units = "in")



# gama

YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,"gama"]



data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,"gamma"]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,"gamma"])))
)

# Create the density plot with ggplot2
ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("dotted","solid")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = gama, color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(gamma), x =  "", y = "Density") +  # Single heading
  theme_minimal() +
  theme(
    # legend.position = "bottom",  # Place legend at the bottom
    legend.position = "none",
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


# ggsave("sim_full_HMCvsVB_gamma.png", dpi = 150, width = 8, height = 7, units = "in")


### convergence plots




colnames(HMC.YJSEMGau.theta.chains)[1:2] <- c("intercept", "beta1")
head(HMC.YJSEMGau.theta.chains)
HMC.YJSEMGau.theta.chains.to.plot<-HMC.YJSEMGau.theta.chains[,c("intercept", "beta1","sigma2","rho","gamma")]


head(HMC.YJSEMGau.theta.chains.to.plot)  

#For HMC

# Example Greek labels (adjust as needed for your parameters)
greek_labels <- list(expression(beta[0]), expression(beta[1]), expression(sigma^2),expression(rho),expression(gamma))

# Set up plot layout
n <- ncol(HMC.YJSEMGau.theta.chains.to.plot)
layout_rows <- ceiling(sqrt(n))
layout_cols <- ceiling(n / layout_rows)

par(mfrow = c(layout_rows, layout_cols), mar = c(3, 3, 2, 1))  # smaller margins

for (i in 1:n) {
  ts_data <- HMC.YJSEMGau.theta.chains.to.plot[, i]
  plot(ts_data, type = "l", col = "red",
       main = greek_labels[[i]], xlab = "Iteration", ylab = "Value")
}


#For VB


nbetas<-ncol(x)+1
n<-ncol(w)
VB.trojectry=(fit.YJ.SEM.Gau$all_paras)[,c(1,2,nbetas+1,nbetas+2,nbetas+3)]
tail(VB.trojectry)
dim(VB.trojectry)
VB.trojectry[,3]<-exp(VB.trojectry[,3])
VB.trojectry[,4]<-((exp(VB.trojectry[,4]) -1)/(exp(VB.trojectry[,4]) +1))
VB.trojectry[,5]<-2*exp(VB.trojectry[,5])/(exp(VB.trojectry[,5])+1)


VB.trojectry<-as.data.frame(VB.trojectry)
VB.trojectry$Time<-1:10000
vb1.trojectry_melted <- melt(VB.trojectry, id.vars = "Time")

legend_labels <- c("intercept",expression(beta[1]),
                   expression(sigma[e]^2), expression(rho),expression(gamma))

line_types <- c("solid", "dashed", "dotted", "dotdash","longdash")

pp.con.YJSEMGau<- ggplot(data = vb1.trojectry_melted, aes(x = Time, y = value, color = variable, linetype = variable)) +
  geom_line(size = 0.25) +
  scale_color_manual(values = rep("red", length(legend_labels)), name = NULL, labels = legend_labels) +
  scale_linetype_manual(values = line_types, name = NULL, labels = legend_labels) + # Specify line types and remove linetype legend title
  labs(title = "",
       x = "Iteration", y = "Value") +
  theme_minimal() +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.25, 0.6),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+ggtitle("YJ-SEM-Gau")+ylim(-1.7,1.5)

pp.con.YJSEMGau



# ggsave(
#   "con_plots/con_VB_vs_HMC_VB_full.png" #ev means evaluate SEMs
#   ,plot = pp.con.YJSEMGau,width = 1500,height = 1000,units = "px")









