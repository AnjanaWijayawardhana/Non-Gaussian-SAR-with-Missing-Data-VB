
# This script compares parameter estimates of the YJ-SEM-Gau using VB and HMC methods-missing data case

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


##### Generate missing values


psi=c(-1,rep(0.5,1),-0.1)
# psi=c(-2.5,rep(0.5,1),0.1)
 #set psi's (changing these values, we can control the missing value %)

# xm<-as.matrix(x[,1])
n<-ncol(w)
xm<-rlnorm(n)

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
nu      #The number of missing values
no<-nrow(w)-nu
no



stan.data.list<-make.SATN.data(x,xm,w,m,yo) # Data to be included in the SATN function to implement HMC algorithm


###################### Fit YJ-SEM-Gau using HMC

filepath<-"YJ_SEM_Gau_miss.stan"

N<-2500

tic()
fit.YJSEMGau.MNNR.HMC<-stan(file = filepath,
                            data = stan.data.list,chains = 1,iter = N) # default warm upsize=iter/2.
toc()




traceplot(fit.YJSEMGau.MNNR.HMC, main="Traceplot of MCMC Chains")


mcmc.YJSEMGau.chains<-as.matrix(fit.YJSEMGau.MNNR.HMC)
HMC.YJSEMGau.theta.dash.chains<-mcmc.YJSEMGau.chains[,(1:(ncol(mcmc.YJSEMGau.chains)-nu-1))]
HMC.YJSEMGau.theta.chains<-HMC.YJSEMGau.theta.dash.chains
HMC.YJSEMGau.theta.chains[,ncol(x)+2]<-exp(HMC.YJSEMGau.theta.chains[,ncol(x)+2])
HMC.YJSEMGau.theta.chains[,ncol(x)+3]<-(exp(HMC.YJSEMGau.theta.chains[,ncol(x)+3]) -1)/(exp(HMC.YJSEMGau.theta.chains[,ncol(x)+3]) +1)

HMC.YJSEMGau.theta.chains[,ncol(x)+4]<-2*exp(HMC.YJSEMGau.theta.chains[,ncol(x)+4])/(exp(HMC.YJSEMGau.theta.chains[,ncol(x)+4])+1)

HMC.YJSEMGau.yu.chains<-mcmc.YJSEMGau.chains[,-(1:(ncol(mcmc.YJSEMGau.chains)-nu-1))]
HMC.YJSEMGau.yu.chains<-HMC.YJSEMGau.yu.chains[,-(nu+1)]
plot(yu,apply(HMC.YJSEMGau.yu.chains,2,mean))

round(apply(as.matrix(HMC.YJSEMGau.theta.chains)[,1:(length(para)+4)], 
            2,mean),4) # Posterior mean.
round(apply(as.matrix(HMC.YJSEMGau.theta.chains)[,1:(length(para)+4)], 
            2,sd),4) # Posterior standard deviations.




###################### Fit YJ-SEM-Gau using HVB-NoB







########### finding starting values

#Set starting parameters

woo<-w[1:no,1:no]
xo_<-x[1:no,]
obs.SEM.Gau<-errorsarlm(yo~xo_,listw = mat2listw(woo),zero.policy = T,method = "LU")


start.rho.dash=log(((1+obs.SEM.Gau$lambda)/(1-obs.SEM.Gau$lambda)))
start.beta<-obs.SEM.Gau$coefficients
start.sigma2.dash<-log(obs.SEM.Gau$s2)
start.gama<-1 # assume no Asymmetry
start.gama.dash<-log((start.gama/(2-start.gama))+0.000001) # Add 0.000001, to overcome numerical issues.
start.psi<-psi



#Set starting missing values


start.theta=c(start.beta,obs.SEM.Gau$s2,
              obs.SEM.Gau$lambda,start.gama,start.psi) 

bsize.set.start<-50
N1.set.start<-50
start.yu.set.start<-yu

start.yu<-generate_posterior_yu_YJSEMGau(theta=start.theta,yo,start.yu=start.yu.set.start,
                                         x,xm,w,m,bsize.set.start,N1.set.start)
plot(yu,start.yu)
abline(a=0,b=1,col="red")



start.theta.dash<-c(start.beta,start.sigma2.dash,start.rho.dash,
               start.gama.dash,psi) #Starting values of parameters 
                                    #in the transformed space

## Fit YJ-SEM-Gau


N=10000 # Number of HVB iterations.
N1=5   # Number of MCMC iterations.
tic()
fit.YJ.SEM.Gau.NoB.n.info<-YJ_SEM_Gau_MNAR_NoB(x,xm,yo,m,N,w,start.theta.dash,
                                               start.yu=yu,p=4,N1)
toc()




# initial analyse


mean(fit.YJ.SEM.Gau.NoB.n.info$accepts)/5 # Acceptance rate in the MCMC step
n=ncol(w)
p<-ncol(as.matrix(x))+1
YJ.SEM.Gau.NoB.estimates<-(fit.YJ.SEM.Gau.NoB.n.info$mu)
YJ.SEM.Gau.NoB.estimates[p+1]<-exp(YJ.SEM.Gau.NoB.estimates[p+1])
YJ.SEM.Gau.NoB.estimates[p+2]<-(exp(YJ.SEM.Gau.NoB.estimates[p+2])-1)/(exp(YJ.SEM.Gau.NoB.estimates[p+2])+1)
YJ.SEM.Gau.NoB.estimates[p+3]<-2*exp(YJ.SEM.Gau.NoB.estimates[p+3])/(exp(YJ.SEM.Gau.NoB.estimates[p+3])+1)
row.names(YJ.SEM.Gau.NoB.estimates)<-c(rep("betai",p),"sigma2","rho","gama",rep("psi",3))
round(YJ.SEM.Gau.NoB.estimates,4) # Estimated variational means

plot(ts((fit.YJ.SEM.Gau.NoB.n.info$all_paras)[,-(1:(p))]),
     main="Trajectory of parameters YJ-SEM-Gau-NoB") # Trajectories of variational means of some of the parameters
                                                      #over VB iterations.


### sample from the posterior of model parameters.


nparas=length(fit.YJ.SEM.Gau.NoB.n.info$mu)
sample.size=1000
YJ.SEM.Gau.theta.sample<-matrix(rep(0,sample.size*(nparas)),ncol = nparas) # samples in transformed parameter space

for(i in 1:sample.size){ 
  
  
  theta.samplei<-simulate_theta_YJ_SEM_Gau(mu=fit.YJ.SEM.Gau.NoB.n.info$mu,
                                           B=fit.YJ.SEM.Gau.NoB.n.info$B,
                                           d=fit.YJ.SEM.Gau.NoB.n.info$d,
                                           x=x)
  YJ.SEM.Gau.theta.sample[i,]<-theta.samplei$theta.dash #  # samples in transformed parameter space
  YJ.SEM.Gau.theta.sample[i,]<-theta.samplei$theta #  # samples in actual parameter space
  
}

YJ.SEM.Gau.theta.sample<-YJ.SEM.Gau.theta.sample
colnames(YJ.SEM.Gau.theta.sample)<-names(YJ.SEM.Gau.NoB.estimates)
head(YJ.SEM.Gau.theta.sample)


### Comparison of posterior densities of model parameters estimated via VB and HMC.

#beta_0



YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,1]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,1]),
  method = rep(c("HVB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,1])))
)

# Create the density plot with ggplot2
B0<-ggplot(data, aes(x = value, color = method, linetype = method)) +
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



# s2



p<-ncol(x)+1
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,p+1]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,p+1]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,p+1])))
)



# Create the density plot with ggplot2
s2<-ggplot(data, aes(x = value, color = method, linetype = method)) +
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


# rho
p<-ncol(x)+1
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,p+2]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,p+2]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,p+1])))
)


# Create the density plot with ggplot2
rho<-ggplot(data, aes(x = value, color = method, linetype = method)) +
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




# gama


p<-ncol(x)+1
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,p+3]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,p+3]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,p+3])))
)





# Create the density plot with ggplot2
gamma<-ggplot(data, aes(x = value, color = method, linetype = method)) +
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




#Psi_0


p<-ncol(x)+1
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,p+4]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,p+4]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,p+4])))
)





# Create the density plot with ggplot2
psi0<-ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("dotted","solid")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = psi[1], color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(psi[0]), x =  "", y = "Density") +  # Single heading
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




#Psi_y


p<-ncol(x)+1
YJ.SEM.Gau.theta.i.sample<-YJ.SEM.Gau.theta.sample[,p+6]


data <- data.frame(
  value = c(YJ.SEM.Gau.theta.i.sample, HMC.YJSEMGau.theta.chains[,p+6]),
  method = rep(c("VB", "HMC"), times = c(length(YJ.SEM.Gau.theta.i.sample), length(HMC.YJSEMGau.theta.chains[,p+6])))
)


psiy<-ggplot(data, aes(x = value, color = method, linetype = method)) +
  geom_density(alpha = 0,size = 2) +
  scale_color_manual(values = c("red", "red")) +   # Purple color for both lines
  scale_linetype_manual(values = c("dotted","solid")) + # Solid line for VB, dotted for HMC
  geom_vline(xintercept = psi[3], color = "orange", size = 2, linetype = "dashed") +  # Vertical line at beta[1]
  # facet_wrap(~ group, scales = "free", ncol = 2) +  # 2x2 panel layout
  labs(title = expression(psi[y]), x =  "", y = "Density") +  # Single heading
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




p.all<-grid.arrange(B0, s2, rho, gamma, psi0, psiy, nrow = 2, ncol = 3)
p.all






### Sample from the posterior of missing values.

N1<-10
bsize<-250
tic()
YJ.SEM.Gau.post.sample<-apply(YJ.SEM.Gau.theta.sample[1:500,], 1, generate_posterior_yu_YJSEMGau,yo,start.yu=fit.YJ.SEM.Gau.NoB.n.info$yunote,x,
                                   xm,w,m,bsize,N1)
toc()

plot(apply(YJ.SEM.Gau.post.sample,1,mean), # Compare posterior mean of yu using HMC and HVB.
     apply(HMC.YJSEMGau.yu.chains,2,mean))



### Convergence plots

#For HMC

kkk<-as.matrix(HMC.YJSEMGau.theta.chains)[,1:length(start.theta)] ## extract parameters 

colnames(kkk)[1:2] <- c("intercept", "beta1")
head(kkk)
HMC.YJSEMGau.theta.chains.to.plot<-kkk[,c("intercept", "beta1","omegadash","rhodash","gamadash")]


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



#For HVB


nbetas<-ncol(x)+1
n<-ncol(w)
VB.trojectry=(fit.YJ.SEM.Gau.NoB.n.info$all_paras)[,c(1,2,nbetas+1,nbetas+2,nbetas+3)]
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

# Define line types
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
  )+ggtitle("YJ-SEM-Gau")+ylim(-3.2,1.5)

pp.con.YJSEMGau













