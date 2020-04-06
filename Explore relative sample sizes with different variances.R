######################################################################
####  New simulation settings to accomondate reviewer's comments  ####
####                  Jianan Hui, Apr 3, 2020                     ####
######################################################################

# Varying relative sample sizes (ratio of the two sample sizes from adult to pediatric). Specifically, add scenarios where adult sample size = 500 and pediatric sample size varies from 500, 400, 300, 200, 100, 50 and 25.
# The following assumptions will be adopted:
#   n_a=500,n_p=500,400,300,200,100,50,25 
# Under null hypothesis:
#   
# H00-1. Under H0: (mu_a=1, var_a=7^2, mu_p=0, var_p=5^2)
# H00-2. Under H0: (mu_a=1, var_a=7^2, mu_p=0, var_p=7^2)
# H00-3. Under H0: (mu_a=1, var_a=7^2, mu_p=0, var_p=10^2)
# 
# Under alternative hypothesis:
#   
# H11-1. Under H1: (mu_a=1, var_a=7^2, mu_p=0.8, var_p=7^2)
# H11-2. Under H1: (mu_a=1, var_a=7^2, mu_p=2, var_p=7^2)
# H11-3. Under H1: (mu_a=1, var_a=7^2, mu_p=1.5, var_p=7^2)

setwd("/Users/jianan/Dropbox/Side projects/Profile Bayesian/Profile-Bayesian/")
###Simulating function###
set.seed(1000)
library(data.table)
library(ggplot2)
library(RBesT)

###Continuous endpoints
Bayes_continuous=function(mu_a,var_a,n_a,n_p,mu_p,var_p,n.samples,alpha,rep,gamma=0.5){
  
  #initialize res vectors
  reject_mixture55=reject_mixture19=reject_minimax=reject_regular=reject_freq=double(rep)
  
  #simulate adult data
  set.seed(25)
  data_a=rnorm(n_a,mu_a,sqrt(var_a))
  #sd_a=sd(data_a)
  sd_a=sqrt(var_a)
  mean_a=mean(data_a)
  mean_a
  for (i in 1:rep){
    
    #simulate pediatric data
    #n_p=ceiling(p*n_a)
    data_p=rnorm(n_p,mu_p,sqrt(var_p))
    mean_p=mean(data_p)
    sd_p=sqrt(var_p)
    #sd_p=sd(data_p)
    
    #simulate parameter of interest theta
    #minimax
    mu_theta_minimax=mean_p
    #regular
    #w1=(var_a/n_a)/((var_a/n_a)+(var_p/n_p))
    w1=(n_a/sd_a^2)/(n_a/sd_a^2+n_p/sd_p^2)
    #w1=1/(n_a/var_a+n_p/var_p)*n_a/var_a
    w2=1-w1
    mu_theta_regular=w2*mean_p+w1*mean_a
    #common sd
    sd_theta=sqrt(1/(n_a/sd_a^2+n_p/sd_p^2))
    
    #simulating
    if (mean_p>gamma*mean_a){
      theta_minimax=rnorm(n.samples,mu_theta_minimax,sd_theta)
      reject_minimax[i]=ifelse(mean(theta_minimax<0)<=alpha,1,0)
    }else{
      reject_minimax[i]=ifelse(mean_p*sqrt(n_p)/sd_p>qnorm(1-alpha),1,0)
    }
    
    theta_regular=rnorm(n.samples,mu_theta_regular,sd_theta)
    
    nm <- mixnorm(adult=c(1, mean_a, sd_a/sqrt(n_a)), sigma=sd_a)
    rnMix55 <- robustify(nm, weight=0.5, mean=0, n=1, sigma=sd_a)
    posterior.sum55 <- postmix(rnMix55, m=mean_p, n=n_p, sigma=sd_p)
    components55 <- sample(1:2,size=n.samples,prob=posterior.sum55[1,],replace=TRUE)
    mus <- posterior.sum55[2,]
    sds <- posterior.sum55[3,]
    theta_mixture55 <- rnorm(n.samples,mean=mus[components55],sd=sds[components55])
    #estimate probability of getting theta estimate that is greater than zero
    reject_mixture55[i]=ifelse(mean(theta_mixture55<0)<=alpha,1,0)
    
    rnMix19 <- robustify(nm, weight=0.1, mean=0, n=1, sigma=sd_a)
    posterior.sum19 <- postmix(rnMix19, m=mean_p, n=n_p, sigma=sd_p)
    components19 <- sample(1:2,size=n.samples,prob=posterior.sum19[1,],replace=TRUE)
    mus <- posterior.sum19[2,]
    sds <- posterior.sum19[3,]
    theta_mixture19 <- rnorm(n.samples,mean=mus[components19],sd=sds[components19])
    #estimate probability of getting theta estimate that is greater than zero
    reject_mixture19[i]=ifelse(mean(theta_mixture19<0)<=alpha,1,0)
    
    reject_regular[i]=ifelse(mean(theta_regular<0)<=alpha,1,0)
    reject_freq[i]=ifelse(mean_p*sqrt(n_p)/sd_p>qnorm(1-alpha),1,0)
    
  }
  
  res=data.frame("mixture55"=reject_mixture55,"mixture19"=reject_mixture19,"minimax"=reject_minimax,"regular"=reject_regular,"freq"=reject_freq)
  return(colMeans(res))
  
}


intensity=5000
list_sample_size=list(c(500,500),c(500,400),c(500,300),c(500,200),c(500,100),c(500,50),c(500,25))

###Under null hypothesis
#Scenario I: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H00_1=data.table(SampleSize_a_p=list_sample_size,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

#Scenario II: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H00_2=data.table(SampleSize_a_p=list_sample_size,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])


#Scenario III: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H00_3=data.table(SampleSize_a_p=list_sample_size,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

###Under alternative hypothesis
#Scenario I: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H11_1=data.table(SampleSize_a_p=list_sample_size,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

#Scenario II: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=1,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H11_2=data.table(SampleSize_a_p=list_sample_size,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

#Scenario III: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=1.5,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H11_3=data.table(SampleSize_a_p=list_sample_size,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])



# Varying r, which is the proportion of effect size p over effect size a.
# r=0, 0.25, 0.5, 0.65, 0.8, 1
# 
# Under null hypothesis:
#   
# H00-1. Under H0: (mu_a=1, var_a=7^2, mu_p=0, var_p=5^2)
# H00-2. Under H0: (mu_a=1, var_a=7^2, mu_p=0, var_p=7^2)
# H00-3. Under H0: (mu_a=1, var_a=7^2, mu_p=0, var_p=10^2)
# 
# Under alternative hypothesis:
#   
# H11-1. Under H1: (mu_a=1, var_a=7^2, mu_p=0.8, var_p=7^2)
# H11-2. Under H1: (mu_a=1, var_a=7^2, mu_p=2, var_p=7^2)
# H11-3. Under H1: (mu_a=1, var_a=7^2, mu_p=1.5, var_p=7^2)

r_gamma = c(0,0.25,0.5,0.65,0.8,1)
intensity=5000

#Under null hypothesis

# H00-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=1000,n_p=600,mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_1=data.table(r=r_gamma,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

# H00-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=1000,n_p=600,mu_p=0,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_2=data.table(r=r_gamma,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

# H00-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=1000,n_p=600,mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_3=data.table(r=r_gamma,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

#Under alternative hypothesis

# H11-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=1000,n_p=600,mu_p=0.8,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_1=data.table(r=r_gamma,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

# H11-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=1000,n_p=600,mu_p=1,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_2=data.table(r=r_gamma,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])

# H11-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=1000,n_p=600,mu_p=1.5,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_3=data.table(r=r_gamma,mixture55=SI_res[,1],mixture91=SI_res[,2],minimax=SI_res[,3],regular=SI_res[,4],frequentist=SI_res[,5])


#Rendering plots
library(ggplot2)
pdf("VaryNVaryR-Images.pdf")
VaryN_H00_1$ss_p=c(500,400,300,200,100,50,25)
data_wide=VaryN_H00_1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p1 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0 and ",sigma[p]," = 5")))+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w = 0.5)","Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist"))
p1


VaryN_H00_2$ss_p=c(500,400,300,200,100,50,25)
data_wide=VaryN_H00_2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p2 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0 and ",sigma[p]," = 7")))+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w = 0.5)","Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist"))
p2

VaryN_H00_3$ss_p=c(500,400,300,200,100,50,25)
data_wide=VaryN_H00_3
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p3 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0 and ",sigma[p]," = 10")))+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w = 0.5)","Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist"))
p3

#alternative
VaryN_H11_1$ss_p=c(500,400,300,200,100,50,25)
data_wide=VaryN_H11_1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p4 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 0.8 and ",sigma[p]," = 7")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w = 0.5)","Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist"))
p4

VaryN_H11_2$ss_p=c(500,400,300,200,100,50,25)
data_wide=VaryN_H11_2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p5 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1 and ",sigma[p]," = 7")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w = 0.5)","Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist"))
p5


VaryN_H11_3$ss_p=c(500,400,300,200,100,50,25)
data_wide=VaryN_H11_3
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p6 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1.5 and ",sigma[p]," = 7")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w = 0.5)","Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist"))
p6


###Rendering plots for r

#Combine plots for R
VaryR_H00_all <- data.frame("r"=VaryR_H00_1$r,"H00-1"=VaryR_H00_1$minimax,"H00-2"=VaryR_H00_2$minimax,"H00-3"=VaryR_H00_3$minimax)

data_wide=VaryR_H00_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p7 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Type I Error", title = "Adult sample size = 1000 and Pediatric sample size = 600")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 7","Sigma = 10"))
p7


VaryR_H11_all <- data.frame("r"=VaryR_H11_1$r,"H00-1"=VaryR_H11_1$minimax,"H00-2"=VaryR_H11_2$minimax,"H00-3"=VaryR_H11_3$minimax)
data_wide=VaryR_H11_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p8 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Power", title = "Adult sample size = 1000 and Pediatric sample size = 600")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3),labels=c("Mean = 0.8","Mean = 1","Mean = 1.5"))
p8

dev.off()