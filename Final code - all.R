#########################################################################
#### Functions used for Inference from Adult to Pediatric data FINAL ####
####                         Jianan Hui                              ####
#########################################################################


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
  
  res=data.frame("mixture19"=reject_mixture19,"minimax"=reject_minimax,"regular"=reject_regular,"freq"=reject_freq,"mixture55"=reject_mixture55)
  return(colMeans(res))
  
}

intensity=100
list_sample_size=list(c(1000,1000),c(1000,800),c(1000,600),c(1000,400),c(1000,200),c(1000,100),c(1000,50),c(1000,25))

###Under null hypothesis
#Scenario 0: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=-0.05,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_0=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario I: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_1=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])


#Scenario II: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)

SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_2=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])


#Scenario III: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)

### investigating on small sample sizes of pediatric data
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_3=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario I: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_1=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])


#Scenario IV: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_4=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario V: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_5=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])


#Scenario VI: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_6=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario VIII: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_8=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario IX: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_9=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

result_H0.1 <- result_H0_0
result_H0.2 <- result_H0_1
result_H0.3 <- result_H0_2
result_H0.4 <- result_H0_3
result_H0.5 <- result_H1_1

result_H1.1 <- result_H1_4
result_H1.2 <- result_H1_5
result_H1.3 <- result_H1_6
result_H1.4 <- result_H1_8
result_H1.5 <- result_H1_9

####Simulation to address comments on Vary N and Vary R

intensity=5000
list_sample_size=list(c(500,500),c(500,400),c(500,300),c(500,200),c(500,100),c(500,50),c(500,25))

###Under null hypothesis
#Scenario I: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H00_1=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario II: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H00_2=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])


#Scenario III: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H00_3=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

###Under alternative hypothesis
#Scenario I: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H11_1=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario II: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=1,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H11_2=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

#Scenario III: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=7^2,n_a=x[1],n_p=x[2],mu_p=1.5,var_p=7^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
VaryN_H11_3=data.table(SampleSize_a_p=list_sample_size,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4],mixture55=SI_res[,5])

# Varying r, which is the proportion of effect size p over effect size a.
# r=0, 0.25, 0.5, 0.65, 0.8, 1
# 
# Under null hypothesis:
#   
# H00-1. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=5^2)
# H00-2. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=10^2)
# H00-3. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=15^2)
# 
# Under alternative hypothesis:
#   
# H11-1. Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=10^2)
# H11-2. Under H1: (mu_a=1, var_a=10^2, mu_p=2, var_p=10^2)
# H11-3. Under H1: (mu_a=1, var_a=10^2, mu_p=1.5, var_p=10^2)


r_gamma <- c(0,0.25,0.5,0.65,0.8,1)
intensity <- 5000

n_p <- c(600,400,200,100)
#Under null hypothesis

n_pp <- n_p[1]
# H00-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_1<- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

#Under alternative hypothesis

# H11-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0.8,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_1 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

VaryR_600_H00_all <- data.frame("r"=VaryR_H00_1$r,"H00-1"=VaryR_H00_1$minimax,"H00-2"=VaryR_H00_2$minimax,"H00-3"=VaryR_H00_3$minimax)
VaryR_600_H11_all <- data.frame("r"=VaryR_H11_1$r,"H00-1"=VaryR_H11_1$minimax,"H00-2"=VaryR_H11_2$minimax,"H00-3"=VaryR_H11_3$minimax)

####
n_pp <- n_p[2]
# H00-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_1<- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

#Under alternative hypothesis

# H11-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0.8,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_1 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

VaryR_400_H00_all <- data.frame("r"=VaryR_H00_1$r,"H00-1"=VaryR_H00_1$minimax,"H00-2"=VaryR_H00_2$minimax,"H00-3"=VaryR_H00_3$minimax)
VaryR_400_H11_all <- data.frame("r"=VaryR_H11_1$r,"H00-1"=VaryR_H11_1$minimax,"H00-2"=VaryR_H11_2$minimax,"H00-3"=VaryR_H11_3$minimax)

####

####
n_pp <- n_p[3]
# H00-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_1<- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

#Under alternative hypothesis

# H11-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0.8,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_1 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

VaryR_200_H00_all <- data.frame("r"=VaryR_H00_1$r,"H00-1"=VaryR_H00_1$minimax,"H00-2"=VaryR_H00_2$minimax,"H00-3"=VaryR_H00_3$minimax)
VaryR_200_H11_all <- data.frame("r"=VaryR_H11_1$r,"H00-1"=VaryR_H11_1$minimax,"H00-2"=VaryR_H11_2$minimax,"H00-3"=VaryR_H11_3$minimax)

####

####
n_pp <- n_p[4]
# H00-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_1<- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H00-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H00_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

#Under alternative hypothesis

# H11-1
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=0.8,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_1 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-2
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_2 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

# H11-3
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=1000,n_p=n_pp,mu_p=1.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity,gamma = x)
SI <- sapply(r_gamma,res_I)
SI_res <- t(SI)
VaryR_H11_3 <- data.table(r=r_gamma,mixture19=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,5],mixture55=SI_res[,1])

VaryR_100_H00_all <- data.frame("r"=VaryR_H00_1$r,"H00-1"=VaryR_H00_1$minimax,"H00-2"=VaryR_H00_2$minimax,"H00-3"=VaryR_H00_3$minimax)
VaryR_100_H11_all <- data.frame("r"=VaryR_H11_1$r,"H00-1"=VaryR_H11_1$minimax,"H00-2"=VaryR_H11_2$minimax,"H00-3"=VaryR_H11_3$minimax)

####
save.image("FinalData.RData")

#######################
####Rendering plots####
#######################


####Producing plots for the main simulation section
#producing plots
library(ggplot2)

result_H0.1$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0.1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p1 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = -0.05 and ",sigma[p]," = 5")))+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p1

result_H0.2$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0.2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p2 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p2

result_H0.3$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0.3
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p3

result_H0.4$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0.4
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p4 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0 and ",sigma[p]," = 15")))+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p4

result_H0.5$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0.5
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p5 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0.5 and ",sigma[p]," = 10")))+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p5


#alternative
result_H1.1$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1.1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p6 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 0.8 and ",sigma[p]," = 5")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p6

result_H1.2$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1.2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p7 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 0.8 and ",sigma[p]," = 10")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p7

result_H1.3$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1.3
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p8 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 0.8 and ",sigma[p]," = 15")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p8

result_H1.4$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1.4
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p9 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1 and ",sigma[p]," = 10")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p9

result_H1.5$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1.5
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:6]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p10 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1.5 and ",sigma[p]," = 10")))+
  # theme(
  #   legend.position = c(.95, .05),
  #   legend.justification = c("right", "bottom"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p10


library(ggpubr)
ggarrange(p2, p3, p4, p5, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
ggarrange(p6, p8, p9, p10, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


####Producing plots for the discussion and addressing comments for Vary N and Vary R
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
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
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4,5),labels=c("Robust mixture prior (w=0.9)","Profile Bayesian","Regular Bayesian","Frequentist","Robust mixture prior (w=0.5)"))
p6


###Rendering plots for varying r


###Pediatric sample size = 600
data_wide=VaryR_600_H00_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p7 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Type I Error", title = "Adult sample size = 1000 and Pediatric sample size = 600")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Standard deviation",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 10","Sigma = 15"))
p7


data_wide=VaryR_600_H11_all
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
  scale_colour_discrete(name="Mean",breaks=c(1,2,3),labels=c("Mean = 0.8","Mean = 1","Mean = 1.5"))
p8

###Pediatric sample size = 400
data_wide=VaryR_400_H00_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p9 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Type I Error", title = "Adult sample size = 1000 and Pediatric sample size = 400")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Standard deviation",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 10","Sigma = 15"))
p9


data_wide=VaryR_400_H11_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p10 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Power", title = "Adult sample size = 1000 and Pediatric sample size = 400")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Standard deviation",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 10","Sigma = 15"))
p10

###Pediatric sample size = 200
data_wide=VaryR_200_H00_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p11 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Type I Error", title = "Adult sample size = 1000 and Pediatric sample size = 200")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Standard deviation",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 10","Sigma = 15"))
p11


data_wide=VaryR_200_H11_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p12 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Power", title = "Adult sample size = 1000 and Pediatric sample size = 200")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Mean",breaks=c(1,2,3),labels=c("Mean = 0.8","Mean = 1","Mean = 1.5"))
p12


###Pediatric sample size = 100
data_wide=VaryR_100_H00_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p13 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Type I Error", title = "Adult sample size = 1000 and Pediatric sample size = 100")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Standard deviation",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 10","Sigma = 15"))
p13


data_wide=VaryR_100_H11_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p14 = ggplot(data=data,aes(x=r,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Power", title = "Adult sample size = 1000 and Pediatric sample size = 100")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_colour_discrete(name="Mean",breaks=c(1,2,3),labels=c("Mean = 0.8","Mean = 1","Mean = 1.5"))
p14


library(ggpubr)

#VaryN
ggarrange(p1, p2, p3, ncol=3, common.legend = TRUE, legend="bottom")

ggarrange(p4, p5, p6, ncol=3, common.legend = TRUE, legend="bottom")

#VaryR
ggarrange(p13, p11, p9, p7, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

ggarrange(p14, p12, p10, p8, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


ggarrange(p11, p9, ncol=2, common.legend = TRUE, legend="bottom")

ggarrange(p12, p10, ncol=2, common.legend = TRUE, legend="bottom")


#####line type VaryR####
###Pediatric sample size = 400
data_wide=VaryR_400_H00_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p9 = ggplot(data=data,aes(x=r,y=Value,group=time))+
  geom_line(aes(linetype=time))+
  geom_point()+
  geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Type I Error", title = "Adult sample size = 1000 and Pediatric sample size = 400")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_linetype_discrete(name="Standard deviation",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 10","Sigma = 15"))
p9


data_wide=VaryR_400_H11_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p10 = ggplot(data=data,aes(x=r,y=Value,group=time))+
  geom_line(aes(linetype=time))+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Power", title = "Adult sample size = 1000 and Pediatric sample size = 400")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_linetype_discrete(name="Mean",breaks=c(1,2,3),labels=c("Mean = 0.8","Mean = 1","Mean = 1.5"))
p10

###Pediatric sample size = 200
data_wide=VaryR_200_H00_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p11 = ggplot(data=data,aes(x=r,y=Value,group=time))+
  geom_line(aes(linetype=time))+
  geom_point()+
  geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Type I Error", title = "Adult sample size = 1000 and Pediatric sample size = 200")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_linetype_discrete(name="Standard deviation",breaks=c(1,2,3),labels=c("Sigma = 5","Sigma = 10","Sigma = 15"))
p11


data_wide=VaryR_200_H11_all
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="r")
data$time=factor(data$time)
p12 = ggplot(data=data,aes(x=r,y=Value,group=time))+
  geom_line(aes(linetype=time))+
  geom_point()+
  #geom_hline(yintercept=0.025,linetype="dashed",color="darkgrey")+
  labs(x = "r", y = "Power", title = "Adult sample size = 1000 and Pediatric sample size = 200")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )+
  # theme(legend.position="right","top")+
  scale_linetype_discrete(name="Mean",breaks=c(1,2,3),labels=c("Mean = 0.8","Mean = 1","Mean = 1.5"))
p12


ggarrange(p11, p9, ncol=2, common.legend = TRUE, legend="bottom")

ggarrange(p12, p10, ncol=2, common.legend = TRUE, legend="bottom")