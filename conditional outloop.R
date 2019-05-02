
######################################################################
#### Functions used for Inference from Adult to Pediatric data V2 ####
####                       Jianan Hui                             ####
######################################################################

##Three updates were made: 
# .	Change weight of non-informative prior from 0.5 to 0.1. 
# .	Change sample variance to population variance. 

# Under null hypothesis:
#   
# H0-0. Under H0: (mu_a=1, var_a=10^2, mu_p=-0.05, var_p=5^2)
# H0-1. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=5^2)
# H0-2. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=10^2)
# H0-3. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=15^2)
# H0-4. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=20^2)
# 
# Under alternative hypothesis:
#   
# H1-0. Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=5^2)
# H1-1. Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=10^2)
# H1-2. Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=15^2)
# H1-3. Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=20^2)
# H1-4. Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=5^2)
# H1-5. Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=10^2)
# H1-6. Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=15^2)
# H1-7. Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=20^2)
# H1-8. Under H1: (mu_a=1, var_a=10^2, mu_p=1, var_p=10^2)
# H1-9. Under H1: (mu_a=1, var_a=10^2, mu_p=1.5, var_p=10^2)
setwd("/Users/jianan/Dropbox/Side projects/Profile Bayesian/Profile-Bayesian/")
###Simulating function###
set.seed(1000)
library(data.table)
library(ggplot2)
library(RBesT)
#setwd("//eu.boehringer.com/users/rdg/users4/jhui/Documents/Side projects/Minimax Adult to Pediatric/Minimax simulation/Add adjusted alpha and visualization on rejection region")

###Continuous endpoints
Bayes_continuous=function(mu_a,var_a,n_a,n_p,mu_p,var_p,n.samples,alpha,rep){
  
  #initialize res vectors
  reject_mixture=reject_minimax=reject_regular=reject_freq=double(rep)
  
  #simulate adult data
  set.seed(25)
  data_a=rnorm(n_a,mu_a,sqrt(var_a))
  sd_a=sd(data_a)
  #sd_a=sqrt(var_a)
  mean_a=mean(data_a)
  mean_a
  for (i in 1:rep){
	gamma=1/2
    
    #simulate pediatric data
    #n_p=ceiling(p*n_a)
    data_p=rnorm(n_p,mu_p,sqrt(var_p))
    mean_p=mean(data_p)
    #sd_p=sqrt(var_p)
    sd_p=sd(data_p)
    
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
    #rnMix <- robustify(nm, weight=0.1, mean=0, n=1, sigma=sqrt(var_a))
    rnMix <- robustify(nm, weight=0.5, mean=0, n=1, sigma=sd_a)
    posterior.sum <- postmix(rnMix, m=mean_p, n=n_p, sigma=sd_p)
    #posterior.sum <- postmix(rnMix, data=data_p)
    # indicator <- runif(n.samples) <= posterior.sum[1,1]
    # theta_mixture <- double(n.samples)
    # for (j in 1:n.samples){
    #   if (indicator[i]){
    #     theta_mixture[j] <- rnorm(1,posterior.sum[2,1],posterior.sum[3,1])
    #   }else{
    #     theta_mixture[j] <- rnorm(1,posterior.sum[2,2],posterior.sum[3,2])
    #   }
    # }
    components <- sample(1:2,size=n.samples,prob=posterior.sum[1,],replace=TRUE)
    mus <- posterior.sum[2,]
    sds <- posterior.sum[3,]
    theta_mixture <- rnorm(n.samples,mean=mus[components],sd=sds[components])
    
    #estimate probability of getting theta estimate that is greater than zero
    reject_mixture[i]=ifelse(mean(theta_mixture<0)<=alpha,1,0)
    
    reject_regular[i]=ifelse(mean(theta_regular<0)<=alpha,1,0)
    reject_freq[i]=ifelse(mean_p*sqrt(n_p)/sd_p>qnorm(1-alpha),1,0)
    
  }
  
  res=data.frame("mixture"=reject_mixture,"minimax"=reject_minimax,"regular"=reject_regular,"freq"=reject_freq)
  return(colMeans(res))
  
}

###Collecting results
#3.1 
intensity=5000
list_sample_size=list(c(1000,1000),c(1000,800),c(1000,600),c(1000,400),c(1000,200),c(1000,100),c(1000,50),c(1000,25))

###Under null hypothesis
#Scenario 0: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=-0.05,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_0=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario I: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_1=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario II: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)

SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_2=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario III: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)

### investigating on small sample sizes of pediatric data
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_3=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario IV: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=20^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_4=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


###Under alternative hypothesis
#Scenario 0: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_0=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario I: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_1=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario II: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_2=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario III: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=20^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_3=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario IV: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_4=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario V: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_5=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario VI: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_6=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario VII: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=20^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_7=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario VIII: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_8=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario IX: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_9=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

result_H0_0
result_H0_1
result_H0_2
result_H0_3
result_H0_4
result_H1_0
result_H1_1
result_H1_2
result_H1_3
result_H1_4
result_H1_5
result_H1_6
result_H1_7
result_H1_8
result_H1_9

save.image("conditional_dataout_half_55.Rdata")


#################
###Continuous endpoints
Bayes_continuous=function(mu_a,var_a,n_a,n_p,mu_p,var_p,n.samples,alpha,rep){
  
  #initialize res vectors
  reject_mixture=reject_minimax=reject_regular=reject_freq=double(rep)
  
  #simulate adult data
  set.seed(25)
  data_a=rnorm(n_a,mu_a,sqrt(var_a))
  sd_a=sd(data_a)
  #sd_a=sqrt(var_a)
  mean_a=mean(data_a)
  mean_a
  for (i in 1:rep){
	gamma=1
    
    #simulate pediatric data
    #n_p=ceiling(p*n_a)
    data_p=rnorm(n_p,mu_p,sqrt(var_p))
    mean_p=mean(data_p)
    #sd_p=sqrt(var_p)
    sd_p=sd(data_p)
    
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
    #rnMix <- robustify(nm, weight=0.1, mean=0, n=1, sigma=sqrt(var_a))
    rnMix <- robustify(nm, weight=0.5, mean=0, n=1, sigma=sd_a)
    posterior.sum <- postmix(rnMix, m=mean_p, n=n_p, sigma=sd_p)
    #posterior.sum <- postmix(rnMix, data=data_p)
    # indicator <- runif(n.samples) <= posterior.sum[1,1]
    # theta_mixture <- double(n.samples)
    # for (j in 1:n.samples){
    #   if (indicator[i]){
    #     theta_mixture[j] <- rnorm(1,posterior.sum[2,1],posterior.sum[3,1])
    #   }else{
    #     theta_mixture[j] <- rnorm(1,posterior.sum[2,2],posterior.sum[3,2])
    #   }
    # }
    components <- sample(1:2,size=n.samples,prob=posterior.sum[1,],replace=TRUE)
    mus <- posterior.sum[2,]
    sds <- posterior.sum[3,]
    theta_mixture <- rnorm(n.samples,mean=mus[components],sd=sds[components])
    
    #estimate probability of getting theta estimate that is greater than zero
    reject_mixture[i]=ifelse(mean(theta_mixture<0)<=alpha,1,0)
    
    reject_regular[i]=ifelse(mean(theta_regular<0)<=alpha,1,0)
    reject_freq[i]=ifelse(mean_p*sqrt(n_p)/sd_p>qnorm(1-alpha),1,0)
    
  }
  
  res=data.frame("mixture"=reject_mixture,"minimax"=reject_minimax,"regular"=reject_regular,"freq"=reject_freq)
  return(colMeans(res))
  
}

###Collecting results
#3.1 
intensity=5000
list_sample_size=list(c(1000,1000),c(1000,800),c(1000,600),c(1000,400),c(1000,200),c(1000,100),c(1000,50),c(1000,25))

###Under null hypothesis
#Scenario 0: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=-0.05,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_0=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario I: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_1=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario II: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)

SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_2=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario III: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)

### investigating on small sample sizes of pediatric data
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_3=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario IV: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=20^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_4=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


###Under alternative hypothesis
#Scenario 0: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_0=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario I: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_1=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario II: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_2=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario III: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=20^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_3=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario IV: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=5^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_4=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario V: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_5=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario VI: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=15^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_6=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario VII: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=20^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_7=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario VIII: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_8=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario IX: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1.5,var_p=10^2,n.samples=intensity,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_9=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

result_H0_0
result_H0_1
result_H0_2
result_H0_3
result_H0_4
result_H1_0
result_H1_1
result_H1_2
result_H1_3
result_H1_4
result_H1_5
result_H1_6
result_H1_7
result_H1_8
result_H1_9

save.image("conditional_dataout_one_55.Rdata")
#################


###Discrete endpoint
###Continuous endpoints
Bayes_discrete=function(n_a,p_a,n_p,p_p,n.samples,p_target,alpha,rep){
  
  #simulate adult data
  data_a=rbinom(n_a,1,p_a)
  mean_a=mean(data_a)
  
  #initialize res vectors
  reject_minimax=reject_regular=reject_freq=double(rep)
  
  for (i in 1:rep){
    
    #simulate pediatric data
    data_p=rbinom(n_p,1,p_p)
    mean_p=mean(data_p)
    
    
    #simulate parameter of interest theta
    #minimax
    p_alpha_minimax=(n_a+n_p-1)*mean_p
    
    p_beta_minimax=(n_a+n_p-1)*(1-mean_p)
    
    #regular
    p_alpha_regular=(n_a-1)*mean_a+(sum(data_p))
    p_beta_regular=(n_a-1)*(1-mean_a)+n_p-sum(data_p)
    
    
    #simulating
    p_minimax=rbeta(n.samples,p_alpha_minimax,p_beta_minimax)
    p_regular=rbeta(n.samples,p_alpha_regular,p_beta_regular)
    
    #estimate probability of getting theta estimate that is greater than zero
    reject_minimax[i]=ifelse(mean(p_minimax<p_target)<=alpha,1,0)
    reject_regular[i]=ifelse(mean(p_regular<p_target)<=alpha,1,0)
    reject_freq[i]=ifelse((mean_p-p_target)/sqrt(mean_p*(1-mean_p)/n_p)>1.645,1,0)
    
  }
  
  res=data.frame("minimax"=reject_minimax,"regular"=reject_regular,"freq"=reject_freq)
  return(colMeans(res))
  
}


###Collecting results
#Bayes_discrete(n_a=1000,p_a=0.7,n_p=1000,p_p=0.7,n.samples=10^4,p_target=0.5,alpha=0.05,rep=10^4)
###Under alternative hypothesis
#Scenario I: there is treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_discrete(n_a=1000,p_a=0.6,n_p=x,p_p=0.6,n.samples=10^4,p_target=0.5,alpha=0.05,rep=10^4)
ss_p=c(1000,800,600,400,200,100,50,25)
SI=sapply(ss_p,res_I)
SI_res=do.call(rbind,SI)
result_SI=data.table(ss_p=ss_p,minimax=SI_res[,1],regular=SI_res[,2],frequentist=SI_res[,3])


#Scenario II: there is treatment effect for pediatric population, computes Type I error
res_II <- function(x)Bayes_discrete(n_a=1000,p_a=0.7,n_p=x,p_p=0.6,n.samples=10^4,p_target=0.5,alpha=0.05,rep=10^4)
ss_p=c(1000,800,600,400,200,100,50,25)
SII=sapply(ss_p,res_II)
SII_res=do.call(rbind,SII)
result_SII=data.table(ss_p=ss_p,minimax=SII_res[,1],regular=SII_res[,2],frequentist=SII_res[,3])


#Scenario III: there is treatment effect for pediatric population, computes Type I error
res_III <- function(x)Bayes_discrete(n_a=1000,p_a=0.7,n_p=x,p_p=0.5,n.samples=10^4,p_target=0.5,alpha=0.05,rep=10^4)
ss_p=c(1000,800,600,400,200,100,50,25)
SIII=sapply(ss_p,res_III)
SIII_res=do.call(rbind,SII)
result_SIII=data.table(ss_p=ss_p,minimax=SIII_res[,1],regular=SIII_res[,2],frequentist=SIII_res[,3])

res_all=cbind(t(SI),t(SII),t(SIII))

write.csv(res_all,"res_all.csv")
