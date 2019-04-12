
######################################################################
#### Functions used for Inference from Adult to Pediatric data V2 ####
####                 25 Mar 2019, Jianan Hui                      ####
######################################################################

##Three updates were made: 
# .	Change weight of non-informative prior from 0.5 to 0.1. 
# .	Change sample variance to population variance. 
# .	I added one scenario (5) in addition to what we have:
# 1.       Under H0: (mu_a=1, var_a=10, mu_p=0, var_p=5)
# 2.       Under H0: (mu_a=1, var_a=10, mu_p=-0.05, var_p=5)
# 3.       Under H1: (mu_a=1, var_a=10, mu_p=0.5, var_p=10)
# 4.       Under H1: (mu_a=1, var_a=10, mu_p=0.8, var_p=15)
# 5.       Under H1: (mu_a=1, var_a=10, mu_p=1, var_p=10)
# 6.       Under H1: (mu_a=1, var_a=10, mu_p=1.5, var_p=10)

###Simulating function###
set.seed(2019)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RBesT)
setwd("//eu.boehringer.com/users/rdg/users4/jhui/Documents/Side projects/Minimax Adult to Pediatric/Minimax simulation/Add adjusted alpha and visualization on rejection region")

###Continuous endpoints
Bayes_continuous=function(mu_a,var_a,n_a,n_p,mu_p,var_p,n.samples,alpha,rep){
  
  #initialize res vectors
  reject_mixture=reject_minimax=reject_regular=reject_freq=double(rep)
  
  for (i in 1:rep){
    #simulate adult data
    data_a=rnorm(n_a,mu_a,sqrt(var_a))
    sd_a=sqrt(var_a)
    mean_a=mean(data_a)
    
    #simulate pediatric data
    #n_p=ceiling(p*n_a)
    data_p=rnorm(n_p,mu_p,sqrt(var_p))
    mean_p=mean(data_p)
    sd_p=sqrt(var_p)
    
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
    theta_minimax=rnorm(n.samples,mu_theta_minimax,sd_theta)
    theta_regular=rnorm(n.samples,mu_theta_regular,sd_theta)
    
    nm <- mixnorm(adult=c(mean_a, sd_a, sd_a/sqrt(n_a)), sigma=sd_a)
    rnMix <- robustify(nm, weight=0.1, mean=0, n=1, sigma=sqrt(var_a))
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
    components <- sample(1:2,prob=posterior.sum[1,])
    mus <- posterior.sum[2,]
    sds <- posterior.sum[3,]
    theta_mixture <- rnorm(n.samples,mean=mus[components],sd=sds[components])
    
    #estimate probability of getting theta estimate that is greater than zero
    reject_mixture[i] <- ifelse(mean(theta_mixture<0)<=alpha,1,0)
    reject_minimax[i] <- ifelse(mean(theta_minimax<0)<=alpha,1,0)
    reject_regular[i] <- ifelse(mean(theta_regular<0)<=alpha,1,0)
    reject_freq[i] <- ifelse(mean_p*sqrt(n_p)/sd_p>qnorm(1-alpha),1,0)
    
  }
  
  res <- data.frame("mixture"=reject_mixture,"minimax"=reject_minimax,"regular"=reject_regular,"freq"=reject_freq)
  theta <- data.frame("mixture"=theta_mixture,"minimax"=theta_minimax,"regular"=theta_regular)
  theta_freq <- data_p
  return(list("res"=colMeans(res),"theta"=theta,"theta_freq"=theta_freq))

  
  
}


###Collecting results
#3.1 
intensity=1
n.samples=10000
list_sample_size=list(c(1000,1000),c(1000,800),c(1000,600),c(1000,400),c(1000,200))

###@ Yaoshi: Suggestion: use same sets of sample sizes across different scenarios.

###Under null hypothesis
#Scenario I: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=5^2,n.samples=n.samples,alpha=0.025,rep=intensity)
object <- lapply(list_sample_size,res_I)[[5]]
SI <- object$theta
SII <- data.frame("frequentist"=object$theta_freq)
p1 <- ggplot(SI,aes(x=mixture))+geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
  geom_density(alpha=.2, fill="#FF6666")+geom_vline(aes(xintercept=0, color="red"),
                                                    linetype="dashed") + theme(legend.position="none")
p2 <- ggplot(SI,aes(x=minimax))+geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
  geom_density(alpha=.2, fill="#FF6666")+geom_vline(aes(xintercept=0, color="red"),
                                                    linetype="dashed") + theme(legend.position="none")
p3 <- ggplot(SI,aes(x=regular))+geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
  geom_density(alpha=.2, fill="#FF6666")+geom_vline(aes(xintercept=0, color="red"),
                                                    linetype="dashed") + theme(legend.position="none")
p4 <- ggplot(SII,aes(x=frequentist))+geom_histogram(aes(y=..density..), colour="black", fill="white",bins=5)+
  geom_density(alpha=.2, fill="#FF6666")+geom_vline(aes(xintercept=10*qnorm(0.975), color="red"),
                                                    linetype="dashed") + theme(legend.position="none")
grid.arrange(p1, p2, p3, p4, nrow=2)+ theme(legend.position="none")
# SI_res=do.call(rbind,SI)# result_H0_1=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario II: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=-0.05,var_p=5^2,n.samples=n.samples,alpha=0.025,rep=intensity)

SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_2=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario III: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0,var_p=10^2,n.samples=n.samples,alpha=0.025,rep=intensity)

### investigating on small sample sizes of pediatric data
list_sample_size=list(c(1000,1000),c(1000,800),c(1000,600),c(1000,400),c(1000,200),c(1000,100),c(1000,50),c(1000,25))
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_3=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario IV: there is no treatment effect for pediatric population, computes Type I error
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=-0.05,var_p=10^2,n.samples=n.samples,alpha=0.025,rep=intensity)

SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H0_4=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


###Under alternative hypothesis
#Scenario I: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.5,var_p=10^2,n.samples=n.samples,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_0=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario II: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=10^2,n.samples=n.samples,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_1=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario III: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=0.8,var_p=15^2,n.samples=n.samples,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_2=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])

#Scenario IV: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1,var_p=10^2,n.samples=n.samples,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_3=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


#Scenario V: there is treatment effect for pediatric population, computes power
res_I <- function(x)Bayes_continuous(mu_a=1,var_a=10^2,n_a=x[1],n_p=x[2],mu_p=1.5,var_p=10^2,n.samples=n.samples,alpha=0.025,rep=intensity)
SI=lapply(list_sample_size,res_I)
SI_res=do.call(rbind,SI)
result_H1_4=data.table(SampleSize_a_p=list_sample_size,mixture=SI_res[,1],minimax=SI_res[,2],regular=SI_res[,3],frequentist=SI_res[,4])


result_H0_1
result_H0_2
result_H0_3
result_H0_4
result_H1_0
result_H1_1
result_H1_2
result_H1_3
result_H1_4



save.image("datain.Rdata")