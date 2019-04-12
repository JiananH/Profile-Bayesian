
######################################################################
#### Functions used for Inference from Adult to Pediatric data V2 ####
####                 25 Mar 2019, Jianan Hui                      ####
######################################################################

##Three updates were made: 
# .	Change weight of non-informative prior from 0.5 to 0.1. 
# .	Change sample variance to population variance. 
# .	I added one scenario (5) in addition to what we have:
# H0-1.       Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=5^2)
# H0-2.       Under H0: (mu_a=1, var_a=10^2, mu_p=-0.05, var_p=5^2)
# H0-3.       Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=10^2)
# H0-4.       Under H0: (mu_a=1, var_a=10^2, mu_p=-0.05, var_p=10^2)
# H1-0.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=10^2)
# H1-1.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=15^2)
# H1-2.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=15^2)
# H1-3.       Under H1: (mu_a=1, var_a=10^2, mu_p=1, var_p=10^2)
# H1-4.       Under H1: (mu_a=1, var_a=10^2, mu_p=1.5, var_p=10^2)

###Simulating function###
set.seed(2019)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RBesT)

###Continuous endpoints
#Bayes_continuous=function(mu_a,var_a,n_a,n_p,mu_p,var_p,n.samples,alpha,rep){

#hyperparameters
rep=1000
n.samples=10000
alpha=0.025
mu_a=1
var_a=20^2
mu_p=0.8
var_p=20^2
n_a=1000
n_p=200

  #initialize res vectors
  reject_mixture=reject_minimax=reject_regular=reject_freq=double(rep)
  actual_info=matrix(data=0,nrow=rep,ncol=4)
  
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
    
    actual_info[i,] <- c(mean(data_a),sd(data_a),mean(data_p),sd(data_p))
    
  }
  
  data <- data.frame("mean_a"=actual_info[,1],"sd_a"=actual_info[,2],"mean_p"=actual_info[,3],"sd_p"=actual_info[,4],"mixture"=reject_mixture,"minimax"=reject_minimax,"regular"=reject_regular,"freq"=reject_freq)


p = ggplot() + 
    geom_point(data=data,aes(x=mean_a,y=mean_p),color=mixture) + labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0")))
  
