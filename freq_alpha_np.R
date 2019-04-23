###Continuous endpoints
load("datain.Rdata")
freq_reject=function(mu_a,var_a,n_a,n_p,mu_p,var_p,n.samples,alpha,rep){
  
  #initialize res vectors
  reject_freq=double(rep)
  
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
    
    
    reject_freq[i]=ifelse(mean_p*sqrt(n_p)/sd_p>qnorm(1-alpha),1,0)
    
  }
  return(mean(reject_freq))
}
ss_p=c(1000,800,600,400,200,100,50,25)
minimax_alpha_5=result_H0_1$minimax
minimax_alpha_10=result_H0_2$minimax
minimax_alpha_15=result_H0_3$minimax
minimax_alpha_20=result_H0_4$minimax
intensity=5000

#H0-3. Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=10^2)

#H1-0. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=5^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_0=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_5)

#H1-1. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_1=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

#H1-2. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=15^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_2=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_15)

#H1-3. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=20^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_3=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_20)

#H1-4. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=5^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_4=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_5)

#H1-5. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_5=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

#H1-6. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=15^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_6=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_15)

#H1-7. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=20^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_7=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_20)

#H1-8. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=1,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_8=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

#H1-9. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=1.5,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_9=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

result_H1_0$frequentist_inflate=H1_0
result_H1_1$frequentist_inflate=H1_1
result_H1_2$frequentist_inflate=H1_2
result_H1_3$frequentist_inflate=H1_3
result_H1_4$frequentist_inflate=H1_4
result_H1_5$frequentist_inflate=H1_5
result_H1_6$frequentist_inflate=H1_6
result_H1_7$frequentist_inflate=H1_7
result_H1_8$frequentist_inflate=H1_8
result_H1_9$frequentist_inflate=H1_9

save.image("add frequentist_inflate_inloop.Rdata")

rm(ls())
###Continuous endpoints
load("dataout.Rdata")
freq_reject=function(mu_a,var_a,n_a,n_p,mu_p,var_p,n.samples,alpha,rep){
  
  #initialize res vectors
  reject_freq=double(rep)
  
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
    
    
    reject_freq[i]=ifelse(mean_p*sqrt(n_p)/sd_p>qnorm(1-alpha),1,0)
    
  }
  return(mean(reject_freq))
}
ss_p=c(1000,800,600,400,200,100,50,25)
minimax_alpha_5=result_H0_1$minimax
minimax_alpha_10=result_H0_2$minimax
minimax_alpha_15=result_H0_3$minimax
minimax_alpha_20=result_H0_4$minimax
intensity=5000

#H1-0. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=5^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_0=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_5)

#H1-1. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_1=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

#H1-2. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=15^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_2=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_15)

#H1-3. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.5,var_p=20^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_3=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_20)

#H1-4. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=5^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_4=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_5)

#H1-5. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_5=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

#H1-6. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=15^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_6=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_15)

#H1-7. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=0.8,var_p=20^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_7=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_20)

#H1-8. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=1,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_8=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

#H1-9. 
freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=1,var_a=10^2,n_a=1000,n_p,mu_p=1.5,var_p=10^2,n.samples=intensity,alpha,rep=intensity)
} 

H1_9=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha_10)

result_H1_0$frequentist_inflate=H1_0
result_H1_1$frequentist_inflate=H1_1
result_H1_2$frequentist_inflate=H1_2
result_H1_3$frequentist_inflate=H1_3
result_H1_4$frequentist_inflate=H1_4
result_H1_5$frequentist_inflate=H1_5
result_H1_6$frequentist_inflate=H1_6
result_H1_7$frequentist_inflate=H1_7
result_H1_8$frequentist_inflate=H1_8
result_H1_9$frequentist_inflate=H1_9

save.image("add frequentist_inflate_outloop.Rdata")