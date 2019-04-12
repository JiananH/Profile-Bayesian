###Continuous endpoints
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
minimax_alpha=result_H0_3$minimax

freq_reject_alpha=function(n_p,alpha){
  freq_reject(mu_a=0,var_a=10,n_a=1000,n_p,mu_p=0,var_p=10,n.samples=10000,alpha,rep=10000)
} 

res=mapply(function(x,y){freq_reject_alpha(n_p=x,alpha=y)},ss_p,minimax_alpha)
