#Create plot for illustration on rejection region



#data summarizing
#sample size of pediatric group
result_H0_1$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0_1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p1 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))
p1

#s2
result_H0_2$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0_2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p2 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = -0.05")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))p2

#s3
result_H1_1$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1_1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p3 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))p3

#s4
result_H1_2$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1_2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p4 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))p4


p = ggplot() + 
  geom_line(data=data,aes(x=ss_p,y=minimax),color="royalblue2") +
  geom_point(data=data,aes(x=ss_p,y=minimax),color="royalblue2") +
  geom_line(data=data,aes(x=ss_p,y=regular),color="darkorange1") +
  geom_point(data=data,aes(x=ss_p,y=regular),color="darkorange1") +
  geom_line(data=data,aes(x=ss_p,y=frequentist),color="forestgreen") +
  geom_point(data=data,aes(x=ss_p,y=frequentist),color="forestgreen") +
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0")))

p

p = ggplot(data=data,aes())

# #Scenario IV: there is treatment effect for pediatric population, computes Type I error
# res_IV <- function(x)Bayes_discrete(n_a=1000,p_a=0.7,n_p=x,p_p=0.6,n.samples=10^4,p_target=0.5,alpha=0.05,rep=10^4)
# ss_p=c(1000,800,600,400,200,100,50,25)
# SIV=sapply(ss_p,res_IV)
# SIV_res=do.call(rbind,SIV)
# result_SIV=data.table(ss_p=ss_p,minimax=SIV_res[,1],regular=SIV_res[,2],frequentist=SIV_res[,3])



###inloop###
#sample size of pediatric group
result_H0_1$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0_1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p1 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))
p1

#s2
result_H0_2$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H0_2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p2 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = -0.05")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))p2

#s3
result_H1_1$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1_1
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p3 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))p3

#s4
result_H1_2$ss_p=c(1000,800,600,400,200,100,50,25)
data_wide=result_H1_2
data=reshape(data_wide,direction="long",varying=list(names(data_wide)[2:4]),v.names="Value",idvar="ss_p")
data$SampleSize_a_p=NULL
data$time=factor(data$time)
p4 = ggplot(data=data,aes(x=ss_p,y=Value,group=time,color=time))+
  geom_line()+
  geom_point()+
  #geom_hline(yintercept=0.05,linetype="dashed",color="darkgrey")+
  labs(x = "Pediatric Sample size", y = "Power", title = expression(paste("Under alternative hypothesis that ",mu[p]," = 1")))+
  #theme(legend.position="none")+
  scale_colour_discrete(name="Method",breaks=c(1,2,3,4),labels=c("Robust prior","Conditional Bayesian","Bayesian","Frequentist"))p4


p = ggplot() + 
  geom_line(data=data,aes(x=ss_p,y=minimax),color="royalblue2") +
  geom_point(data=data,aes(x=ss_p,y=minimax),color="royalblue2") +
  geom_line(data=data,aes(x=ss_p,y=regular),color="darkorange1") +
  geom_point(data=data,aes(x=ss_p,y=regular),color="darkorange1") +
  geom_line(data=data,aes(x=ss_p,y=frequentist),color="forestgreen") +
  geom_point(data=data,aes(x=ss_p,y=frequentist),color="forestgreen") +
  labs(x = "Pediatric Sample size", y = "Type I Error", title = expression(paste("Under null hypothesis that ",mu[p]," = 0")))

p

p = ggplot(data=data,aes())