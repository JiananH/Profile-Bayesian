#Produce plots instead of tables per reviewer comment
setwd("~/Dropbox/Side projects/Profile Bayesian/Profile-Bayesian")
load("conditional_dataout_half_19.Rdata")
#load("~/Side projects/Minimax Adult to Pediatric/Minimax simulation/Profile-Bayesian-master/Profile-Bayesian-master/conditional_dataout_half_19.Rdata")

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

#load("~/Side projects/Minimax Adult to Pediatric/Minimax simulation/Profile-Bayesian-master/Profile-Bayesian-master/conditional_dataout_half_55.Rdata")
load("conditional_dataout_half_55.Rdata")

result_H0.1$mixture55 <- result_H0_0$mixture
result_H0.2$mixture55 <- result_H0_1$mixture
result_H0.3$mixture55 <- result_H0_2$mixture
result_H0.4$mixture55 <- result_H0_3$mixture
result_H0.5$mixture55 <- result_H1_1$mixture

result_H1.1$mixture55 <- result_H1_4$mixture
result_H1.2$mixture55 <- result_H1_5$mixture
result_H1.3$mixture55 <- result_H1_6$mixture
result_H1.4$mixture55 <- result_H1_8$mixture
result_H1.5$mixture55 <- result_H1_9$mixture


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
#png(filename = "null.png",width = 860, height = 600)
ggarrange(p2, p3, p4, p5, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
#dev.off()
ggarrange(p6, p8, p9, p10, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
