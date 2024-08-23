library(survRM2)
library(survival)
library(survminer)
library(tidyverse)
library(OpenMx)
library(eventTrack)
setwd("/Users/HuaKimi/Documents/Kaiyuan Hua PhD Dissertation/Topic 1/codes/github_RMST_ENRICH_DESIGN")
source("source_simulation.R")
source("source_cutpoint.R")



lambda0 <- c(1.2,1.2)
lambda1 <- c(1.2,0.5)
beta1 <- 1.0
beta0 <- 0
t.change <- 0.25
t.enroll1 <- 1
t.enrich <- 1
t.enroll2 <- 2
t.FA <- 4
tau <- 2
n1 <- 300
n2 <- 300 # per arm per stage
alpha1 <- 0.025
alpha2 <- 0.025



true.cutoff <- true.cut.piece0(tau,lambda0,lambda1,beta0,beta1,t.change)
(true.cutoff <- as.numeric(true.cutoff[1]))
(RMSTD_pos <- rmst.true(lambda1,beta1,t.change,tau,true.cutoff, 1) - 
    rmst.true(lambda0,beta0,t.change,tau,true.cutoff, 1))
(RMSTD_all <- rmst.true(lambda1,beta1,t.change,tau, 0, 1) - 
    rmst.true(lambda0,beta0,t.change,tau, 0, 1))


CP <- function(d_pos, se_pos){
  d_l <- d_pos - qnorm(0.975)*se_pos
  d_h <- d_pos + qnorm(0.975)*se_pos
  cp <- (d_l <= RMSTD_pos) * (d_h >= RMSTD_pos)
  return(mean(cp, na.rm=T))
}





# 1) Alternative Setting
seed = 513
B = 10000
res_a_s1 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, alpha1, alpha2, true.cutoff, design = 1, pred.method = 1)
#mean(as.numeric(res_a_s1$cutpoint.IA[,1]), na.rm = T)
mean(as.numeric(res_a_s1$cutpoint.FA[,1]), na.rm = T)
sd(as.numeric(res_a_s1$cutpoint.FA[,1]), na.rm = T)

mean(res_a_s1$stat.FA1$d_pos, na.rm = T)
CP(res_a_s1$stat.FA1$d_pos, res_a_s1$stat.FA1$se_pos)

mean(res_a_s1$rej1)
mean(res_a_s1$N.neg)


res_a_s2 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, alpha1, alpha2, true.cutoff, design = 2, pred.method = 1)
mean(as.numeric(res_a_s2$cutpoint.IA[,1]), na.rm = T)
sd(as.numeric(res_a_s2$cutpoint.IA[,1]), na.rm = T)
mean(as.numeric(res_a_s2$cutpoint.FA[,1]), na.rm = T)
sd(as.numeric(res_a_s2$cutpoint.FA[,1]), na.rm = T)

mean(res_a_s2$stat.FA1$d_pos, na.rm = T)
CP(res_a_s2$stat.FA1$d_pos, res_a_s2$stat.FA1$se_pos)

mean(res_a_s2$rej1)
mean(res_a_s2$N.neg)



res_a_s3 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, alpha1, alpha2, true.cutoff, design = 2, pred.method = 2)
mean(as.numeric(res_a_s3$cutpoint.IA[,1]), na.rm = T)
sd(as.numeric(res_a_s3$cutpoint.IA[,1]), na.rm = T)
mean(as.numeric(res_a_s3$cutpoint.FA[,1]), na.rm = T)
sd(as.numeric(res_a_s3$cutpoint.FA[,1]), na.rm = T)

mean(res_a_s3$stat.FA1$d_pos, na.rm = T)
CP(res_a_s3$stat.FA1$d_pos, res_a_s3$stat.FA1$se_pos)

mean(res_a_s3$rej1)
mean(res_a_s3$N.neg)


res_a_s4 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, alpha1, alpha2, true.cutoff, design = 2, pred.method = 3)
mean(as.numeric(res_a_s4$cutpoint.IA[,1]), na.rm = T)
sd(as.numeric(res_a_s4$cutpoint.IA[,1]), na.rm = T)
mean(as.numeric(res_a_s4$cutpoint.FA[,1]), na.rm = T)
sd(as.numeric(res_a_s4$cutpoint.FA[,1]), na.rm = T)

mean(res_a_s4$stat.FA1$d_pos, na.rm = T)
CP(res_a_s4$stat.FA1$d_pos, res_a_s4$stat.FA1$se_pos)

mean(res_a_s4$rej1)
mean(res_a_s4$N.neg)





# 2) Null Setting
alpha1.list <- c(0.015,0.020,0.025,0.030)
alpha2.list <- c(0.020,0.023,0.025)
seed = 333
B = 10000
res_0_s1 <- my.design.null(seed, B, n1, n2, lambda0, lambda0, beta0, beta0,
                           t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, alpha1.list, alpha2.list, 
                           design = 1, pred.method = 1)
colMeans(res_0_s1$rej1, na.rm = T)

res_0_s2 <- my.design.null(seed, B, n1, n2, lambda0, lambda0, beta0, beta0,
                           t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, alpha1.list, alpha2.list, 
                           design = 2, pred.method = 1)
colMeans(res_0_s2$rej1, na.rm = T)


res_0_s4 <- my.design.null(seed, B, n1, n2, lambda0, lambda0, beta0, beta0,
                           t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, alpha1.list, alpha2.list, 
                           design = 2, pred.method = 3)
colMeans(res_0_s4$rej1, na.rm = T)
