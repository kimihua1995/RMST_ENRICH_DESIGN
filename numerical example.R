library(survRM2)
library(survival)
library(survminer)
library(tidyverse)
library(OpenMx)
library(eventTrack)
library(ggpubr)
library(ggplot2)


setwd("/Users/HuaKimi/Documents/Kaiyuan Hua PhD Dissertation/Topic 1/codes/DCC3")
source("source_numerical example.R")



##########################################################
# 1) parameters
lambda0 <- c(log(2)*2.5,log(2)*2.5)
lambda1 <- c(log(2)*6,log(2)*2.0)
beta0 <- 0
beta1 <- -0.8
t.change <- 1/6
n <- 250
tau <- 1.5
seed <- 212
t.enroll1 <- 0.5 # stage I
t.enrich <- 0.5
t.enroll2 <- 1 # stage II
t.FA <- 2.5 # follow-up
######################################################




# 2) Figure 2, KM curves of the simulated data
KM_plot_all <- KM_plot(seed, n, 
                       lambda0, lambda1, beta0, beta1, 
                       t.change, tau, 0.01, "PD-L1 >= 1%")

KM_plot_pos <- KM_plot(seed, n, 
                       lambda0, lambda1, beta0, beta1, 
                       t.change, tau, 0.5, "PD-L1 >= 50%")

pdf("KM_curves.pdf", width = 8, height = 6, onefile = T)
KM_plot_all$curve
KM_plot_pos$curve
dev.off()




################################################################
# 3) True biomarker cutpoint, RMST difference
true.cutoff <- true.cut.piece(tau,lambda0,lambda1,beta0,beta1,t.change)
(true.cutoff <- as.numeric(true.cutoff[1]))
(RMSTD_pos <- rmst.diff.true(tau,lambda0,lambda1,beta0,beta1,t.change,true.cutoff,1))
(RMSTD_all <- rmst.diff.true(tau,lambda0,lambda1,beta0,beta1,t.change,0.01,1))
#########################################################



################################################################
# 4) calculate asymptotic variance in Equation 4.4
################################################################
stat_pos <- est.sigma(seed, B=10000, lambda0, lambda1, beta0, beta1,
                   t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                   tau, true.cutoff)
sigma_pos <- mean(stat_pos$se_pos)*sqrt(10000)
  
stat_all <- est.sigma(seed, B=10000, lambda0, lambda1, beta0, beta1,
                     t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                    tau, 0.01)
sigma_all <- mean(stat_all$se_pos)*sqrt(10000)

true_RMSTreg <- true.RMSTreg(seed, n=5000, B=10000,
                             lambda0, lambda1, beta0, beta1, 
                              t.change, tau)
beta3 <- mean(true_RMSTreg$betas[,4])
sigma_beta3 <- mean(true_RMSTreg$vcovs[,4])*sqrt(10000)
  
  
  
  
  
################################################################
# 5) critical values q by Monte Carlo Method
################################################################
B <- 10000
res.null <- my.design.null(seed, B, lambda0, beta0, t.change,
                           t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, alpha0=0.025, true.cutoff)
  
alpha_list <- seq(0.015, 0.025, by = 0.001)
q_list <- qnorm(1 - alpha_list)
stat_null <- res.null$stat.FA$z_pos
rej <- rep(NA, length(alpha_list))
for (i in 1:length(q_list)){
  rej[i] <- sum(stat_null > q_list[i], na.rm = T)/B
}
### Pick the smallest one in alpha_list with rej >= 90%  





################################################################
# 6) Figure 3, plot the global power by required total sample size
################################################################
n_list <- seq(700,1000,by=5)
q <- qnorm(1-0.023)
q0 <- qnorm(1-0.025)
n_pos_list1 <- ceiling(n_list*(1-true.cutoff)/0.99)
n_pos_list2 <- ceiling(n_list/2*(1-true.cutoff)/0.99 + n_list/2)
eta_list <- 1 - pnorm(q0 - beta3/sigma_beta3*sqrt(n_list))
power1_N <- (1 - pnorm(q - RMSTD_pos/sigma_pos*sqrt(n_pos_list1)))*eta_list +
  (1 - pnorm(q - RMSTD_all/sigma_all*sqrt(n_list)))*(1-eta_list)
power2_N <- (1 - pnorm(q - RMSTD_pos/sigma_pos*sqrt(n_pos_list2)))*eta_list +
  (1 - pnorm(q - RMSTD_all/sigma_all*sqrt(n_list)))*(1-eta_list)

  
  
plot_dat <- data.frame(n = rep(n_list,2),
                       power = c(power1_N, power2_N),
                       Design = rep(c("All-Comer", "Enrich"), each = length(n_list)))
plot_dat$design <- factor(plot_dat$Design,
                          levels = c("All-Comer","Enrich"),
                          labels = c("All-Comer","Enrich"))
  
  
  library(RColorBrewer)
  f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
  (cols <- f("Dark2"))
  (cols <- f("Set1"))
  
  p <- ggplot(data = plot_dat, aes(x = n, y = power, color = Design)) +
    geom_line(size=1) +
    #scale_color_brewer(palette = "PuOr") +
    scale_color_manual(values=c("#377EB8","#FF7F00"),
                       breaks = c("All-Comer","Enrich"))+
    geom_hline(yintercept=0.9, linetype="dashed", color = "red", size = 1) +
    scale_x_continuous(breaks = seq(700,1000,50)) +
    labs(x = "Total Sample Size n", y = "Global Power") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  
  pdf("sample size.pdf", width = 6, height = 4, onefile = T)
  p
  dev.off()





