KM_plot <- function(seed, n, 
                    lambda0, lambda1, beta0, beta1, 
                    t.change, tau, threshold, title){
  set.seed(seed)
  X0 <- runif(n,threshold,1); X1 <- runif(n,threshold,1)
  u0 <- runif(n,0,1); u1 <- runif(n,0,1)
  time0 <- -log(u0)*exp(-X0*beta0)/lambda0[1]
  pos0 <- time0 > t.change
  time0[pos0] <- -log(u0[pos0])*exp(-X0[pos0]*beta0)/lambda0[2] + (1 - lambda0[1]/lambda0[2])*t.change
  time1 <- -log(u1)*exp(-X1*beta1)/lambda1[1]
  pos1 <- time1 > t.change
  time1[pos1] <- -log(u1[pos1])*exp(-X1[pos1]*beta1)/lambda1[2] + (1 - lambda1[1]/lambda1[2])*t.change
  time <- c(time0,time1)
  y <- pmin(time, tau)
  X <- c(X0,X1)
  Z <- c(rep(0,n),rep(1,n))
  status <- as.numeric(time <= tau)
  KM_dat <- data.frame(X = X, Z = Z, time = time, time_m = time*12,
                       status = status, y = y)
  fit_KM <- survfit(Surv(time_m, status) ~ Z, data = KM_dat)
  
  
  curve <- ggsurvplot(fit_KM,data = KM_dat,size = 0.5,
                      censor.size = 3,conf.int = F,
                      #conf.int.style = "step",
                      censor.shape = "|",
                      palette = c("blue","red"),
                      pval = F, #pval.coord = c(0,0.15),
                      xlim = c(0,24),
                      break.time.by = 3,
                      xlab = "Month",
                      ylab = "Progression-Free Survival",
                      ylim = c(0,1),
                      break.y.by = 0.1,
                      title = title,
                      surv.plot.height = 5,
                      legend.labs = c("Docetaxel","Avelumab"),
                      legend.title = "",
                      # risk table
                      risk.table = FALSE,
                      risk.table.fontsize = 3,
                      cumevents = F,
                      tables.height = 0.15,
                      cumevents.height = 0.15,
                      fontsize = 3.0)
  #curve$table <- curve$table + labs(x = NULL) + theme(axis.text.y = element_blank())

  
  return(list(curve = curve, KM_dat = KM_dat))
}




true.RMSTreg <- function(seed, n, B,
                         lambda0, lambda1, beta0, beta1, 
                         t.change, tau){
  set.seed(seed)
  true_betas <- true_vcovs <- matrix(NA, B, 4)
  for (b in 1:B){
    X0 <- runif(n,0.01,1); X1 <- runif(n,0.01,1)
    u0 <- runif(n,0,1); u1 <- runif(n,0,1)
    time0 <- -log(u0)*exp(-X0*beta0)/lambda0[1]
    pos0 <- time0 > t.change
    time0[pos0] <- -log(u0[pos0])*exp(-X0[pos0]*beta0)/lambda0[2] + (1 - lambda0[1]/lambda0[2])*t.change
    time1 <- -log(u1)*exp(-X1*beta1)/lambda1[1]
    pos1 <- time1 > t.change
    time1[pos1] <- -log(u1[pos1])*exp(-X1[pos1]*beta1)/lambda1[2] + (1 - lambda1[1]/lambda1[2])*t.change
    time <- c(time0,time1)
    y <- pmin(time, tau)
    X <- c(X0,X1)
    Z <- c(rep(0,n),rep(1,n))
    status <- as.numeric(time <= tau)
    fit_rmst_true <- lm(y ~ Z + X + Z:X)
    true_betas[b,] <- coef(fit_rmst_true)
    true_vcovs[b,] <- sqrt(diag(vcov(fit_rmst_true)))
  }
  
  return(list(betas = true_betas, vcovs = true_vcovs))
}




true.cut.piece <- function(tau,lambda0,lambda1,beta0,beta1,t.change){
  mu0 <- function(x,t){
    if (t <= t.change){
      integrate(function(y){exp(-lambda0[1]*exp(beta0*x)*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda0[1]*exp(beta0*x)*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda0[1]*exp(beta0*x)*t.change-lambda0[2]*exp(beta0*x)*(y-t.change))},
                  lower = t.change, upper = t, rel.tol = 1e-10)$value
    }
  }
  mu1 <- function(x,t){
    if (t <= t.change){
      integrate(function(y){exp(-lambda1[1]*exp(beta1*x)*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda1[1]*exp(beta1*x)*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda1[1]*exp(beta1*x)*t.change-lambda1[2]*exp(beta1*x)*(y-t.change))},
                  lower = t.change, upper = t, rel.tol = 1e-10)$value
    }
  }
  
  if ((mu1(0,tau) - mu0(0,tau))*(mu1(1,tau) - mu0(1,tau)) < 0){
    cut_t <- uniroot(function(x){mu1(x,tau)-mu0(x,tau)},
                     interval = c(0,1), tol = 1e-10)$root
    if ((mu1(0,tau) - mu0(0,tau)) > 0){
      direction <- "left"}else {
        direction <- "right"}
  }
  
  if ((mu1(0,tau) - mu0(0,tau)) > 0 & 
      (mu1(1,tau) - mu0(1,tau)) > (mu1(0,tau) - mu0(0,tau))){
    cut_t <- 0
    direction <- "right"
  }
  if ((mu1(1,tau) - mu0(1,tau)) > 0 & 
      (mu1(1,tau) - mu0(1,tau)) < (mu1(0,tau) - mu0(0,tau))){
    cut_t <- 1
    direction <- "left"
  }
  if ((mu1(0,tau) - mu0(0,tau)) < 0 & 
      (mu1(1,tau) - mu0(1,tau)) < 0){
    cut_t <- NA
    direction <- "all negative"
  }
  
  
  return(c(cut_t, direction))
}




rmst.diff.true <- function(tau,lambda0,lambda1,beta0,beta1,t.change, 
                           lower, upper){
  mu0 <- function(x, tau){
    ifelse(tau <= t.change,
           integrate(function(y){exp(-lambda0[1]*exp(beta0*x)*y)},
                     lower = 0, upper = tau, rel.tol = 1e-10)$value,
           integrate(function(y){exp(-lambda0[1]*exp(beta0*x)*y)},
                     lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
             integrate(function(y){exp(-lambda0[1]*exp(beta0*x)*t.change-lambda0[2]*exp(beta0*x)*(y-t.change))},
                       lower = t.change, upper = tau, rel.tol = 1e-10)$value)
  }
  
  mu1 <- function(x, tau){
    ifelse(tau <= t.change,
           integrate(function(y){exp(-lambda1[1]*exp(beta1*x)*y)},
                     lower = 0, upper = tau, rel.tol = 1e-10)$value,
           integrate(function(y){exp(-lambda1[1]*exp(beta1*x)*y)},
                     lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
             integrate(function(y){exp(-lambda1[1]*exp(beta1*x)*t.change-lambda1[2]*exp(beta1*x)*(y-t.change))},
                       lower = t.change, upper = tau, rel.tol = 1e-10)$value)
  }
  
  integrate(Vectorize(function(x){(mu1(x,tau)-mu0(x,tau))/(upper-lower)}), 
            lower = lower, upper = upper, rel.tol = 1e-10)$value
}





sim.data.piece <- function(X,lambda,beta,t.change,t.enroll.start,t.enroll.end,t.IA){
  n <- length(X)
  u <- runif(n,0.00001,1)
  time.event <- -log(u)*exp(-X*beta)/lambda[1]
  pos <- time.event > t.change
  time.event[pos] <- -log(u[pos])*exp(-X[pos]*beta)/lambda[2] + (1 - lambda[1]/lambda[2])*t.change
  time.enroll <- runif(n, t.enroll.start, t.enroll.end)
  time.loss <- rexp(n, -log(1-0.05)/2)
  time.censor <- pmin(time.loss, t.IA - time.enroll)
  time <- pmin(time.event, time.censor)
  status <- as.numeric(time.event < time.censor)
  
  dat <- data.frame(X = X, time.event = time.event, time.enroll = time.enroll,
                    time.loss = time.loss, time.censor = time.censor, time = time, status = status)
  return(dat)
}




sample_size <- function(seed, B, lambda0, lambda1, beta0, beta1,
                        t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                        tau, alpha, power, true.cutoff){
  n1 <- 2500
  n2 <- 2500
  stat <- data.frame(z_pos = rep(NA,B), d_pos = rep(NA,B), se_pos = rep(NA,B))
  for (b in 1:B){
    if (!b%%1000) print(b)
    X0_1 <- runif(n1, true.cutoff, 1)
    X1_1 <- runif(n1, true.cutoff, 1)
    dat0_1 <- sim.data.piece(X0_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    dat1_1 <- sim.data.piece(X1_1,lambda1,beta1,t.change,0,t.enroll1,t.enrich)
    X0_2 <- runif(n2, true.cutoff, 1)
    X1_2 <- runif(n2, true.cutoff, 1)
    dat0_2 <- sim.data.piece(X0_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat1_2 <- sim.data.piece(X1_2,lambda1,beta1,t.change,t.enrich,t.enroll2,t.FA)
    dat0 <- rbind(dat0_1, dat0_2); dat0$A = 0
    dat1 <- rbind(dat1_1, dat1_2); dat1$A = 1
    
    
    dat0 <- dat.modify(dat0, t.FA)
    dat1 <- dat.modify(dat1, t.FA)
    
    time_pos=c(dat0$time, dat1$time)
    status_pos=c(dat0$status, dat1$status)
    arm_pos=c(dat0$A, dat1$A)
    
    fit_pos=rmst2(time_pos, status_pos, arm_pos, tau=tau)
    rmst1est_pos=fit_pos$RMST.arm1$rmst[1:2]
    rmst0est_pos=fit_pos$RMST.arm0$rmst[1:2]
    d_pos=(rmst1est_pos[1]-rmst0est_pos[1])
    se_pos=sqrt(rmst1est_pos[2]^2+rmst0est_pos[2]^2)
    z_pos=d_pos/se_pos
    
    stat[b,] <- c(z_pos,d_pos,se_pos)
  }
  
  
  RMSTD <- rmst.diff.true(tau,lambda0,lambda1,beta0,beta1,t.change,true.cutoff, 1)
  
  n.sim <- 2*(n1+n2)
  n <- ((qnorm(1-alpha) + qnorm(power))*mean(stat$se_pos)*sqrt(n.sim)/RMSTD)^2/2
  
  return(n)
}




dat.modify <- function(dat,t.IA){
  dat$time.censor <- pmin(dat$time.loss, t.IA - dat$time.enroll)
  dat$time <- pmin(dat$time.event, dat$time.censor)
  dat$status <- as.numeric(dat$time.event < dat$time.censor)
  return(dat)
}



est.sigma <- function(seed, B, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, true.cutoff){
  n1 <- 2500
  n2 <- 2500
  stat1 <- stat3 <- stat4 <- stat5 <- 
    data.frame(z_pos = rep(NA,B), d_pos = rep(NA,B), se_pos = rep(NA,B))
  for (b in 1:B){
    if (!b%%1000) print(b)
    X0_1 <- runif(n1, true.cutoff, 1)
    X1_1 <- runif(n1, true.cutoff, 1)
    dat0_1 <- sim.data.piece(X0_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    dat1_1 <- sim.data.piece(X1_1,lambda1,beta1,t.change,0,t.enroll1,t.enrich)
    X0_2 <- runif(n2, true.cutoff, 1)
    X1_2 <- runif(n2, true.cutoff, 1)
    dat0_2 <- sim.data.piece(X0_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat1_2 <- sim.data.piece(X1_2,lambda1,beta1,t.change,t.enrich,t.enroll2,t.FA)
    dat0 <- rbind(dat0_1, dat0_2); dat0$A = 0
    dat1 <- rbind(dat1_1, dat1_2); dat1$A = 1
    
    
    dat0 <- dat.modify(dat0, t.FA)
    dat1 <- dat.modify(dat1, t.FA)
    
    stat1[b,] <- stat1(dat0,dat1,tau)
  }
  
  return(stat1)
}










my.design.null <- function(seed, B, lambda0, beta0,
                           t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, alpha0, true.cutoff){
  
  cutpoint.FA.list <- data.frame(cutpoint = rep(NA,B), direction = rep(NA,B))
  stat.FA.list <- data.frame(z_pos = rep(NA,B), d_pos = rep(NA,B), se_pos = rep(NA,B))
  gamma3.list <- rep(NA,B)
  n1 <- n2 <- 2500
  set.seed(seed)
  for (b in 1:B){
    if (!b%%1000) {print(b)}
    #print(b)
    #set.seed(b)
    # Stage I
    X0_1 <- runif(n1, 0, 1)
    X1_1 <- runif(n1, 0, 1)
    dat0_1 <- sim.data.piece(X0_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    dat1_1 <- sim.data.piece(X1_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    X0_2 <- runif(n2, true.cutoff, 1)
    X1_2 <- runif(n2, true.cutoff, 1)
    dat0_2 <- sim.data.piece(X0_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat1_2 <- sim.data.piece(X1_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat0 <- rbind(dat0_1, dat0_2); dat0$A = 0
    dat1 <- rbind(dat1_1, dat1_2); dat1$A = 1
   
    
    
    
    #####################################
    ## FA
    dat0 <- dat.modify(dat0, t.FA)
    dat1 <- dat.modify(dat1, t.FA)
    dat <- rbind(dat0, dat1)
    
    ## estimate the cutoff point at FA
    rmst_fit_FA <- fit.rmst.reg(dat0, dat1, tau)
    gamma3 <- coef(rmst_fit_FA)[4]/sqrt(vcov(rmst_fit_FA)[4,4])
    gamma3.list[b] <- gamma3
    cutpoint.FA <- find.cutpoint(rmst_fit_FA)
    cutpoint.FA.list[b,] <- cutpoint.FA
    
  
    if (cutpoint.FA[2] == "all negative"){
      cutpoint.FA.list[b,] <- c(NA, "all negative")
    }else if (gamma3 <= qnorm(1 - alpha0)){
      stat.FA <- stat1(dat0, dat1, tau)
      stat.FA.list[b,] <- stat.FA
    }else {
      #dat0_pos <- dat0[dat0$X >= as.numeric(cutpoint.FA[1]),]
      #dat1_pos <- dat1[dat1$X >= as.numeric(cutpoint.FA[1]),]
      dat0_pos <- dat0[dat0$X >= true.cutoff,]
      dat1_pos <- dat1[dat1$X >= true.cutoff,]
      stat.FA <- stat1(dat0_pos, dat1_pos, tau)
      stat.FA.list[b,] <- stat.FA
    }

    
   
    
  }
  
  return(list(cutpoint.FA = cutpoint.FA.list,
              stat.FA = stat.FA.list,
              gamma3 = gamma3.list))
  
}








my.func_surv <- function(y, d){
  #--input--
  #y=time
  #d=status
  
  #--
  id=order(y)
  y=y[id]
  d=d[id]
  
  #--
  t_idx = unique(c(0,y))
  ny = length(y)
  
  #--
  Y = N = C = S = H = D = E = rep(0,length(t_idx))
  
  #i=1
  Y[1] = ny
  N[1] = 0
  C[1] = 0
  S[1] = 1
  H[1] = 0
  D[1] = 0
  E[1] = 0
  
  #i>=2
  for(i in 2:length(t_idx)){
    Y[i] = Y[i-1] - N[i-1] - C[i-1]
    N[i] = ifelse(sum(y==t_idx[i] & d==1)>0, sum(y==t_idx[i] & d==1), 0)
    C[i] = ifelse(sum(y==t_idx[i] & d==0)>0, sum(y==t_idx[i] & d==0), 0)
    
    if(Y[i]<0){Y[i] = 0}
    
    S[i] = ifelse(Y[i]==0, S[i-1], S[i-1]*(1-(N[i]/Y[i])))
    H[i] = ifelse(Y[i]*(Y[i]-N[i])==0, 0, N[i]/(Y[i]*(Y[i]-N[i])))
    
    if(S[i]<0){S[i] = 0}
    
    D[i] = sum(H[2:i])
    E[i] = sqrt((S[i]**2)*D[i])
    
    if(is.na(S[i])){S[i] = 0}
    if(is.na(E[i])){E[i] = 0}
  }
  
  #--output--
  out           = as.data.frame(cbind(t_idx, Y, N, C, S, E))
  colnames(out) = c("t_idx", "n_risk", "n_event", "n_censor", "surv", "se")
  
  #--to match the output of survfit--
  out2 = out[t_idx!=0,]
  
  #--
  Z2 = list()
  Z2$out      = out2
  Z2$t_idx    = out2[,"t_idx"]
  Z2$n_risk   = out2[,"n_risk"]
  Z2$n_event  = out2[,"n_event"]
  Z2$n_censor = out2[,"n_censor"]
  Z2$surv     = out2[,"surv"]
  Z2$se       = out2[,"se"]
  
  return(Z2)
}



my.rmst2reg=function(y, delta, arm, x, tau, w=rep(1,length(y))){
  
  n=length(y)
  x=as.matrix(cbind(1, x))
  p=length(x[1,])
  
  y0=pmin(y, tau)
  d0=delta
  d0[y0==tau]=1
  
  d10=d0[arm==1]
  d00=d0[arm==0]
  y10=y0[arm==1]
  y00=y0[arm==0]
  x1=x[arm==1,]
  x0=x[arm==0,]
  n1=length(d10)
  n0=length(d00)
  
  id1=order(y10)
  y10=y10[id1]
  d10=d10[id1]
  x1=x1[id1,]
  
  id0=order(y00)
  y00=y00[id0]
  d00=d00[id0]
  x0=x0[id0,]
  
  fitc1=my.func_surv(y10, 1-d10)
  fitc0=my.func_surv(y00, 1-d00)
  
  weights1=d10/rep(pmax(fitc1$surv,0.001), table(y10))
  weights0=d00/rep(pmax(fitc0$surv,0.001), table(y00))
  
  w1=w[arm==1]
  w0=w[arm==0]
  w1=w1[id1]
  w0=w0[id0]
  weights=c(weights1, weights0)*c(w1,w0)
  
  
  fitt=lm(c(y10,y00)~ rbind(x1, x0)-1, weights=weights)
  
  return(fitt)
}




fit.rmst.reg <- function(dat0,dat1,tau){
  tau.temp <- min(max(dat0$time), max(dat1$time))
  tau.temp <- min(tau.temp, tau)
  
  dat0$A=0
  dat1$A=1
  dat <- rbind(dat0,dat1)
  cov <- data.frame(A=dat$A, X=dat$X, AX=dat$A*dat$X)
  rmst_fit <- my.rmst2reg(y = dat$time,
                          delta = dat$status,
                          x = cov,
                          arm = dat$A,
                          tau = tau.temp)
  return(rmst_fit)
}


find.cutpoint <- function(rmst_fit){
  
  gamma <- coef(rmst_fit)
  cutpoint <- -gamma[2]/gamma[4]
  cutpoint <- min(max(cutpoint,0),1)
  sign <- gamma[4]
  
  
  if ((cutpoint == 0 & sign < 0) | (cutpoint == 1 & sign > 0)){
    direction <- "all negative"
  }else if (sign > 0){
    direction <- "right"
  }else if (sign < 0){
    direction <- "left"
  }
  
  
  
  return(c(cutpoint, direction))
}



stat1 <- function(dat0,dat1,tau){
  time_pos=c(dat0$time, dat1$time)
  status_pos=c(dat0$status, dat1$status)
  arm_pos=c(dat0$A, dat1$A)
  
  fit_pos=rmst2(time_pos, status_pos, arm_pos, tau=tau)
  rmst1est_pos=fit_pos$RMST.arm1$rmst[1:2]
  rmst0est_pos=fit_pos$RMST.arm0$rmst[1:2]
  d_pos=(rmst1est_pos[1]-rmst0est_pos[1])
  se_pos=sqrt(rmst1est_pos[2]^2+rmst0est_pos[2]^2)
  z_pos=d_pos/se_pos
  
  return(as.numeric(c(z_pos, d_pos, se_pos)))
}








