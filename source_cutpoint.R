#################################################################################
########################## RMST Regression Method ###############################

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



#################################################################################
########################## Prediction Method  ###############################

# alternative setting #

true.cut.piece <- function(tau,lambda0,lambda1,
                           beta0,beta1,t.change0,t.change1){
  mu <- function(x,t,t.change,lambda,beta){
    if (t <= t.change[1]){
      integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change[1] & t <= t.change[2]){
      integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                lower = 0, upper = t.change[1], rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*t.change[1] -
                                    lambda[2]*exp(beta*(1-x))*(y-t.change[1]))},
                  lower = t.change[1], upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change[2]){
      integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                lower = 0, upper = t.change[1], rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*t.change[1] -
                                    lambda[2]*exp(beta*(1-x))*(y-t.change[1]))},
                  lower = t.change[1], upper = t.change[2], rel.tol = 1e-10)$value +
        integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*t.change[1] -
                                    lambda[2]*exp(beta*(1-x))*(t.change[2]-t.change[1]) - 
                                    lambda[3]*exp(beta*(1-x))*(y-t.change[2]))},
                  lower = t.change[2], upper = t, rel.tol = 1e-10)$value
    }
  }
  
  if (beta1 == beta0) {
    cut_t <- NA
    direction <- "all negative"
  }else if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))*
            (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) < 0){
    cut_t <- uniroot(function(x){mu(x,tau,t.change1,lambda1,beta1) - mu(x,tau,t.change0,lambda0,beta0)},
                     interval = c(0,1), tol = 1e-10)$root
    if (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0) > 0){
      direction <- "left"}else {
        direction <- "right"}
  }else if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0)) >= 0 & 
            (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) > 
            (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))){
    cut_t <- 0
    direction <- "right"
  }else if ((mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) >= 0 & 
            (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) <
            (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))){
    cut_t <- 1
    direction <- "left"
  }else if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0)) < 0 & 
            (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) < 0){
    cut_t <- NA
    direction <- "all negative"
  }
  
  
  
  return(c(cut_t, direction))
}




fit.piece.0 <- function(dat,k,t.change,t.analysis){
  if (k == 0){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ 1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }else if (k >= 1){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.change,t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ -1 + interval + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }
  return(fit)
}


fit.piece.1 <- function(dat,k,t.change,t.analysis){
  dat$X1 <- 1 - dat$X
  if (k == 0){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ X1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }else if (k >= 1){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.change,t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ -1 + interval + X1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }
  return(fit)
}



cut.pred1 <- function(dat0, dat1, tau, t.change, t.analysis){
  # when the number of cutpoint and location are known
  fit0 <- fit.piece.0(dat0,0,NA,t.analysis)
  fit1 <- fit.piece.1(dat1,1,t.change,t.analysis)
  
  beta1 <- ifelse(coef(fit1)[3]/sqrt(vcov(fit1)[3,3]) <= qnorm(0.975), 0, coef(fit1)[3])
  cutpoint0 <- true.cut.piece(tau,
                              lambda0 = rep(exp(coef(fit0)[1]),3),
                              lambda1 = exp(coef(fit1)[c(1,2,2)]),
                              beta0 = 0,
                              beta1 = beta1,
                              t.change0 = c(tau,tau),
                              t.change1 = c(t.change,tau))
  return(list(cutpoint0 = cutpoint0, beta1 = beta1))
}




cut.pred2 <- function(dat0, dat1, tau, t.analysis){
  # when the number of cutpoint is known, but location is unknown
  fit_piece1 <- piecewiseExp_MLE(dat1$time, dat1$status, 1)
  t.change <- fit_piece1$tau
  
  fit0 <- fit.piece.0(dat0,0,NA,t.analysis)
  fit1 <- fit.piece.1(dat1,1,t.change,t.analysis)
  
  beta1 <- ifelse(coef(fit1)[3]/sqrt(vcov(fit1)[3,3]) <= qnorm(0.975), 0, coef(fit1)[3])
  cutpoint0 <- true.cut.piece(tau,
                              lambda0 = rep(exp(coef(fit0)[1]),3),
                              lambda1 = exp(coef(fit1)[c(1,2,2)]),
                              beta0 = 0,
                              beta1 = beta1,
                              t.change0 = c(tau,tau),
                              t.change1 = c(t.change,tau))
  return(list(cutpoint0 = cutpoint0, 
              t.change = t.change, beta1 = beta1))
}





cut.pred3 <- function(dat0, dat1, tau, t.analysis){
  # when the number of cutpoint and locations are both unknown
  rej0 <- 1
  k0 <- 0
  while (rej0 == 1 & k0 < 3){
    k0 <- k0 + 1
    fit_piece0 <- piecewiseExp_MLE(dat0$time, dat0$status, k0)
    test0 <- piecewiseExp_test_changepoint(fit_piece0)
    rej0 <- prod(test0$reject)
  }
  k0 <- k0 - 1
  fit_piece0 <- piecewiseExp_MLE(dat0$time, dat0$status, k0)
  t.change0 <- fit_piece0$tau
  fit0 <- fit.piece.0(dat0,k0,t.change0,t.analysis)
  
  if (k0 == 0){
    lambda0 <- rep(exp(coef(fit0)[1]),3)
    beta0 <- 0
    t.change0 <- c(tau,tau)
  }else if (k0 == 1){
    lambda0 <- exp(coef(fit0)[c(1,2,2)])
    beta0 <- 0
    t.change0 <- c(t.change0, tau)
  }else if (k0 == 2){
    lambda0 <- exp(coef(fit0)[1:3])
    beta0 <- 0
    t.change0 <- pmin(t.change0, tau)
  }
  
  
  rej1 <- 1
  k1 <- 0
  while (rej1 == 1 & k1 < 3){
    k1 <- k1 + 1
    fit_piece1 <- piecewiseExp_MLE(dat1$time, dat1$status, k1)
    test1 <- piecewiseExp_test_changepoint(fit_piece1)
    rej1 <- prod(test1$reject)
  }
  k1 <- k1 - 1
  fit_piece1 <- piecewiseExp_MLE(dat1$time, dat1$status, k1)
  t.change1 <- fit_piece1$tau
  fit1 <- fit.piece.1(dat1,k1,t.change1,t.analysis)
  
  if (k1 == 0){
    lambda1 <- rep(exp(coef(fit1)[1]),3)
    beta1 <- ifelse(coef(fit1)[2]/sqrt(vcov(fit1)[2,2]) <= qnorm(0.975), 0, coef(fit1)[2])
    t.change1 <- c(tau,tau)
  }else if (k1 == 1){
    lambda1 <- exp(coef(fit1)[c(1,2,2)])
    beta1 <- ifelse(coef(fit1)[3]/sqrt(vcov(fit1)[3,3]) <= qnorm(0.975), 0, coef(fit1)[3])
    t.change1 <- c(t.change1, tau)
  }else if (k1 == 2){
    lambda1 <- exp(coef(fit1)[1:3])
    beta1 <- ifelse(coef(fit1)[4]/sqrt(vcov(fit1)[4,4]) <= qnorm(0.975), 0, coef(fit1)[4])
    t.change1 <- pmin(t.change1, tau)
  }
  
  
  
  cutpoint0 <- true.cut.piece(tau,
                              lambda0, lambda1,
                              beta0, beta1,
                              t.change0, t.change1)
  
  
  return(list(cutpoint0 = cutpoint0, 
              t.change0 = t.change0, 
              t.change1 = t.change1, 
              k0 = k0, k1 = k1, beta1 = beta1))
}







# Null Setting #


true.cut.piece.null <- function(tau,lambda0,lambda1,
                                beta0,beta1,t.change0,t.change1){
  mu <- function(x,t,t.change,lambda,beta){
    if (t <= t.change){
      integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*t.change -
                                    lambda[2]*exp(beta*(1-x))*(y-t.change))},
                  lower = t.change, upper = t, rel.tol = 1e-10)$value
    }
  }
  
  
  if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))*
      (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) < 0){
    cut_t <- uniroot(function(x){mu(x,tau,t.change1,lambda1,beta1) - mu(x,tau,t.change0,lambda0,beta0)},
                     interval = c(0,1), tol = 1e-10)$root
    if (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0) > 0){
      direction <- "left"}else {
        direction <- "right"}
  }else if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0)) >= 0 & 
            (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) >= 
            (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))){
    cut_t <- 0
    direction <- "right"
  }else if ((mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) >= 0 & 
            (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) <= 
            (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))){
    cut_t <- 1
    direction <- "left"
  }else if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0)) < 0 & 
            (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) < 0){
    cut_t <- NA
    direction <- "all negative"
  }
  
  
  
  return(c(cut_t, direction))
}




fit.piece.null.1 <- function(dat,k,t.change,t.analysis){
  dat$X1 <- 1 - dat$X
  if (k == 0){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ X1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }else if (k >= 1){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.change,t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ -1 + interval + X1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }
  return(fit)
}


fit.piece.null.0 <- function(dat,k,t.change,t.analysis){
  if (k == 0){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ 1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }else if (k >= 1){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.change,t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ -1 + interval + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }
  return(fit)
}







# prediction
cut.pred1.null <- function(dat0, dat1, tau, t.analysis){
  # when the number of cutpoint and location are known
  fit0 <- fit.piece.null.0(dat0,0,tau,t.analysis)
  fit1 <- fit.piece.null.1(dat1,0,tau,t.analysis)
  #fit0 <- fit.piece.null(dat0,0,tau,t.analysis)
  #fit1 <- fit.piece.null(dat1,0,tau,t.analysis)
  
  beta1 <- ifelse(coef(fit1)[2]/sqrt(vcov(fit1)[2,2]) <= qnorm(0.975),
                  0, coef(fit1)[2])
  
  cutpoint0 <- true.cut.piece.null(tau,
                                   lambda0 = rep(exp(coef(fit0)[1]),2),
                                   lambda1 = rep(exp(coef(fit1)[1]),2),
                                   beta0 = 0,
                                   beta1 = beta1,
                                   t.change0 = tau,
                                   t.change1 = tau)
  return(list(cutpoint0 = cutpoint0, beta1 = beta1))
}







cut.pred3.null <- function(dat0, dat1, tau, t.analysis){
  # when the number of cutpoint is known, but location is unknown
  rej0 <- 1
  k0 <- 0
  while (rej0 == 1 & k0 < 2){
    k0 <- k0 + 1
    fit_piece0 <- piecewiseExp_MLE(dat0$time, dat0$status, k0)
    test0 <- piecewiseExp_test_changepoint(fit_piece0)
    rej0 <- prod(test0$reject)
  }
  k0 <- k0 - 1
  fit_piece0 <- piecewiseExp_MLE(dat0$time, dat0$status, k0)
  t.change0 <- fit_piece0$tau
  fit0 <- fit.piece.null.0(dat0,k0,t.change0,t.analysis)
  
  if (k0 == 0){
    lambda0 <- rep(exp(coef(fit0)[1]),2)
    beta0 <- 0
    t.change0 <- tau
  }else if (k0 == 1){
    lambda0 <- exp(coef(fit0)[1:2])
    beta0 <- 0
    t.change0 <- pmin(t.change0, tau)
  }
  
  
  rej1 <- 1
  k1 <- 0
  while (rej1 == 1 & k1 < 2){
    k1 <- k1 + 1
    fit_piece1 <- piecewiseExp_MLE(dat1$time, dat1$status, k1)
    test1 <- piecewiseExp_test_changepoint(fit_piece1)
    rej1 <- prod(test1$reject)
  }
  k1 <- k1 - 1
  fit_piece1 <- piecewiseExp_MLE(dat1$time, dat1$status, k1)
  t.change1 <- fit_piece1$tau
  fit1 <- fit.piece.null.1(dat1,k1,t.change1,t.analysis)
  
  if (k1 == 0){
    lambda1 <- rep(exp(coef(fit1)[1]),2)
    beta1 <- ifelse(coef(fit1)[2]/sqrt(vcov(fit1)[2,2]) <= qnorm(0.975),
                    0, coef(fit1)[2])
    t.change1 <- tau
  }else if (k1 == 1){
    lambda1 <- exp(coef(fit1)[1:2])
    beta1 <- ifelse(coef(fit1)[3]/sqrt(vcov(fit1)[3,3]) <= qnorm(0.975),
                    0, coef(fit1)[3])
    t.change1 <- pmin(t.change1, tau)
  }
  
  
  
  cutpoint0 <- true.cut.piece.null(tau,
                                   lambda0, lambda1,
                                   beta0, beta1,
                                   t.change0, t.change1)
  
  
  return(list(cutpoint0 = cutpoint0, 
              t.change0 = t.change0, 
              t.change1 = t.change1, 
              k0 = k0, k1 = k1,
              beta1 = beta1))
}














