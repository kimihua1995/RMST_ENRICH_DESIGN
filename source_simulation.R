#################################################################################
############# simulate piecewise exponential time-to-event data #################
#################################################################################
sim.data.piece <- function(X,lambda,beta,t.change,t.enroll.start,t.enroll.end,t.IA){
  # administrative censoring only
  n <- length(X)
  u <- runif(n,0.00001,1)
  time.event <- -log(u)*exp(-(1-X)*beta)/lambda[1]
  pos <- time.event > t.change
  #time.event[pos] <- (-log(u[pos])-lambda[1]*exp(beta*X[pos])*t.change+lambda[2]*exp(beta*X[pos])*t.change)*
  #  exp(-X[pos]*beta)/lambda[2]
  time.event[pos] <- -log(u[pos])*exp(-(1-X[pos])*beta)/lambda[2] + (1 - lambda[1]/lambda[2])*t.change
  time.enroll <- runif(n, t.enroll.start, t.enroll.end)
  time.loss <- rexp(n)/0.12
  #time.loss <- rep(10000,n)
  time.censor <- pmin(time.loss, t.IA - time.enroll)
  time <- pmin(time.event, time.censor)
  status <- as.numeric(time.event < time.censor)
  
  dat <- data.frame(X = X, time.event = time.event, time.enroll = time.enroll,
                    time.loss = time.loss, time.censor = time.censor, time = time, status = status)
  return(dat)
}



#################################################################################
###### True biomarker cutpoint given piecewise exponential hazard model #########
#################################################################################
true.cut.piece0 <- function(tau,lambda0,lambda1,beta0,beta1,t.change){
  mu0 <- function(x,t){
    if (t <= t.change){
      integrate(function(y){exp(-lambda0[1]*exp(beta0*(1-x))*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda0[1]*exp(beta0*(1-x))*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda0[1]*exp(beta0*(1-x))*t.change-lambda0[2]*exp(beta0*(1-x))*(y-t.change))},
                  lower = t.change, upper = t, rel.tol = 1e-10)$value
    }
  }
  mu1 <- function(x,t){
    if (t <= t.change){
      integrate(function(y){exp(-lambda1[1]*exp(beta1*(1-x))*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda1[1]*exp(beta1*(1-x))*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda1[1]*exp(beta1*(1-x))*t.change-lambda1[2]*exp(beta1*(1-x))*(y-t.change))},
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




#################################################################################
### True marginal RMST within certain biomarker subgroup ####
#################################################################################
rmst.true <- function(lambda, beta, t.change, tau, lower, upper){
  mu <- function(x, tau){
    ifelse(tau <= t.change,
           integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                     lower = 0, upper = tau, rel.tol = 1e-10)$value,
           integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                     lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
             integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*t.change-lambda[2]*exp(beta*(1-x))*(y-t.change))},
                       lower = t.change, upper = tau, rel.tol = 1e-10)$value)
  }
  
  integrate(Vectorize(function(x){mu(x,tau)/(upper-lower)}), lower = lower, upper = upper, rel.tol = 1e-10)$value
}





#################################################################################
############### function for modifying censoring time ###########################
#################################################################################
dat.modify <- function(dat,t.IA){
  dat$time.censor <- pmin(dat$time.loss, t.IA - dat$time.enroll)
  dat$time <- pmin(dat$time.event, dat$time.censor)
  dat$status <- as.numeric(dat$time.event < dat$time.censor)
  return(dat)
}





#################################################################################
######## function for simulation (alternative setting) ####################
#################################################################################
my.design <- function(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, alpha1, alpha2, true.cutoff, design = 1, pred.method = 1){
  
  stat.FA1 <- data.frame(z_pos = rep(NA,B), d_pos = rep(NA,B), se_pos = rep(NA,B))
  stat.FA1.true <- data.frame(z_pos = rep(NA,B), d_pos = rep(NA,B), se_pos = rep(NA,B))
  
  cutpoint.IA.list <- cutpoint.FA.list <- data.frame(cutpoint = rep(NA,B), direction = rep(NA,B))
  
  change.point <- rep(NA,B)
  change.point0 <- change.point1 <- matrix(NA, B, 2)
  k0.list <- k1.list <- rep(NA,B)
  
  rej1 <- rej1.true <- rep(FALSE,B)
  gamma3.FA.list <- beta1.list <- rep(NA,B)
  
  N.neg <- rep(NA,B)
  
  set.seed(seed)
  for (b in 1:B){
    if (!b%%1000) {print(b)}
    # Stage I
    X0_1 <- runif(n1, 0, 1)
    X1_1 <- runif(n1, 0, 1)
    dat0_1 <- sim.data.piece(X0_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    dat1_1 <- sim.data.piece(X1_1,lambda1,beta1,t.change,0,t.enroll1,t.enrich)
    dat0_1$S = 1
    dat1_1$S = 1
    
    
    if (design == 1){
      X0_2 <- runif(n2, 0, 1)
      X1_2 <- runif(n2, 0, 1)
    }else if (design == 2){
      if (pred.method == 1){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred1(dat0_1, dat1_1, tau, t.change, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {next}
      }else if (pred.method == 2){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred2(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {
          next
        }else {
          change.point[b] <- cut_pred$t.change
        }
      }else if (pred.method == 3){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred3(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {
          next
        }else {
          change.point0[b,] <- cut_pred$t.change0
          change.point1[b,] <- cut_pred$t.change1
          k0.list[b] <- cut_pred$k0
          k1.list[b] <- cut_pred$k1
        }
      }
      
      cutpoint.IA <- cut_pred$cutpoint0
      cutpoint.IA.list[b,] <- cutpoint.IA
      beta1.list[b] <- cut_pred$beta1
      
      
      # Stage II
      if (cutpoint.IA[2] == "all negative") {
        #cutpoint.IA.list[b,] <- c(NA, "all negative")
        X0_2 <- runif(n2, 0, 1)
        X1_2 <- runif(n2, 0, 1)
      }else if (cutpoint.IA[2] == "left"){
        X0_2 <- runif(n2, 0, as.numeric(cutpoint.IA[1]))
        X1_2 <- runif(n2, 0, as.numeric(cutpoint.IA[1]))
      }else if (cutpoint.IA[2] == "right" & as.numeric(cutpoint.IA[1]) > 0.9){
        X0_2 <- runif(n2, 0, 1)
        X1_2 <- runif(n2, 0, 1)
      }else if (cutpoint.IA[2] == "right" & as.numeric(cutpoint.IA[1]) < 0.9){
        X0_2 <- runif(n2, as.numeric(cutpoint.IA[1]), 1)
        X1_2 <- runif(n2, as.numeric(cutpoint.IA[1]), 1)
      }
      
    }
    
    
    dat0_2 <- sim.data.piece(X0_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat1_2 <- sim.data.piece(X1_2,lambda1,beta1,t.change,t.enrich,t.enroll2,t.FA)
    dat0_2$S = 2
    dat1_2$S = 2
    dat0 <- rbind(dat0_1, dat0_2); dat0$A = 0
    dat1 <- rbind(dat1_1, dat1_2); dat1$A = 1
    
    
    
    #####################################
    ## FA
    dat0 <- dat.modify(dat0, t.FA)
    dat1 <- dat.modify(dat1, t.FA)
    dat <- rbind(dat0, dat1)
    
    ## estimate the cutoff point at FA
    rmst_fit_FA <- fit.rmst.reg(dat0, dat1, tau)
    gamma3_FA <- coef(rmst_fit_FA)[4]/sqrt(vcov(rmst_fit_FA)[4,4])
    gamma3.FA.list[b] <- gamma3_FA
    cutpoint.FA <- find.cutpoint(rmst_fit_FA)
    cutpoint.FA.list[b,] <- cutpoint.FA
    
    
    if (cutpoint.FA[2] == "all negative"){
      cutpoint.FA.list[b,] <- c(NA, NA, "all negative")
    }else if (gamma3_FA <= qnorm(1 - alpha1)){
      stat.FA1[b,] <- stat1(dat0,dat1,0,"right",tau)
      rej1[b] <- stat.FA1[b,1] > qnorm(1 - alpha2)
    }else{
      stat.FA1[b,] <- stat1(dat0,dat1,cutpoint.FA[1],cutpoint.FA[2],tau)
      rej1[b] <- stat.FA1[b,1] > qnorm(1 - alpha2)
      
    }
    
    
    N.neg[b] <- sum(dat$X < true.cutoff)
    
    
    # true.cutoff
    stat.FA1.true[b,] <- stat1(dat0,dat1,true.cutoff,"right",tau)
    rej1.true[b] <- stat.FA1.true[b,1] > qnorm(1 - alpha2)
    
    
    
  }
  
  return(list(stat.FA1 = stat.FA1, stat.FA1.true = stat.FA1.true, 
              cutpoint.IA = cutpoint.IA.list, cutpoint.FA = cutpoint.FA.list,
              rej1 = rej1, rej1.true = rej1.true,
              gamma3.FA = gamma3.FA.list,
              N.neg = N.neg, beta1 = beta1.list,
              change.point = change.point,
              change.point0 = change.point0, change.point1 = change.point1,
              k0.list = k0.list, k1.list = k1.list))
  
}


#####################################################################
######## function for simulation (null setting) #####################
#####################################################################
my.design.null <- function(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                           t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, alpha1, alpha2, pred.method=1, design=1){
  
  cutpoint.IA.list <- cutpoint.FA.list <- data.frame(cutpoint = rep(NA,B), direction = rep(NA,B))
  
  change.point <- rep(NA,B)
  change.point0 <- change.point1 <- matrix(NA, B, 2)
  k0.list <- k1.list <- rep(NA,B)
  
  a1 <- length(alpha1); a2 <- length(alpha2)
  rej1 <- matrix(FALSE,B,a1*a2)
  gamma3.FA.list <- beta1.list <- rep(NA,B)
  
  set.seed(seed)
  for (b in 1:B){
    if (!b%%1000) {print(b)}
    #set.seed(b)
    # Stage I
    X0_1 <- runif(n1, 0, 1)
    X1_1 <- runif(n1, 0, 1)
    dat0_1 <- sim.data.piece(X0_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    dat1_1 <- sim.data.piece(X1_1,lambda1,beta1,t.change,0,t.enroll1,t.enrich)
    dat0_1$S = 1
    dat1_1$S = 1
    
    if (design == 1){
      X0_2 <- runif(n2, 0, 1)
      X1_2 <- runif(n2, 0, 1)
    }else if (design == 2){
      if (pred.method == 1){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred1.null(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {next}
      }else if (pred.method == 3){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred3.null(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {
          next
        }else {
          change.point0[b,] <- cut_pred$t.change0
          change.point1[b,] <- cut_pred$t.change1
          k0.list[b] <- cut_pred$k0
          k1.list[b] <- cut_pred$k1
        }
      }
      
      cutpoint.IA <- cut_pred$cutpoint0
      cutpoint.IA.list[b,] <- cutpoint.IA
      beta1.list[b] <- cut_pred$beta1
      
      # Stage II
      if (cutpoint.IA[2] == "all negative") {
        #cutpoint.IA.list[b,] <- c(NA, "all negative")
        X0_2 <- runif(n2, 0, 1)
        X1_2 <- runif(n2, 0, 1)
      }else if (cutpoint.IA[2] == "left"){
        X0_2 <- runif(n2, 0, as.numeric(cutpoint.IA[1]))
        X1_2 <- runif(n2, 0, as.numeric(cutpoint.IA[1]))
      }else if (cutpoint.IA[2] == "right" & as.numeric(cutpoint.IA[1]) > 0.9){
        X0_2 <- runif(n2, 0, 1)
        X1_2 <- runif(n2, 0, 1)
      }else if (cutpoint.IA[2] == "right" & as.numeric(cutpoint.IA[1]) < 0.9){
        X0_2 <- runif(n2, as.numeric(cutpoint.IA[1]), 1)
        X1_2 <- runif(n2, as.numeric(cutpoint.IA[1]), 1)
      }
      
    }
    
    
    
    dat0_2 <- sim.data.piece(X0_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat1_2 <- sim.data.piece(X1_2,lambda1,beta1,t.change,t.enrich,t.enroll2,t.FA)
    dat0_2$S = 2
    dat1_2$S = 2
    dat0 <- rbind(dat0_1, dat0_2); dat0$A = 0
    dat1 <- rbind(dat1_1, dat1_2); dat1$A = 1
    
    
    
    #####################################
    ## FA
    dat0 <- dat.modify(dat0, t.FA)
    dat1 <- dat.modify(dat1, t.FA)
    dat <- rbind(dat0, dat1)
    
    ## estimate the cutoff point at FA
    rmst_fit_FA <- fit.rmst.reg(dat0, dat1, tau)
    gamma3_FA <- coef(rmst_fit_FA)[4]/sqrt(vcov(rmst_fit_FA)[4,4])
    gamma3.FA.list[b] <- gamma3_FA
    cutpoint.FA <- find.cutpoint(rmst_fit_FA)
    cutpoint.FA.list[b,] <- cutpoint.FA
    
    
    
    
    #######################################
    stat.FA1.all <- stat1(dat0,dat1,0,"right",tau)
    
    for (i in 1:a1){
      if (cutpoint.FA[2] == "all negative"){
        cutpoint.FA.list[b,] <- c(NA, "all negative")
      }else if (gamma3_FA <= qnorm(1 - alpha1[i])){
        for (j in 1:a2){
          rej1[b,a2*(i-1)+j] <- stat.FA1.all[1] > qnorm(1 - alpha2[j])
        }
      }else{
        stat.FA1 <- stat1(dat0,dat1,cutpoint.FA[1],cutpoint.FA[2],tau)
        for (j in 1:a2){
          rej1[b,a2*(i-1)+j] <- stat.FA1[1] > qnorm(1 - alpha2[j])
        }
      }
    }
    
    
    
    
    
  }
  
  return(list(cutpoint.IA = cutpoint.IA.list, cutpoint.FA = cutpoint.FA.list,
              rej1 = rej1,
              gamma3.FA = gamma3.FA.list,
              beta1 = beta1.list,
              change.point = change.point,
              change.point0 = change.point0, change.point1 = change.point1,
              k0.list = k0.list, k1.list = k1.list))
  
}






#################################################################################
#################### RMST estimator in positive group ###########################
#################################################################################
stat1 <- function(dat0,dat1,cutpoint,pos,tau){
  if (pos == "right"){
    dat0_pos <- dat0[dat0$X >= cutpoint,]
    dat1_pos <- dat1[dat1$X >= cutpoint,]
  }else if (pos == "left"){
    dat0_pos <- dat0[dat0$X <= cutpoint,]
    dat1_pos <- dat1[dat1$X <= cutpoint,]
  }
  
  time_pos=c(dat0_pos$time, dat1_pos$time)
  status_pos=c(dat0_pos$status, dat1_pos$status)
  arm_pos=c(rep(0,nrow(dat0_pos)), rep(1,nrow(dat1_pos)))
  
  if (min(max(dat0_pos$time), max(dat1_pos$time)) < tau) {
    tau.temp <- min(max(dat0_pos$time), max(dat1_pos$time))
  }else {tau.temp <- tau}
  
  fit_pos=rmst2(time_pos, status_pos, arm_pos, tau=tau.temp)
  rmst1est_pos=fit_pos$RMST.arm1$rmst[1:2]
  rmst0est_pos=fit_pos$RMST.arm0$rmst[1:2]
  d_pos=(rmst1est_pos[1]-rmst0est_pos[1])
  se_pos=sqrt(rmst1est_pos[2]^2+rmst0est_pos[2]^2)
  z_pos=d_pos/se_pos
  
  
  
  return(as.numeric(c(z_pos,d_pos,se_pos)))
}







