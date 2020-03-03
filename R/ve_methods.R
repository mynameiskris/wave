### VE funtions

# Estimate VE using method from Durham et al. 1988 (AJE)
durham_ve <- function(x, df = 2, nsmo = 40,var,...){
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
  temp <- c(pred.x, xx)
  lmat <- ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:nsmo, ]
  xmat <- lmat[-(1:nsmo), ]
  qmat <- qr(xmat)
  if (x$transform!="identity") 
    stop("please re-fit the Cox model with the identity transform")
  if (qmat$rank < df) 
    stop("Spline fit is singular, try a smaller degrees of freedom")
  # se
  bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
  xtx <- bk %*% t(bk)
  seval <- ((pmat %*% xtx) * pmat) %*% rep(1, df)
  
  # ve estimate  
  yhat.beta <- pmat %*% qr.coef(qmat, yy)
  yhat <- 1-exp(yhat.beta) # VE estimate
    
  # se
  temp <- 2 * sqrt(x$var[1,1] * seval)
  ylow.beta <- yhat.beta - temp
  yup.beta <- yhat.beta + temp
  ylow <- 1-exp(yup.beta)
  yup <- 1-exp(ylow.beta)
  yr <- range(yhat, yup, ylow)
  
  # output  
  rtn <- tibble(time = round(pred.x), 
                ve = as.numeric(yhat),
                ve_lower = as.numeric(ylow),
                ve_upper = as.numeric(yup))
  return(rtn)
}

# Estimate VE using method from Ferdinands et al. 2017 (CID?)
ferdinands_ve <- function(dat, splines = FALSE){
  
  if (splines == FALSE){ # fit logistic regression model with dichotomous vaccination variable
    ve.logr1<- glm(FARI ~ V + DINF, family=binomial("logit"), data = dat)
    exp(cbind(OR = coef(ve.logr1), confint(ve.logr1)))
  }
  else{ # fit splines
    vacc.dat <- dat %>% filter(V == 1)
    qus <- quantile(vacc.dat$DINF, probs = c(0.10,0.98,0.99,0.80,0.20,0.50))
    names(qus) <- c('qu0','qu1','qu2','qu3','qu4','qu5')
  
    # Spline model - days from vaccination to onset modeled as spline
    # Specification of spline knots based on comparison of model fit
    # for numerous alternative specifications (not shown here)
  
    ve.logr2<- glm(FARI ~ V +  ns(DINF, knots=c(qus["qu1"], qus["qu2"])), family=binomial("logit"), data = dat)
    #exp(cbind(OR = coef(ve.logr2), confint(ve.logr2)))
    summary(ve.logr2)
  
    aOR.avg <- exp(ve.logr1$coefficient[5])
    aVE.avg <- (1 - aOR.avg)*100; aVE.avg
  
    LH1 <- logLik(ve.logr1)
    LH2 <- logLik(ve.logr2)
  
    delta <- 2*(LH2[1] - LH1[1])
  }
  rtn <- list(lh1 = LH1, lh2 = LH2, delta = delta)
  return(rtn)
  
}

# Estimate VE using method from Tian et al. 2005 
tian_ve <- function(dat, n_timepoint_breaks = 40, alpha = 0.05){
  reject_h0 <- 0
  # fit timecox model
  fit = timecox(Surv(DINF_new, FARI) ~ V, data = dat, n.sim = 500, max.time = 700)
  KS_pvalue = fit$pval.testBeqC[2]
  if (KS_pvalue < alpha){reject_h0 = 1} 
  # cummulative estimates
  # rtn_cum <- tibble(alpha0 = fit$cum[,2], beta_t = fit$cum[,-(1:2)]) %>% 
  #               mutate(ve = 1 - exp(beta_t))
  # de-cummulative estimates
    predicted.timepts <- seq(from = min(fit$cum[,1]), to = max(fit$cum[,1]), length = n_timepoint_breaks)
    decumulated.ests = CsmoothB(fit$cum, predicted.timepts, b = 1)
    timecox.hazardrate = exp(decumulated.ests[,2])
    timecox.var = decumulated.ests[,-(1:2)]
    rtn <- tibble(time = predicted.timepts,hazardrate = timecox.hazardrate) %>% 
              mutate(ve = 1 - hazardrate)
  # output
  return(list(output = rtn, reject_h0 = reject_h0))
}

# Estimate VE using method from Ainslie et al. 2017 (SIM)
 ainslie_ve <- function(dat, n_days){
   # calculate prevalence
   N <- length(unique(dat$ID))
   prev <- numeric(n_days)
   for (d in 1:n_days){
     prev[d] <- length(which(dat$DINF == d))/N
   }
   prev <- ifelse(prev == 0, 0.001, prev)
   x <- list(n_days = n_days, prev = prev, dinf = dat$DINF_new, v = dat$V)
   
   logLik <- function(x, par){
     gamma_0 = par[1]
     gamma_1 = par[2]
     
    # conditional probabilities
       pi_00 <- 1 - gamma_0 * x$prev
       pi_01 <- 1 - gamma_1 * x$prev
       pi_10 <- gamma_0 * x$prev
       pi_11 <- gamma_1 * x$prev
    # unconditional probabilities
       phi_00 <- phi_01 <- phi_10 <- phi_11 <- c(1,rep(0,x$n_days-1))
       for (d in 2:x$n_days){
        phi_00[d] <- pi_00[d] * phi_00[d-1]
        phi_01[d] <- pi_01[d] * phi_01[d-1]
        phi_10[d] <- pi_10[d] * phi_00[d-1]
        phi_11[d] <- pi_11[d] * phi_01[d-1]
       }
    # personal contribution to the likelihood
      Li <- numeric(N)
      for (i in 1:N){
        if (x$dinf[i] < 999){
          if (x$v[i] == 0){Li[i] <- phi_10[x$dinf[i]]
          } else {Li[i] <- phi_11[x$dinf[i]]}
        } else if (x$dinf[i] == 999){
          if (x$v[i] == 0){Li[i] <- phi_00[x$n_days]
          } else {Li[i] <- phi_01[x$n_days]}
        }
      }

    # likelihood 
      L = prod(Li)
      
      return(-sum(log(L)))
   }
   # maximum likelihood estimates
    mle <- optim(par = c(0.1, 0.1), logLik, x = x, method = "L-BFGS-B", lower = 0, upper = 1, hessian = TRUE)
   
   # calculate lamdas for each individual
    lambda <- numeric(N)
    for (i in 1:N){
      if (x$v[i] == 0) {lambda[i] <- sum(mle$par[1]*prev*c(phi_00[-1],phi_00[1]))
      } else {lambda[i] <- sum(mle$par[2]*prev*c(phi_01[-1],phi_01[1]))}
    }
    ve = 1 - (sum(lambda[which(x$v == 1)])/sum(lambda[which(x$v == 0)]))
    return(ve)
 }

