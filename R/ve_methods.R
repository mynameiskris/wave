### VE funtions

# Estimate VE using method from Durham et al. 1988 (AJE)
durham_ve <- function(x, df = 2, n_days, n_periods, n_days_period, var,...){
  xx <- x$x
  yy <- x$y
  d <- nrow(yy)
  df <- max(df)
  nvar <- ncol(yy)
  if (length(n_days) == 1) {
    pred.x <- seq(from = min(xx), to = max(xx), length = n_days)
  } else {pred.x <- n_days}
  temp <- c(pred.x, xx)
  lmat <- ns(temp, df = df, intercept = TRUE)
  pmat <- lmat[1:length(pred.x), ]
  xmat <- lmat[-(1:length(pred.x)), ]
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
  periods <- rep(1:n_periods, each = n_days_period)
  ve_dat <- tibble(day = round(pred.x), period = periods, ve = as.numeric(yhat)) %>%
    select(-day) %>%
    group_by(period) %>%
    summarise_all(.funs = mean)
  
  return(ve_dat)
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
tian_ve <- function(dat, n_days, n_periods, n_days_period, alpha = 0.05){
  reject_h0 <- 0
  # fit timecox model
  fit = timecox(Surv(DINF_new, FARI) ~ V, data = dat, max.time = n_days, n.sim = 500)
  KS_pvalue = fit$pval.testBeqC[2]
  if (KS_pvalue < alpha){reject_h0 = 1} 
   # de-cummulative estimates
    if (length(n_days) == 1) {
      predicted.timepts <- 1:n_days #seq(from = min(fit$cum[,1]), to = max(fit$cum[,1]), length = n_days)
    } else {predicted.timepts <- n_days}
    
  # cummulative estimates - we don't want to use these because the are 
  # estimates of the regression coefficients from t=0 to t=j, where j 
  # is the time point, j = 0, ..., T (T is the last time point)
  # rtn_cum <- tibble(time = predicted.timepts, 
  #                   beta_t = fit$cum[,3][which(fit$cum[,1] %in% time)],
  #                   ve = 1 - exp(beta_t))
  # we need to calculate the decumulated regression coefficient estimates
  # then we can calculate the hazard rate ratio
  #   decumulated.ests = CsmoothB(fit$cum, predicted.timepts, b = 1)
  #   timecox.hazardrate = exp(decumulated.ests[,2])
  #   timecox.var = decumulated.ests[,-(1:2)]
  #   rtn <- tibble(time = predicted.timepts, hazardrate = timecox.hazardrate) %>% 
  #             mutate(ve = 1 - hazardrate)
  # # output
  # periods <- rep(1:n_periods, each = n_days_period)
  # ve_dat <- tibble(day = round(pred.timepoints), period = periods, ve = as.numeric(yhat)) %>%
  #   select(-day) %>%
  #   group_by(period) %>%
  #   summarise_all(.funs = mean)
  return(reject_h0)
  #list(output = ve_dat, reject_h0 = reject_h0)
}

# Estimate VE using method from Ainslie et al. 2017 (SIM)
 ainslie_ve <- function(dat, n_days, n_periods, n_days_period, latent_period = 1, infectious_period = 4){
   # calculate prevalence
   N <- length(unique(dat$ID))
   prev <- numeric(n_days)
   for (d in 1:n_days){
     # calculate which days individuals got infected to be infectious on day d
     possible_day_of_infection <- (d - latent_period - infectious_period):(d - latent_period)
     prev[d] <- length(which(dat$DINF_new %in% possible_day_of_infection))/N
   }
   prev <- ifelse(prev == 0, 0.001, prev)
   x <- list(n_days = n_days, prev = prev, dinf = dat$DINF_new, v = dat$V)
   
   logLik <- function(x, par){
     beta = par[1]
     theta_0 = par[2]
     lambda = par[3]
     
      pi_00 <- pi_01 <- pi_10 <- pi_11 <- c(1,rep(0,x$n_days-1)) 
      phi_00 <- phi_01 <- phi_10 <- phi_11 <- c(1,rep(0,x$n_days-1)) 
    # loop over days 
     for (d in 2:x$n_days){
    # theta
      theta <- theta_0 + lambda * d
    # conditional probabilities (pi_jv, where j = infection status, v = vaccination status)
       pi_00[d] <- 1 - beta * x$prev[d]                    # not infected, unvaccinated
       pi_01[d] <- 1 - beta * theta * x$prev[d]            # not infected, vaccinated
       pi_10[d] <- beta * x$prev[d]                        # infected, unvaccinated
       pi_11[d] <- beta * theta * x$prev[d]                # infected, vaccinated
    # unconditional probabilities (phi_jv, j = infection status, v = vaccination status)
       if (d == 1){
         phi_00[d] <- pi_00[d]
         phi_01[d] <- pi_01[d]
         phi_10[d] <- pi_10[d]
         phi_11[d] <- pi_11[d]
       } else {
        phi_00[d] <- pi_00[d] * phi_00[d-1]          # not infected, unvaccinated
        phi_01[d] <- pi_01[d] * phi_01[d-1]          # not infected, vaccinated
        phi_10[d] <- pi_10[d] * phi_00[d-1]          # infected, unvaccinated
        phi_11[d] <- pi_11[d] * phi_01[d-1]          # infected, vaccinated
       }
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
  
      return(-sum(log(Li)))
   }
   # maximum likelihood estimates
    mle <- optim(par = c(0.1, 0.4, 0.5), logLik, x = x, method = "L-BFGS-B", lower = c(0.0001, 0.0001, 0.0001), 
                 upper = c(1, 1, 1), hessian = TRUE)
    
    se <- sqrt(diag(solve(mle$hessian)))
    
    param_est <- tibble(param = c("beta", "theta_0", "lambda"), mle = mle$par, se = se,
                        lower = mle - 1.96 * se, upper = mle + 1.96 * se)
    
    periods <- rep(1:n_periods, each = n_days_period)
    ve_dat <- tibble(day = 1:n_days, period = periods, ve = 1-(mle$par[2] + mle$par[3] * day)) %>%
                  select(-day) %>%
                  group_by(period) %>%
                  summarise_all(.funs = mean)
    
   # output
     rtn <- list(param_est = param_est, ve_dat = ve_dat)

     return(rtn)
 }

