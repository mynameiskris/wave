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
   
   N <- length(unique(dat$ID))
   prev <- numeric(n_days)
   
   for (d in 1:n_days){
     # calculate which days individuals got infected to be infectious on day d
     possible_day_of_infection <- (d  - latent_period - infectious_period):(d - latent_period)
     prev[d] <- length(which(dat$DINF_new %in% possible_day_of_infection))/N
   }
   prev <- ifelse(prev == 0, 0.001, prev)
   x <- list(n_days = n_days, 
             prev = prev, 
             dinf = dat$DINF_new, 
             v = dat$V
   )
   
   logLik <- function(x, par){
     alpha = par[1]
     theta_0 = par[2]
     lambda = par[3]
     
     pi_0u <- pi_0v <- pi_1u <- pi_1v <- c(1,rep(0,x$n_days-1)) 
     psi_0u <- psi_0v <- psi_1u <- psi_1v <- pi_0u 
     # loop over days 
     for (d in 2:x$n_days){
       # theta
       theta_d <- theta_0 + lambda * d
       alpha_d <- alpha * x$prev[d]
       # conditional probabilities: pi_ju & pi_jv, where 
       #   j = infection status, 
       #   u = unvaccinated, 
       #   v = vaccinated
       pi_0u[d] <- ifelse(1 - alpha_d < 0, 0.0001, 1 - alpha_d)           
       pi_0v[d] <- ifelse(1 - alpha_d * theta_d < 0, 0.0001, 1 - alpha_d * theta_d)
       pi_1u[d] <- alpha_d       
       pi_1v[d] <- ifelse(alpha_d * theta_d > 1, 1, alpha_d * theta_d)
       # unconditional probabilities: psi_ju & psi_jv, where 
       #   j = infection status,
       #   u = unvaccinated,
       #   v = vaccinated
       if (d == 1){
         psi_0u[d] <- pi_0u[d]
         psi_0v[d] <- pi_0v[d]
         psi_1u[d] <- pi_1u[d]
         psi_1v[d] <- pi_1v[d]
       } else {
         psi_0u[d] <- pi_0u[d] * psi_0u[d-1]          # not infected, unvaccinated
         psi_0v[d] <- pi_0v[d] * psi_0v[d-1]          # not infected, vaccinated
         psi_1u[d] <- pi_1u[d] * psi_0u[d-1]          # infected, unvaccinated
         psi_1v[d] <- pi_1v[d] * psi_0v[d-1]          # infected, vaccinated
       }
     }
     # personal contribution to the likelihood ---------------------------------
     Li <- numeric(N)
     
     for (i in 1:N){
       if( x$dinf[i] == 999 ){
         if( x$v[i] == 0 ){ Li[i] <- psi_0u[x$n_days] }
         if( x$v[i] == 1 ){ Li[i] <- psi_0v[x$n_days] }
       }
       if ( x$dinf[i] != 999 ){
         if( x$v[i] == 0 ){ Li[i] <- psi_1u[x$dinf[i]] }
         if( x$v[i] == 1 ){ Li[i] <- psi_1v[x$dinf[i]] }
       }
     }
     
     return(-sum(log(Li)))
   }
   
   # use DE optim to get initial values
   initial <- DEoptim(fn=logLik, 
                      x = x,
                      lower = c(0.0001, 0.0001, 0.0001), 
                      upper = c(1, 1, 1),
                      control = list(trace = FALSE)
   )
   # maximum likelihood estimates ----------------------------------------------
   mle <- optim(par = initial$optim$bestmem, 
                fn = logLik, 
                x = x, 
                method = "L-BFGS-B", 
                lower = c(0.0001, 0.0001, 0.0001), 
                upper = c(1, 1, 1), 
                hessian = TRUE
   )  
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

