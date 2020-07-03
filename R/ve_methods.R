# VE funtions
# ----------------------------------------------------------


# ----------------------------------------------------------
#' Estimate VE using method from Durham et al. 1988 (AJE)
#'
#' This method is based on smoothing scaled residuals from the Cox proportional hazard regression model.
#' It consists of four steps. First, an ordinary proportional hazard model is fitted using a partial likelihood
#' function. Second, Schoenfeld residuals are calculated. These residuals are used to test the independence
#' between residuals and time. Third, the residuals are scaled and added to the coefficient from the ordinary
#' proportional hazard regression model. Fourth, after smoothing we can get the estimated hazard ratios as
#' function of time by recovering the time-varying regression coefficient. This allows testing the hypothesis of
#' no VE waning and the estimation of TVE at each time point.
#'
#' @param x survival object
#' @param df degrees of freedom
#' @param n_days number of days in the study
#' @param n_periods number of periods in the study
#' @param n_days_period number of days per period
#' @param var name of variable to assess
#' @return tibble of VE estimates for each period (the estimate for each period is the average VE over the days
#' within each period.
#' @keywords wave
#' @import splines
#' @import dplyr
#' @import tidyr
#' @export
durham_ve <- function(x, df = 2, n_days, n_periods, n_days_period, var){
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
# ----------------------------------------------------------

#' Estimate VE using method from Ferdinands et al. 2017
#'
#' @param dat data set.
#' @param splines logical. whether or not to include calendar time as a spline.
#' @return estimates of the likelihood ratio under each model and the difference between the two models.
#' @keywords wave
#' @import splines
#' @import dplyr
#' @import tidyr
#' @export
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
# ----------------------------------------------------------

# ----------------------------------------------------------
#' Estimate VE using method from Tian et al. 2005
#'
#' In this method, kernel-weighted partial likelihood approach is used to estimate the time-dependent coefficient
#' in the generalized Cox model [REF]. At each time point, the estimate is obtained by maximizing a smooth
#' concave function of a p x 1 vector of parameters, where p is the dimension of the vector of covariates. The
#' (1-alpha) confidence bands for the time-dependent coefficient beta(t) can be obtained by beta(t) ±c_alpha w [(t)]^(-1
#' ), where t is in the time interval of interest (b_1 <=t<=b_2), c_alpha is the 100(1-alpha)th percentile the approximate
#' distribution of beta(t)_hat, and w(t)_hat is a positive weighting function. To test the hypothesis of no VE waning
#' (i.e., the proportional hazards assumption is met), confidence bands from the distribution of the estimated
#' cumulative function of beta(t)_hat) can be obtained. If the line (0,0) is not contained within the confidence band
#' s over the entire time interval, then the proportional hazards assumption is not met and the null hypothesis
#' of no VE waning is rejected.
#' @param dat data set
#' @param n_days number of days in the study
#' @param n_periods number of periods in the study
#' @param n_days_period number of days per period
#' @param alpha hypothesis critical value
#' @return number of times the null hypothesis is rejected
#' @keywords wave
#' @import survival
#' @import timereg
#' @import dplyr
#' @import tidyr
#' @export
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
# ----------------------------------------------------------

# ----------------------------------------------------------
#' Estimate VE using a maximum likelihood (ML) method
#'
#' This method is similar to that of Ainslie et al. 2017 (SIM). The ML method is based on calculating the
#' contribution to the likelihood function of each study participant. Let t (t=1,2,…T) denote the number
#' of days since vaccination or receipt of the placebo. Denote lambda_d = hazard of infection of an unvaccinated
#' susceptible person on day t.P(t) denotes the probability that a randomly selected population member
#' is infectious on day t  (prevalence of infection). We assume that the hazard of infection is proportional
#' to the prevalence: lambda_t=lambda*P(t). A method for estimatingon the prevalence is described later. Next, let theta_t
#' be defined such that the hazard of a vaccinated susceptible person on day t is lambda_t*theta_t. Then the TVE
#' on day t, defined as one minus the ratio of the hazards of infection for an unvaccinated and a
#' vaccinated person on that day is [1-theta]_t. Next, we define a status variable for each study participant
#' i (i=1,2,…N) on each day t:
#' * Y_it=0 if this person remained susceptible by the end of the day,
#' * Y_it=1 if this person became infected on this day, and
#' * Y_it=2 if this person became infected before this day.
#' (A person with Y_it=2 may still be infectious to others; we use this notation to indicate that this person
#' is no longer at risk of becoming infected). We assume that on day 0 everyone is susceptible, i.e., Y_i0=0.
#'
#' Next, let  pi_itj=Pra({Y_it=j|Y_i(t-1)=0)}  be the conditional probability that person i is in status j
#' (j=0,1,2) given that s/he was susceptible at the end of the previous day. If person i is unvaccinated then
#' pi_it0=1-lambda_t,  pi_it1=lambda_t,  pi_it2=0. If person i is vaccinated then lambda_t in the last equations is replaced by
#' lambda_t*theta_t. Denoting by psi_itj=Pra(a{Y_it=j)} the unconditional probabilities of Y_it, it is easy to see tha
#' t psi_it0= psi_(i(t-1)0)*pi_it0, psi_it1=psi_(i(t-1)0)pi_it1, psi_it2= 1-(psi_it0+psi_it1).
#'
#' We now can write the contribution of each study subject to the likelihood function in terms of the
#' unconditional probabilities psi_itj. The contribution of a person who became infected on day t is
#' [psi]_it1, while the contribution of a person who did not become infected during the entire study
#' is psi_iT0, where T is the last day of the study. For a participant who dropped out uninfected, T is replaced
#' by the day s/he dropped out. We can assume that the contributions of different study participants are
#' independent because we condition on the daily prevalences when we determine the unconditional distribution
#' of each Y_it. Therefore, the final likelihood function is the product of the contributions of all the study
#' participants.
#'
#' To evaluate waning of TVE, we recall that TVE on day t is [1-theta]_t. Therefore, temporal changes in TVE
#' are determined by how theta_t depends on t. A simple approach to modeling TVE as a function of time assumes a
#' fixed daily increase phi>=0 in theta_t (which amounts to a fixed daily decrease in TVE), i.e. theta_t=theta_1+phi*(t-1)
#' for t=1,2,…T. Then the likelihood function has 3 parameters: alpha,theta_1,and phi. One can calculate the maximum
#' likelihood estimated of these parameters along with their standard errors by maximizing the likelihood
#' function. One can then test the hypothesis phi=0 to determine the presence of waning and obtain a confidence
#' interval for phi, the daily rate of TVE waning.
#' @param dat data set
#' @param n_days number of days in the study
#' @param n_periods number of periods in the study
#' @param n_days_period number of days per period
#' @param latent_period length of latent period
#' @param infectious_period length of infectious period
#' @return list with a tibble of VE estimates for each period (the estimate for each period is the average VE over the days
#' within each period) and the maximu likelihood estimates of the parameters.
#' @keywords wave
#' @import dplyr
#' @import tidyr
#' @import xlsx
#' @export
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
     phi = par[3]

     pi_0u <- pi_0v <- pi_1u <- pi_1v <- c(1,rep(0,x$n_days-1))
     psi_0u <- psi_0v <- psi_1u <- psi_1v <- pi_0u

     #initialise period
     period <- 1
     period_start_days <- seq(1, params$ND, by = params$NDJ)
     # loop over days
     for (d in 2:x$n_days){
       if(d %in% period_start_days){period <- period + 1}
       #print(period)
       lambda <- phi - 1
       theta_d <- theta_0 + lambda * period
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
   # initial <- DEoptim(fn=logLik,
   #                    x = x,
   #                    lower = c(0.0001, 0.0001, 0.0001),
   #                    upper = c(1, 1, 2),
   #                    control = list(itermax = 100, trace = FALSE)
   # )
   # print(initial$optim$bestmem)
   # maximum likelihood estimates ----------------------------------------------
   #tryCatch({

   mle <- optim(par = c(0.3, 0.4, 1),
                fn = logLik,
                x = x,
                method = "L-BFGS-B",
                lower = c(0.0001, 0.0001, 0.0001),
                upper = c(1, 1, 2),
                hessian = TRUE
                # control = list(trace = 3,
                #                maxit = 1000,
                #                ndeps = 1e-4)
   )
   #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    se <- sqrt(diag(solve(mle$hessian)))

    param_est <- tibble(param = c("alpha", "theta_0", "lambda"), mle = mle$par - c(0,0,1), se = se,
                        lower = mle - 1.96 * se, upper = mle + 1.96 * se)

    periods <- rep(1:n_periods, each = n_days_period)
    ve_dat <- tibble(day = 1:n_days, period = periods, ve = 1-(mle$par[2] + (mle$par[3] - 1) * day)) %>%
                  select(-day) %>%
                  group_by(period) %>%
                  summarise_all(.funs = mean)

   # output
     rtn <- list(param_est = param_est, ve_dat = ve_dat)

     return(rtn)
 }
# ----------------------------------------------------------

# ----------------------------------------------------------
#' Apply different VE estimation methods
#'
#' This function applies teh Durham method, Tian method, and ML method to data from a vaccine
#' efficacy/effectiveness study.
#' @param dat data set
#' @param params list of input parameters
#' @param write_to_file logical. If true, outputs are written to file.
#' @param path path where files are to be written. Defaults to working directory
#' @return list of VE estimates from each method, maximum likelihood parameter estimates from the ML method for
#' each simulation, the proportion of times the null hypothesis was rejected for each method, mean VE estimate over simulations
#' for each time period, and the mean maximum likelihood parameters over all simulations.
#' @keywords wave
#' @import survival
#' @import dplyr
#' @import tidyr
#' @export
estimate_ve <- function(dat = outcomes_dat, params, write_to_file = TRUE, path = getwd()){

# initialise count of number of simulations in which H0 is rehected --------------------------------------
   reject_h0_durham <- reject_h0_tian <- reject_h0_ainslie <- 0

# loop through simulations and apply each method ---------------------------------------------------------
   for (i in 1:max(outcomes_dat$Sim)){
     print(i)

    # subset data for each simulation---------------------------------------------------------------------
      outcomes_dat1 <- outcomes_dat %>% filter(Sim == i)

# method from Durham et al. 1988 -------------------------------------------------------------------------

# fit ordinary Cox propotional hazards model -------------------------------------------------------------
  flu_coxmod <- coxph(Surv(DINF_new,FARI) ~ V, data=outcomes_dat1)

# test the proportional hazards assumption and compute the Schoenfeld residuals ($y) ---------------------
  flu_zph <- cox.zph(fit = flu_coxmod, transform = "identity")
  reject_h0_durham <- reject_h0_durham + ifelse(flu_zph$table[1,3] < 0.05, 1, 0)

# calculate VE -------------------------------------------------------------------------------------------
# the nsmo argument indicates the number of time points to calculate VE at
  temp <- durham_ve(flu_zph, n_days = params$ND, n_periods = params$NJ,
                    n_days_period = params$NDJ,var = "V") %>%
    mutate(Sim = i, Method = "Durham")

# bind results together ----------------------------------------------------------------------------------
  if (i > 1){
    ve_est <- bind_rows(ve_est,temp)
  } else {ve_est <- temp}


# method from Tian et al. 2005 ---------------------------------------------------------------------------

     # calculate VE
     # n_timepoint_breaks argument specifies the number of time points to calculate VE for
     temp2 <- tian_ve(outcomes_dat1, n_days = params$ND, n_periods = params$NJ,
                      n_days_period = params$NDJ)
     #temp2a <- temp2$output %>% mutate(Sim = i, Method = "Tian")
     #ve_est <- bind_rows(ve_est,temp2a)
     # proportion of sims where null hypothesis is rejected
     reject_h0_tian <- reject_h0_tian + temp2 #temp2$reject_h0


# method from Ainslie et al. 2017 ------------------------------------------------------------------------

     temp3 <- ainslie_ve(outcomes_dat1, n_days = params$ND, n_periods = params$NJ,
                         n_days_period = params$NDJ, latent_period = 1, infectious_period = 4)
     temp3a <- temp3$ve_dat %>% mutate(Sim = i, Method = "Ainslie")
     ve_est <- bind_rows(ve_est,temp3a)

     # proportion of sims where null hypothesis is rejected
     reject_h0_ainslie <- reject_h0_ainslie + ifelse(temp3$param_est$lower[3] > 0, 1, 0)

     # mle parameter estimates
     temp3b <- temp3$param_est %>% mutate(Sim = i, Method = "Ainslie")
     if (i > 1){
       mle_param_est <- bind_rows(mle_param_est,temp3b)
     } else {mle_param_est <- temp3b}
   }

# outputs ------------------------------------------------------------------------------------------------
   prop_reject_h0 <- tibble(method = c('Durham','Tian', 'Ainslie'),
                            count = c(reject_h0_durham, reject_h0_tian,reject_h0_ainslie),
                            proportion = count/params$sim)
   mean_ve <- ve_est %>%
     group_by(Method, period) %>%
     summarise_at("ve", c(mean, sd)) %>%
     rename(ve_mean = fn1, ve_sd = fn2)

   mean_mle_params <- mle_param_est %>%
     group_by(param) %>%
     summarise_at(.vars = "mle", .funs = "mean")

   rtn <- list(ve_est = ve_est,
               mle_param_est = mle_param_est,
               prop_reject_h0 = prop_reject_h0,
               mean_ve = mean_ve,
               mean_mle_params = mean_mle_params)

   if(write_to_file){
     write.xlsx(ve_est, file="output.xlsx",sheetName="ve_estimates", row.names=F )
     write.xlsx(ve_est, file="output.xlsx",sheetName="mle_parameter_estimates_",append = T, row.names=F )
     write.xlsx(ve_est, file="output.xlsx",sheetName="reject_h0_prop_",append = T, row.names=F )
     write.xlsx(ve_est, file="output.xlsx",sheetName="mean_ve",append = T, row.names=F )
     write.xlsx(ve_est, file="output.xlsx",sheetName="mle_parameter_estimates_",append = T, row.names=F )
   }

   return(rtn)
 }
# ----------------------------------------------------------
