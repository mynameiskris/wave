
# ----------------------------------------------------------
#' Estimate VE using Time-dependent covariates (TDC) method
#'
#' This method fits a Cox Proportional hazard regression model with a time-dependent covariate as follows:
#' lambda(t) = lambda_0(t)*Exp(Beta*V+gamma*V*g(t)), where g(t) is an arbitary function of time; usually the
#' function g(t) = log(t), V is the vaccine status (0 = unvaccinated, 1 = vaccinated)
#'
#'
#' @param dat data set
#' @param n_days number of days in the study
#' @param alpha hypothesis critical value
#' @return number of times the null hypothesis is rejected
#' @keywords wave
#' @import survival
#' @import timereg
#' @import dplyr
#' @import tidyr
#' @export
tdc_ve <- function(dat, n_days, n_periods, n_days_period, alpha = 0.05){
  reject_h0 <- 0
  # fit timecox model
  gt = log(DINF_new)
  fit = timecox(Surv(DINF_new, FARI) ~ const(V)+ V*gt, data = dat, max.time = n_days, n.sim = 500)
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
