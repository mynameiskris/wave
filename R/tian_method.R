# ----------------------------------------------------------
#' Estimate VE using method from Tian et al. 2005
#'
#' In this method, kernel-weighted partial likelihood approach is used to estimate the time-dependent coefficient
#' in the generalized Cox model [REF] .  At each time point, the estimate is obtained by maximizing a smooth
#' concave function of a p x 1 vector of parameters, where p is the dimension of the vector of covariates. The
#' (1-α) confidence bands for the time-dependent coefficient β(t) can be obtained by β(t) ±c_α w ̂〖(t)〗^(-1
#' ), where t is in the time interval of interest (b_1≤t≤b_2), c_α is the 100(1-α)th percentile the approximate
#' distribution of β ̂(t), and w ̂(t) is a positive weighting function. To test the hypothesis of no VE waning
#' (i.e., the proportional hazards assumption is met), confidence bands from the distribution of the estimated
#' cumulative function of β ̂(t) can be obtained. If the line (0,0) is not contained within the confidence band
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
