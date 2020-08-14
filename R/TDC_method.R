
# ----------------------------------------------------------
#' Estimate VE using Time-dependent covariates (TDC) method
#'
#' This method fits a Cox Proportional hazard regression model with a time-dependent covariate as follows:
#' lambda(t) = lambda_0(t)*Exp(Beta*V+t*V*g(t)), where g(t) is an arbitary function of time; usually the
#' function g(t) = log(t), V is the vaccine status (0 = unvaccinated, 1 = vaccinated)
#'
#' @param dat data set
#' @param alpha hypothesis critical value
#' @return whether or not the null hypothesis is rejected
#' @keywords wave
#' @import survival
#' @import timereg
#' @import dplyr
#' @export
tdc_ve <- function(dat, alpha = 0.05){

  # fit coxph model
  res.logt = coxph(Surv(DINF_new, FARI) ~ V + tt(V), data = dat, tt=function(x,t,...)x*log(t))
  reject_ho = ifelse(summary(res.logt)$coefficient[2,5]<alpha,1,0)
  # reject_h0 = 1 (has TVE),reject_h0=0 has no TVE
  return(reject_h0)

}
# ----------------------------------------------------------

