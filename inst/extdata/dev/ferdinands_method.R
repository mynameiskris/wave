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
