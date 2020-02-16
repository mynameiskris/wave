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
  rtn <- tibble(variable = rep(var,length(pred.x)),
                pred_x = round(pred.x), 
                y_hat = as.numeric(yhat),
                y_lower = as.numeric(ylow),
                y_upper = as.numeric(yup))
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
  
  LH1 <- logLik(ve.logr1);LH1
  LH2 <- logLik(ve.logr2);LH2
  
  delta <- 2*(LH2[1] - LH1[1]);delta
  
  
}  
}