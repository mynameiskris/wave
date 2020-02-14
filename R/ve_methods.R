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
