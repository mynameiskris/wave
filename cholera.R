# load necessary packages
library(foreign)
library(survival)
library(splines)
# read in data
data.restore("~/Dropbox/Kylie/Projects/VE Waning/code/chol4/datadmp1")
data.restore("~/Dropbox/Kylie/Projects/VE Waning/code/chol4/datadmp2")
data.restore("~/Dropbox/Kylie/Projects/VE Waning/code/chol4/towork.s")
chol.dat <- read.csv(file="~/Dropbox/Kylie/Projects/VE Waning/code/chol4/choldat.csv", na.strings="")

chol.coxmod <- coxph(Surv(cpdays) ~ cpwc + cpbswc + cpageind, data=chol.dat)
chol.tdhaz <- cox.zph(chol.coxmod, transform = "identity")

plot.cox.zph.ve2 <- function(x, df = 4, nsmo = 40, var,...){
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
    #se
        bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
        xtx <- bk %*% t(bk)
        seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
    
    ylab <- paste("VE(t) for", dimnames(yy)[[2]])
    if (missing(var)) 
        var <- 1:nvar
    else {
        if (is.character(var)) 
            var <- match(var, dimnames(yy)[[2]])
        if (any(is.na(var)) || max(var) > nvar || min(var) < 
            1) 
            stop("Invalid variable requested")
    }
    for (i in var) {
        y <- yy[, i]
        yhat.beta <- pmat %*% qr.coef(qmat, y)
        yhat <- 1-exp(yhat.beta)
        
        #se
            temp <- 2 * sqrt(x$var[i, i] * seval)
            ylow.beta <- yhat.beta - temp
            yup.beta <- yhat.beta + temp
            ylow <- 1-exp(yup.beta)
            yup <- 1-exp(ylow.beta)
            yr <- range(yhat, yup, ylow)
        
            plot(range(xx), c(-0.4,1), type = "n", xlab = "Time", xaxt="n",yaxt="n",
            ylab = ylab[i], ...)
            
            axis(1,at = c(0, 184, 365, 549, 730, 914, 1096, 1280, 
            1461, 1645), label = c("May`85", "Nov`85", "May`86", 
            "Nov`86", "May`87", "Nov`87", "May`88", "Nov`88", 
            "May`89", "Nov`89"))

            axis(2, at = seq(-0.4, 1, 0.2), srt = 90)

        lines(pred.x, yhat)
        abline(h=0,lty=3)
        #se
            lines(pred.x, yup, lty = 2)
            lines(pred.x, ylow, lty = 2)
        
    }
}

plot.cox.zph.ve3 <- function(x, df = 4, nsmo = 40, 
        xaxis.pos = c(0, 184, 365, 549, 730, 914, 1096, 1280, 
            1461, 1645),
        xaxis.label = c("May`85", "Nov`85", "May`86", 
            "Nov`86", "May`87", "Nov`87", "May`88", "Nov`88", 
            "May`89", "Nov`89"), 
            yaxis.lower = -0.4, yaxis.upper = 1, yaxis.ticks = 0.2,   
        var,...){
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
    #se
        bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
        xtx <- bk %*% t(bk)
        seval <- d * ((pmat %*% xtx) * pmat) %*% rep(1, df)
    
    ylab <- paste("VE(t) for", dimnames(yy)[[2]])
    if (missing(var)) 
        var <- 1:nvar
    else {
        if (is.character(var)) 
            var <- match(var, dimnames(yy)[[2]])
        if (any(is.na(var)) || max(var) > nvar || min(var) < 
            1) 
            stop("Invalid variable requested")
    }
    for (i in var) {
        y <- yy[, i]
        yhat.beta <- pmat %*% qr.coef(qmat, y)
        yhat <- 1-exp(yhat.beta)
        
        #se
            temp <- 2 * sqrt(x$var[i, i] * seval)
            ylow.beta <- yhat.beta - temp
            yup.beta <- yhat.beta + temp
            ylow <- 1-exp(yup.beta)
            yup <- 1-exp(ylow.beta)
            yr <- range(yhat, yup, ylow)
        
        
        
            plot(range(xx), c(yaxis.lower,yaxis.upper), type = "n", xlab = "Time", xaxt="n",yaxt="n",
            ylab = ylab[i], ...)
            
            axis(1,at = xaxis.pos, label = xaxis.label)

            axis(2, at = seq(yaxis.lower, yaxis.upper, yaxis.ticks), srt = 90)

        lines(pred.x, yhat)
        abline(h=0,lty=3)
        #se
            lines(pred.x, yup, lty = 2)
            lines(pred.x, ylow, lty = 2)
        return(data.frame(cbind(round(pred.x),round(yhat,digits=2),
        round(ylow,digits=2),round(yup,digits=2))))
    }
}
