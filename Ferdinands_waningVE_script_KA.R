##Examining waning immunity in US FLU VE NETWORK data
library(psych)
require(ggplot2)
library(splines)
require(stats)
require(lattice)
library(arm)
library(interplot)
require(MatchIt)
require(survival)
require(matrixStats)

#input dataset
fluve <- read.csv("//cdc.gov/private/M310/zdn5/Flu VE AMBULATORY Network/Waning 2014/pooled/h3pool_22MAR16.csv", head=TRUE)

#define factors
fluve <- fluve %>% mutate(case.f = factor(case, levels = c(0,1), labels = c("controls", "cases")),
                          SVAC_ID.f = factor(SVAC_ID14, levels = c(0,1), labels = c("unvaccinated", "vaccinated")),
                          season.f = factor(season), days = as.numeric(days_flu), days_sq = days_flu*days_flu,
                          days2 = as.numeric(days_sq))

# fluve$case.f <-factor(fluve$case, levels = c(0,1), labels = c("controls", "cases"))
# fluve$SVAC_ID.f <-factor(fluve$SVAC_ID14, levels = c(0,1), labels = c("unvaccinated", "vaccinated"))
# fluve$season.f <-factor(fluve$season)

#Days since start of flu season
# fluve$days <-as.numeric(fluve$days_flu)
# fluve$days_sq <- fluve$days_flu*fluve$days_flu
# fluve$days2 <-as.numeric(fluve$days_sq)


############  CHECKING BALANCE OF CASES AND CONTROLS ON ONSET DATE #############################
histogram(~fluve$onset_percent|fluve$case.f, main="Pooled Influenza A(H3N2) cases and controls", xlab="Onset date percentile (%)", 
          col="grey")


###########  GENERATING RESTRICTED DATASET WITH MATCHING ON ONSET OF ILLNESS ##############################
set.seed(100)
m.out4 = matchit(case ~ onset_percent, data = fluve, method = "nearest" , ratio = 3) 
summary(m.out4)
plot(m.out4, type = "hist") 
m.data4 <- match.data(m.out4)

#Onset percentile (variable onset_percent) is illness onset date of the observation defined 
#relative to percentile of distribution of site-, season- and subtype-specific case onset dates
histogram(~fluve_new$onset_percent|fluve$case.f, main="Pooled Influenza A(H3N2) cases and controls (restricted)", xlab="Onset date percentile (%)", 
          col="grey")

#exporting matchded dataset from matchit
write.csv(m.data4, file = "//cdc.gov/private/M310/zdn5/Flu VE AMBULATORY Network/Waning 2014/pooled/h3pool_nearest3def1.csv")

#inputing matched dataset from matchit
fluve_new <- read.csv("//cdc.gov/private/M310/zdn5/Flu VE AMBULATORY Network/Waning 2014/pooled/h3pool_nearest3def1.csv",head=TRUE)

####################  CASES AND VACCINATION STATUS ######################################
nrow(fluve)
nrow(fluve_new)
table2 <- table(fluve_new$SVAC_ID14, fluve_new$case); table2

#######################################################

# MULTIVARIATE LOGISTIC REGRESSION
# Model with dichotomous vaccination only (no term for time since vaccination)
ve.logr1<- glm(case ~ agele18 + age50to64 + agege65
               + SVAC_ID14
               + any_hr_py + health_not_excel
               + spe2 + spe3
               + onset0to20 + onset20to40 + onset60to80 + onset80to100 
               + UM + UP + SW + GH              
               + female + BLACK + priorvx
               + seas12 + seas13 + seas14 
               + days + days2, family=binomial("logit"), data=fluve_new)
exp(cbind(OR = coef(ve.logr1), confint(ve.logr1)))


###########  DISTRIBUTIONS OF VACCINE INTERVAL FOR DEFINING SPLINES ############################
vacc.all <-subset(fluve_new, SVAC_ID14==1)
quantile(vacc.all$vaccint, probs = seq(0, 1, .1), na.rm = TRUE)
quantile(vacc.all$vaccint, probs = seq(.9, 1, .01), na.rm = TRUE)

qu0 <-quantile(vacc.all$vaccint, probs = 0.10);qu0
qu1 <-quantile(vacc.all$vaccint, probs = 0.98);qu1
qu2 <-quantile(vacc.all$vaccint, probs = 0.99);qu2
qu3 <-quantile(vacc.all$vaccint, probs = 0.80);qu3
qu4 <-quantile(vacc.all$vaccint, probs = 0.20);qu4
qu5 <-quantile(vacc.all$vaccint, probs = 0.50);qu5

#Spline model - days from vaccination to onset modeled as spline
#Specification of spline knots based on comparison of model fit
#for numerous alternative specifications (not shown here)

ve.logr2<- glm(case ~ agele18 + age50to64 + agege65
               + SVAC_ID14 +  ns(vaccint, knots=c(qu1, qu2)) 
               + any_hr_py + health_not_excel
               + spe2 + spe3
               + onset0to20 + onset20to40 + onset60to80 + onset80to100 
               + UM + UP + SW + GH              
               + female + BLACK + priorvx
               + seas12 + seas13 + seas14
               + days + days2, family=binomial("logit"), data=fluve_new)
#exp(cbind(OR = coef(ve.logr2), confint(ve.logr2)))
summary(ve.logr2)

aOR.avg <- exp(ve.logr1$coefficient[5])
aVE.avg <- (1 - aOR.avg)*100; aVE.avg

LH1 <- logLik(ve.logr1);LH1
LH2 <- logLik(ve.logr2);LH2

delta <- 2*(LH2[1] - LH1[1]);delta
#p.val <- pchisq(delta, 4, lower.tail=FALSE);p.val
#p.val <- pchisq(delta, 3, lower.tail=FALSE);p.val

##############################GRAPH#################################
#For VE graph 
newdata2 = data.frame(agele18=0, age50to64=0, agege65=0, any_hr_py=0, health_not_excel=0, spe2=0, spe3=0,
                      onset0to20=0, onset20to40=0, onset60to80=0, onset80to100=0, days=0, days2=0,
                      UM=0, UP=0, SW=0, GH=0, female=0, BLACK=0, priorvx=0, seas12=0, seas13=0, seas14=0, 
                      vaccint=seq(14,180,2), SVAC_ID14=1)

yhat=predict(ve.logr2, newdata2, type="response") 
unvacc <- ve.logr2$coefficient[1]
unvacc.1 <- (exp(unvacc))/(1+exp(unvacc))
b <- rep((unvacc.1), times=84)

aOR <- (yhat/(1-yhat))/(b/(1-b))
aVE <- (1-aOR)*100
plot(newdata2$vaccint, aVE, ylim=c(-50,100), cex.main=0.95, cex.lab = 0.75, cex.sub=0.8, main="VE against influenza A(H3N2) \n NAT CU SPLINE (2 knots, 98th and 99th%) \n n=11,200", 
     type="l", xlab="days between vaccination & onset", ylab="adjusted VE (%)")
max(aVE)
min(aVE)

##Bootstrapped confidence intervals
myfun <- function(){
  srows <- sample(1:nrow(fluve_new),nrow(fluve_new),TRUE)
  glm.out <- glm(case ~ agele18 + age50to64 + agege65
                    + SVAC_ID14 +  ns(vaccint, knots=c(qu1, qu2)) 
                             + any_hr_py + health_not_excel
                             + spe2 + spe3
                             + onset0to20 + onset20to40 + onset60to80 + onset80to100 
                             + UM + UP + SW + GH              
                             + female + BLACK + priorvx
                             + seas12 + seas13 + seas14
                             + days + days2, family=binomial("logit"), data=fluve_new[srows,])
  yhat <- predict(glm.out, newdata = newdata2, "response")
  unvacc <- ve.logr2$coefficient[1]
  unvacc.1 <- (exp(unvacc))/(1+exp(unvacc))
  b <- rep((unvacc.1), times=84)
  
  aOR <- (yhat/(1-yhat))/(b/(1-b))
  aVE <- (1-aOR)*100
  return(aVE)
}
  
set.seed(5)
bootdist1 <- replicate(1000, myfun())
mean(bootdist2)

q <- rowQuantiles(bootdist1, probs=c(0.025, 0.975));q

ci.df <- as.data.frame(q)
ci.df$rownumber = 1:dim(ci.df)[1]
ci.df$vaccint <- 2*(ci.df$rownumber - 1) + 14
names(ci.df)[names(ci.df) == '2.5%'] <- 'ci.lo'
names(ci.df)[names(ci.df) == '97.5%'] <- 'ci.hi'

attach(ci.df)

# Apply loess smoothing using the default span value of 0.8.  
y.loess <- loess(y ~ x, span=0.8, data.frame(x=vaccint, y=ci.lo))
y2.loess <-loess(y2 ~ x, span=0.8, data.frame(x=vaccint, y2=ci.hi))
# Compute loess smoothed values for all points along the curve
y.predict <- predict(y.loess, data.frame(x=vaccint))
y2.predict <- predict(y2.loess, data.frame(x=vaccint))


tiff(file = "Fig 2a H3 VE.tiff", width = 4600, height = 3200, units = "px", res = 800)
par(bty = 'l') 
plot(newdata2$vaccint, aVE, ylim=c(-50,100), lwd=2, cex.main=0.80, cex.lab = 0.75, cex.sub=0.7, 
     main="Vaccine effectiveness against influenza A(H3N2) (n=11,200) \n p=0.004", 
     sub="", type="l", xlab="days between vaccination and onset", 
     ylab="adjusted VE (%)")
lines(vaccint, y2.predict, type="l", col="black", lty=3)
lines(vaccint,y.predict, type="l", col="black", lty=3)
#legend("topleft", inset=0.03, legend=c("VE", "95% CI"),
       #col=c("black", "black"), lty=1:3, cex=0.65, box.lty=0)
legend("bottomleft", inset=0.02, legend=c("VE", "95% CI"),
       col=c("black", "black"), lty=1:3, cex=0.65, box.lty=0)

dev.off() 
detach(ci.df)

