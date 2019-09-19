

# Source functions from other files
source("readParams.R")
set.seed(1234)
# Read parameters from input files
params <- readParams("SimVEE_input.csv")

## Filenames
#output_files = c('Outcomes.csv', 'detailed.csv')

ID <- seq(1, params$N)
X <-  rep(0, params$N)
V <-  rep(0, params$N)
DINF <-  rep(0, params$N)

subject <- data.frame(ID=ID, X=X, V=V, DINF=DINF)

## SubjectY goes from d=0 to d=ND
subjectY <- matrix(NA, nrow=params$N, ncol=params$ND+1)

NDINF = array(0, dim=c(2, 2, params$ND))

### Initialize
subjectY[,1] <- rep(0, params$N) 

# Day
d = 0

params$beta_d11 = params$beta_d01 * params$theta_d
params$beta_d00 = params$beta_d01 * params$phi
params$beta_d10 = params$beta_d01 * params$theta_d * params$phi


## Set value of X for each subject
subject$X = as.numeric(runif(params$N) < params$pai) 

## Set vaccination status for each subject
for (i in ID) {
  if (subject[i,"X"] == 0) {
    subject[i,"V"] = as.numeric(runif(1) < params$alpha_0)
  } 
  if (subject[i,"X"] == 1) {
    subject[i,"V"] = as.numeric(runif(1) < params$alpha_1)
  }
}

## This function returns the beta value that should be applied to a
#  subject based on X and V values.
getBetaForSubject = function (sub) {
  X = sub$X
  V = sub$V
  if (X == 0 && V == 0) return(params$beta_d00)
  if (X == 0 && V == 1) return(params$beta_d01)
  if (X == 1 && V == 0) return(params$beta_d10)
  if (X == 1 && V == 1) return(params$beta_d11)
}


## Iterate over number of days.
while (d < params$ND) {
  d = d + 1
  
  for (i in ID) {
    if (subjectY[i, d] == 0) {
      subjectY[i, (d+1)] = as.numeric(runif(1) < getBetaForSubject(subject[i,])) 
      if (subjectY[i, (d+1)] == 1) {
        subject[i, "DINF"] = d
        NDINF[subject[i,"X"], subject[i,"V"], d] = NDINF[subject[i,"X"], subject[i,"V"], d] + 1
      }
    }
    if (subjectY[i, d] == 1) {
      subjectY[i, (d+1)] = 2
    }
    if (subjectY[i, d] == 2) {
      subjectY[i, (d+1)] = 2
    }
  }
}

if (params$detailed_file == TRUE) {
  detailed <- cbind(subject, subjectY)
  names(detailed) <- c(names(subject), paste0("D",seq(0,params$ND)))
}

if (params$csv == TRUE) {
  if (params$population_report_file == TRUE) {
    write.csv(subject, paste0('Outcomes_',params$title,'.csv'), row.names = FALSE)
  }
  if (params$detailed_file == TRUE) {
    write.csv(detailed, paste0('Detailed_',params$title,'.csv'), row.names = FALSE)
  }
} 

if (params$sas == TRUE) {
  library(haven)
  if (params$population_report_file == TRUE) {
    write_sas(subject, paste0('Outcomes_',params$title,'.sas7bdat'))
  }
  if (params$detailed_file == TRUE) {
    write_sas(detailed, paste0('Detailed_',params$title,'.sas7bdat'))
  }
}


