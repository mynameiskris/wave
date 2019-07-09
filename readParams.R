readParams = function(filename) {
  para <- list()
  Y <- readLines(filename)
  
  para$title <- strsplit(Y[grep("Title",Y)],',')[[1]][2]
  para$sim <- as.numeric(strsplit(Y[grep("Number of simulations",Y)],',')[[1]][2])
  para$rseed <- as.numeric(strsplit(Y[grep("Seed",Y)],',')[[1]][2])
  para$all_or_no_vaccine <- as.numeric(strsplit(Y[grep("Leaky or all-or-none vaccine",Y)],',')[[1]][2])
  
  para$output_format <- as.character(strsplit(Y[grep("Output file format",Y)],',')[[1]][2]) 
  para$inputsfile <- ifelse(strsplit(Y[grep("Input and calculated parameters",Y)],',')[[1]][2]=="yes",TRUE,FALSE) 
  para$detailed_file <- ifelse(strsplit(Y[grep("Detailed for each simulation",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$daily_each_sim_file <- ifelse(strsplit(Y[grep("Incidence-daily for each simulation",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$period_each_sim_file <- ifelse(strsplit(Y[grep("Incidence-period for each simulation",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$seasonal_each_sim_file <- ifelse(strsplit(Y[grep("Incidence-total for each simulation",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$daily_overall_file <- ifelse(strsplit(Y[grep("Incidence-daily means over all simulations",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$monthly_overall_file <- ifelse(strsplit(Y[grep("Incidence-period means over all simulations",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$seasonal_overall_file <- ifelse(strsplit(Y[grep("Incidence-total means over all simulations",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$population_report_file <- ifelse(strsplit(Y[grep("Outcomes for each simulation",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  para$timestamp <- ifelse(strsplit(Y[grep("Add timestamp to output file names",Y)],',')[[1]][2]=="yes",TRUE,FALSE)
  
  para$NDJ <- as.numeric(strsplit(Y[grep("Number of days in each period",Y)],',')[[1]][2])
  para$NJ <- as.numeric(strsplit(Y[grep("Number of periods in the study",Y)],',')[[1]][2])
  para$ND <- as.numeric(strsplit(Y[grep("Number of days in the study",Y)],',')[[1]][2])
  if (para$ND != para$NDJ * para$NJ) {
    stop("Number of days in the study isn't a multiple of number of periods")
  }
  para$N <- as.numeric(strsplit(Y[grep("Number of participants in study",Y)],',')[[1]][2])
  para$pai <- as.numeric(strsplit(Y[grep("Probability of X=1",Y)],',')[[1]][2])
  para$alpha_1 <- as.numeric(strsplit(Y[grep("Probability of vaccination for X=1",Y)],',')[[1]][2])
  para$alpha_0 <- as.numeric(strsplit(Y[grep("Probability of vaccination for X=0",Y)],',')[[1]][2])
  para$beta_d01 <- as.numeric(strsplit(Y[grep("Daily probability of infection for V=0 X=1",Y)],',')[[1]][2])
  para$theta_dx <- as.numeric(strsplit(Y[grep("Multiplier for daily probability of infection for V=1",Y)],',')[[1]][2])
  para$phi <- as.numeric(strsplit(Y[grep("Multiplier for daily probability of infection for X=0",Y)],',')[[1]][2])
  
  
  return(para)
}