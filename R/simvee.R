
## Filenames
#output_files = c('Outcomes.csv', 'detailed.csv')

# main simulation function
simvee <- function(params, simNum) {
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
  d_period = 0
  period = 1
  
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
    d_period = d_period + 1
    
    # print(paste("Day",d, "Period", period, "Day within period", d_period))

    for (i in ID) {
      if (subjectY[i, d] == 0) {
        subjectY[i, (d+1)] = as.numeric(runif(1) < getBetaForSubject(subject[i,])[period]) 
        if (subjectY[i, (d+1)] == 1) {
          subject[i, "DINF"] = d
          NDINF[(subject[i,"X"]+1), (subject[i,"V"]+1), d] = 
            NDINF[(subject[i,"X"]+1), (subject[i,"V"]+1), d] + 1
        }
      }
      if (subjectY[i, d] == 1) {
        subjectY[i, (d+1)] = 2
      }
      if (subjectY[i, d] == 2) {
        subjectY[i, (d+1)] = 2
      }
    }
    if(d_period == params$NDJ){
      d_period = 0
      period = period + 1
    }
  }
  return(list(subject=subject,subjectY=subjectY,NDINF=NDINF))
}

# function to run simulation and generate incidence report
run_simvee <- function(params, path = getwd()){
  for (i in 1:params$sim) {
    results <- simvee(params, i)
  
    if (params$population_report_file == TRUE) {
      temp_population_report <- cbind("Sim"=i, results$subject)
      if (i > 1) {
        population_report <- rbind(population_report, temp_population_report)
      }
      else {
        population_report <- temp_population_report
      }
    }
  
    if (params$detailed_file == TRUE) {
      temp_detailed <- cbind(i, results$subject, results$subjectY)
      colnames(temp_detailed) <- c("sim", names(results$subject), paste0("D",seq(0,params$ND)))
      if (i > 1) {
        detailed = rbind(detailed, temp_detailed)
      }
      else {
        detailed = temp_detailed
      }
    }
  
    day_period <- rep(seq(1:params$NJ), each=params$NDJ)
    temp_daily <- cbind(i, seq(1:params$ND), day_period, t(apply(results$NDINF, 3L, c)))
    colnames(temp_daily) <- c("Sim", "Day", "Period", "X0V0", "X1V0", "X0V1", "X1V1")
    if (i > 1) {
      incidence = rbind(incidence, temp_daily)
    }
    else {
      incidence = temp_daily
    }
  }
  # Tibble for dplyr
  g_inc <- as_tibble(incidence)

  # write simulation results to different output files
  if (params$csv == TRUE) {
    if (params$population_report_file == TRUE) {
      write.csv(population_report, paste0(path,'/Outcomes_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$detailed_file == TRUE) {
      write.csv(detailed, paste0(path,'/Detailed_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$daily_each_sim_file == TRUE) {
      write.csv(incidence, paste0(path,'/Daily_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$daily_overall_file == TRUE) {
      inc_daily_overall <- g_inc %>% group_by(Day,Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(-Sim)
      write.csv(inc_daily_overall, paste0(path,'/Daily_overall_',params$title,'.csv'), row.names = FALSE)
    }
    # if (params$period_each_sim_file == TRUE) {
    #   inc_period_sim <- g_inc %>% group_by(Sim,Period) %>% summarise_all(sum, na.rm = TRUE) %>%
    #                           select(-Day) %>% ungroup() %>% group_by(Period) %>%
    #                           summarise_all(mean, na.rm = TRUE) %>% select(-Sim)
    #   write.csv(inc_period_sim, paste0(path,'Period_',params$title,'.csv'), row.names = FALSE)
    # }
    # if (params$period_overall_file == TRUE) {
    #   inc_period_overall <- g_inc %>% group_by(Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(c(-Day,-Sim))
    #   write.csv(inc_period_overall, paste0(path,'Period_overall_',params$title,'.csv'), row.names = FALSE)
    # }
    if (params$seasonal_each_sim_file == TRUE) {
      inc_seasonal_sim <- g_inc %>% group_by(Sim) %>% summarise_all(sum, na.rm = TRUE) %>% select(c(-Day,-Period))
      write.csv(inc_seasonal_sim, paste0(path,'/Seasonal_',params$title,'.csv'), row.names = FALSE)
    }
    if (params$seasonal_overall_file == TRUE) {
      inc_seasonal_overall <- g_inc %>% group_by(Sim) %>% summarise_all(sum, na.rm = TRUE) %>% 
                              select(c(-Day,-Period)) %>%  summarize_all(mean, na.rm=TRUE) %>%
                              select(c(-Sim))
      write.csv(inc_seasonal_overall, paste0(path,'/Seasonal_overall_',params$title,'.csv'), row.names = FALSE)
    }
  } 
  # return incidence tibble
  return(population_report)
}
# 
# if (params$sas == TRUE) {
#   library(haven)
#   if (params$population_report_file == TRUE) {
#     write_sas(population_report, paste0('Outcomes_',params$title,'.sas7bdat'))
#   }
#   if (params$detailed_file == TRUE) {
#     write_sas(detailed, paste0('Detailed_',params$title,'.sas7bdat'))
#   }
#   if (params$daily_each_sim_file == TRUE) {
#     write_sas(incidence, paste0('Daily_',params$title,'.sas7bdat'))
#   }
#   if (params$daily_overall_file == TRUE) {
#     inc_daily_overall <- g_inc %>% group_by(Day,Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(-Sim)
#     write_sas(data.finc_daily_overall, paste0('Daily_overall_',params$title,'.sas7bdat'))
#   }
#   if (params$period_each_sim_file == TRUE) {
#     inc_period_sim <- g_inc %>% group_by(Sim,Period) %>% summarise_all(sum, na.rm = TRUE) %>%
#       select(-Day) %>% ungroup() %>% group_by(Period) %>%
#       summarise_all(mean, na.rm = TRUE) %>% select(-Sim)
#     write_sas(inc_period_sim, paste0('Period_',params$title,'.sas7bdat'))
#   }
#   if (params$period_overall_file == TRUE) {
#     inc_period_overall <- g_inc %>% group_by(Period) %>% summarise_all(mean, na.rm = TRUE) %>% select(c(-Day,-Sim))
#     write_sas(inc_period_overall, paste0('Period_overall_',params$title,'.sas7bdat'))
#   }
#   if (params$seasonal_each_sim_file == TRUE) {
#     inc_seasonal_sim <- g_inc %>% group_by(Sim) %>% summarise_all(sum, na.rm = TRUE) %>% select(c(-Day,-Period))
#     write_sas(inc_seasonal_sim, paste0('Seasonal_',params$title,'.sas7bdat'))
#   }
#   if (params$seasonal_overall_file == TRUE) {
#     inc_seasonal_overall <- g_inc %>% group_by(Sim) %>% summarise_all(sum, na.rm = TRUE) %>% 
#       select(c(-Day,-Period)) %>%  summarize_all(mean, na.rm=TRUE) %>%
#       select(c(-Sim))
#     write_sas(inc_seasonal_overall, paste0('Seasonal_overall_',params$title,'.sas7bdat'))
#   }
# }



