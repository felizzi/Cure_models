require(flexsurv)
require(MASS)

trial_data <- read.csv(paste0("libraries/example.csv"))
source("functions/hazard_script.R")
source("functions/likelihood_funs.r")
hazard_out <- hazard_time(table_ = trial_data, 
                          evttme = "TIME" , 
                          sex = "SEX", 
                          age = "AGE", 
                          year = "YEAR", 
                          country_trial = "COUNTRY", 
                          country_output = NULL)

trial_data$rate_mod <- hazard_out$rate_vec

models <- c("exponential", "weibull", "llogis", "lognormal", "gompertz", "gamma", "gengamma", "genf")

N_m <- length(models) 

fsr_fits_ <- list()
opt_obj_ <- list()
est_cure <- data.frame(function_type = character(N_m), cure_rate = numeric(N_m), stderr = numeric(N_m), lower.ci = numeric(N_m), upper.ci = numeric(N_m), AIC = numeric(N_m), BIC = numeric(N_m), stringsAsFactors = FALSE)

#define the general purpose ll.mix function 


for (i in 1:length(models)){
  #perform the fits
  print(models[i])
  fsr_fits_[[i]] <- flexsurvreg(Surv(as.numeric(TIME), CNSR)~1, data = trial_data, dist = models[i])
  fsr_use <- fsr_fits_[[i]]$res[,'est'];
  opt_obj_[[i]] <- try(optim(par=c(fsr_use, 0.23), #dummy staring parameter for the cure-fraction 
                            ll.mix, 
                            table = trial_data,
                            time_col = "TIME",
                            cen_col = "CNSR",
                            rate_col = "rate_mod",
                            method = "L-BFGS-B", 
                            obj = fsr_fits_[[i]],
                            hessian = TRUE, 
                            lower=c(rep(-Inf,length(fsr_use)),0), upper = c(rep(Inf,length(fsr_use)),1)) , silent = FALSE)
  
  invHuse <- ginv(opt_obj_[[i]]$hessian)
  est_cure$function_type[i] <- models[i]
  est_cure$cure_rate[i] <- opt_obj_[[i]]$par[length(fsr_use) + 1]
  est_cure$stderr[i] <-  sqrt(invHuse[length(fsr_use) + 1, length(fsr_use) + 1])
  est_cure$lower.ci[i] <- est_cure$cure_rate[i] - qt(0.975, fsr_fits_[[i]]$N)*est_cure$stderr[i]
  est_cure$upper.ci[i] <- est_cure$cure_rate[i] + qt(0.975, fsr_fits_[[i]]$N)*est_cure$stderr[i]
  est_cure$AIC[i] <- 2*(length(fsr_use) + 1) + 2*opt_obj_[[i]]$value
  est_cure$BIC[i] <- (length(fsr_use) + 1)*log(fsr_fits_[[i]]$N)  + 2*opt_obj_[[i]]$value
}

## print tables with AIC, cure proportion, upper and lower confidence intervals 