require(flexsurv)
require(MASS)

trial_data <- read.csv(paste0("libraries/example.csv"))
source("hazard_script.R")
source("functions/likelihood_funs.r")
hazard_out <- hazard_time(table_ = trial_data, 
                          evttme = "TIME" , 
                          sex = "SEX", 
                          age = "AGE", 
                          year = "YEAR", 
                          country_trial = "COUNTRY", 
                          country_output = NULL)

trial_data$rate_mod <- hazard_out$rate_vec

#surv_mean <- colSums(hazard_out$surv_mat, na.rm = T)/length(trial_data$AGE)
surv_mean <- colSums(hazard_out$surv_mat, na.rm = T)/487

models <- c("exponential", "weibull", "llogis", "lognormal", "gompertz", "gamma", "gengamma", "genf")

N_m <- length(models) 
xty <-seq(0,30, by = 0.2) ## timespan -- 30 years 
fsr_fits_ <- list()
opt_obj_ <- list()

St_models <- list()
y_models <- list()
St_base_models <- list()

St_models <- list()
y_models <- list()
St_base_models <- list()
St_models_cures <- list()
y_models_cures <- list() 
St_base_models_cures <- list()

est_cure <- data.frame(function_type = character(N_m), cure_rate = numeric(N_m), stderr = numeric(N_m), lower.ci = numeric(N_m), upper.ci = numeric(N_m), AIC = numeric(N_m), BIC = numeric(N_m), stringsAsFactors = FALSE)

#define the general purpose ll.mix function 


for (i in 1:length(models)){
  #perform the fits
  print(models[i])
  print("pre fit")
  fsr_fits_[[i]] <- flexsurvreg(Surv(as.numeric(TIME), CNSR)~1, data = trial_data, dist = models[i])
  print("post fit")
  fsr_use <- fsr_fits_[[i]]$res[,'est'];
  len <- length(fsr_use)
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
  print("optimization")
  invHuse <- ginv(opt_obj_[[i]]$hessian)
  est_cure$function_type[i] <- models[i]
  est_cure$cure_rate[i] <- opt_obj_[[i]]$par[length(fsr_use) + 1]
  est_cure$stderr[i] <-  sqrt(invHuse[length(fsr_use) + 1, length(fsr_use) + 1])
  est_cure$lower.ci[i] <- est_cure$cure_rate[i] - qt(0.975, fsr_fits_[[i]]$N)*est_cure$stderr[i]
  est_cure$upper.ci[i] <- est_cure$cure_rate[i] + qt(0.975, fsr_fits_[[i]]$N)*est_cure$stderr[i]
  est_cure$AIC[i] <- 2*(length(fsr_use) + 1) + 2*opt_obj_[[i]]$value
  est_cure$BIC[i] <- (length(fsr_use) + 1)*log(fsr_fits_[[i]]$N)  + 2*opt_obj_[[i]]$value
  
  St_ <- do.call(fsr_fits_[[i]]$dfns$p, args = c(as.list(opt_obj_[[i]]$par[1:(len)]),list(q = xty, lower.tail = FALSE)))
  print("pst St")
  St_base <- do.call(fsr_fits_[[i]]$dfns$p, args = c(as.list(fsr_fits_[[i]]$res[,'est']),   list(q = xty, lower.tail = FALSE)))
  y_ <- surv_mean*(est_cure$cure_rate[i] + (1 - est_cure$cure_rate[i])*St_)
  
  St_models_cures[[i]] <- St_
  y_models_cures[[i]] <- y_ 
  St_base_models_cures[[i]] <- St_base
}

## print tables with AIC, cure proportion, upper and lower confidence intervals 
write.csv(est_cure, file = "reports/est_cures.csv")


### create a plot fo the various fits 

pdf("reports/cure_parametric_fits.pdf")
plot(fsr_fits_[[1]], est =F, xlim = c(0,10), lwd.ci = 0, ci = F, lwd.obs = 5, xlab = "Time (years)")
for (i in 1:7){

  lines(xty, y_models_cures[[i]], col = palette()[i+1], lwd = 4)  
}
lines(fsr_fits_[[1]], est =F, lwd.ci = 0, ci = F, lwd.obs = 4)
legend("topright", models[1:7],col = palette()[2:8], lwd = 4 )
dev.off()

St_ <- do.call(fsr_fits_exp$dfns$p, args = c(as.list(obj_$par[1:(len)]),list(q = xty, lower.tail = FALSE)))
St_base <- do.call(fsr_fits_exp$dfns$p, args = c(as.list(fsr_fits_exp$res[,'est']),   list(q = xty, lower.tail = FALSE)))
y_ <- surv_mean_usa*(cure_vals[j] + (1 - cure_vals[j])*St_)

St_models_cures[[i]][[j]] <- St_
y_models_cures[[i]][[j]] <- y_ 
St_base_models_cures[[i]][[j]] <- St_base


## plotting the extrapolations 
plot()
