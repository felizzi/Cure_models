######################################################################
# Worked example of an analysis where the cure fraction is estimable #
# from the trial data                                                #
######################################################################

##############################################################################
# In this script, anonymized and masked data from the BRIM-3 trial, which
#  compared vemurafenib and dacarbazine in the treatment of BRAFV600 mutation-
#  positive melanoma in 12 countries (Chapman et al., 2011, N Engl J Med;
#  Chapman et al., 2017, Ann Oncol), were used to provide a worked example of
#  how to implement cure models if the cure fraction can be estimated from
#  the trial data.
##############################################################################

# Load required packages and scripts with functions
library(flexsurv)
library(MASS)

source("functions/hazard_script.R")
source("functions/likelihood_funs.r")

# Load trial data
trial_data <- read.csv(paste0("libraries/analysis1_clinical_data.csv"))

# Start the analysis

## Calculate the hazard rates and mean survival
hazard_out <- hazard_time(
  table_ = trial_data,
  evttme = "TIME",
  sex = "SEX",
  age = "AGE",
  year = "YEAR",
  country_trial = "COUNTRY",
  country_output = NULL
)

trial_data$rate_mod <- hazard_out$rate_vec
surv_mean <- colSums(hazard_out$surv_mat, na.rm = TRUE) / 487

## Specify the distributions of interest, to be used in parametric modelling
models <- c("exponential", "weibull", "llogis", "lognormal",
            "gompertz", "gamma", "gengamma", "genf")

## Initiate containers for various outputs
N_m <- length(models)
xty <- seq(from = 0, to = 30, by = 0.2) ## 30-year timespan

fsr_fits_ <- list()
opt_obj_ <- list()

St_base_models <- list()
St_base_models_cures <- list()

St_models <- list()
St_models_cures <- list()

y_models <- list()
y_models_cures <- list()

est_cure <- data.frame(function_type = character(N_m),
                       cure_rate = numeric(N_m), stderr = numeric(N_m),
                       lower.ci = numeric(N_m), upper.ci = numeric(N_m),
                       AIC = numeric(N_m), BIC = numeric(N_m),
                       stringsAsFactors = FALSE)

## Create Kaplan-Meier estimates
KM <- survfit(Surv(as.numeric(TIME), CNSR) ~ 1, data = trial_data)

## Fit the various models defined above to the data
for (i in 1:length(models)) {
  # Fit each model and estimate parameters, including the cure fraction, 
  #  using maximum likelihood
  print(models[i])
  print("pre fit")
  fsr_fits_[[i]] <- flexsurvreg(Surv(as.numeric(TIME), CNSR) ~ 1,
                                data = trial_data, dist = models[i])
  print("post fit")
  fsr_use <- fsr_fits_[[i]]$res[, "est"]
  len <- length(fsr_use)
  
  opt_obj_[[i]] <- try(optim(
    par = c(fsr_use, 0.30), # dummy starting parameter for the cure fraction
    ll.mix,
    table = trial_data,
    time_col = "TIME",
    cen_col = "CNSR",
    rate_col = "rate_mod",
    method = "L-BFGS-B",
    obj = fsr_fits_[[i]],
    hessian = TRUE,
    lower = c(rep(-Inf, length(fsr_use)), 0),
    upper = c(rep(Inf, length(fsr_use)), 1)
  ),
  silent = FALSE)
  print("optimization")
  invHuse <- ginv(opt_obj_[[i]]$hessian)
  
  # Store parameters for cure fraction and model goodness-of-fit
  ## Parametric distribution assumed
  est_cure$function_type[i] <- models[i]
  ## Estimated cure fraction and its standard error
  est_cure$cure_rate[i] <- opt_obj_[[i]]$par[length(fsr_use) + 1]
  est_cure$stderr[i] <- sqrt(invHuse[length(fsr_use) + 1, length(fsr_use) + 1])
  ## Lower and upper bound of the 95% confidence interval for the cure fraction
  est_cure$lower.ci[i] <- est_cure$cure_rate[i] -
    qt(0.975, fsr_fits_[[i]]$N) * est_cure$stderr[i]
  est_cure$upper.ci[i] <- est_cure$cure_rate[i] + 
    qt(0.975, fsr_fits_[[i]]$N) * est_cure$stderr[i]
  ## Goodness-of-fit measures: AIC and BIC
  est_cure$AIC[i] <- 2 * (length(fsr_use) + 1) + 2 * opt_obj_[[i]]$value
  est_cure$BIC[i] <- (length(fsr_use) + 1) * log(fsr_fits_[[i]]$N) + 
    2 * opt_obj_[[i]]$value
  
  # Obtain survival estimates
  St_ <- do.call(fsr_fits_[[i]]$dfns$p,
                 args = c(as.list(opt_obj_[[i]]$par[1:(len)]),
                          list(q = xty, lower.tail = FALSE))
  )
  print("pst St")
  St_base <- do.call(fsr_fits_[[i]]$dfns$p,
                     args = c(as.list(fsr_fits_[[i]]$res[, "est"]),
                              list(q = xty, lower.tail = FALSE))
  )
  y_ <- surv_mean * (est_cure$cure_rate[i] + (1 - est_cure$cure_rate[i]) * St_)
  
  St_base_models_cures[[i]] <- St_base
  St_models_cures[[i]] <- St_
  y_models_cures[[i]] <- y_
}

# Generate outputs from the analysis including plots

## Write results to external file
write.csv(est_cure, file = "reports/est_cures.csv")

## Plot the fits from the various models
pdf("reports/cure_parametric_fits.pdf")
### Kaplan-Meier curve
plot(KM, conf.int = FALSE, xlab = "Time (years)", ylab = "Proportion surviving",
     xlim = c(0, 10), lty = 2, lwd = 5)
### Parametric models
for (i in 1:7) {
  lines(xty, y_models_cures[[i]], col = palette()[i + 1], lwd = 2)
}
lines(KM, conf.int = FALSE, xlab = "Time (years)",
      ylab = "Proportion surviving",lty = 2, lwd = 5)
legend("topright", c(models[1:7], "KM"), col = c(palette()[2:8], "black"),
       lwd = 4)
dev.off()

## Plot the extrapolations
plot()
