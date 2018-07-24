##########################################################################
# Worked example of an analysis where the cure fraction is not estimable #
# from the trial data                                                    #
##########################################################################

###############################################################################
# In this script, anonymized and masked data from the coBRIM trial, which
#  compared cobimetinib+vemurafenib and placebo+vemurafenib in the treatment of
#  BRAFV600 mutation-positive unresectable stage IIIC or stage IV melanoma in
#  19 countries (Ascierto et al., 2016, Lancet Oncol), were used to provide a
#  worked example of how to implement cure models if the cure fraction needs to
#  be estimated using an external data source.
###############################################################################

# Load required packages and scripts with functions
library(flexsurv)
library(MASS)

source("functions/hazard_script.R")
source("functions/likelihood_funs.r")

# Load trial data
trial_data <- read.csv(paste0("libraries/analysis2_clinical_data.csv"))

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

surv_mean <- colSums(hazard_out$surv_mat, na.rm = TRUE) /
  length(which(!is.na(hazard_out$surv_mat[, 1])))
xty <- hazard_out$xty

## Specify the distributions of interest, to be used in parametric modelling
models <- c("exponential", "weibull", "llogis", "lognormal", "gompertz",
            "gamma", "gengamma")

## Create vector with estimates/guesses of cure fractions
cure_vals <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2)

## Initialize containers for various outputs
N_m <- length(models)
parms_fit_frame <- list()
init_pos <- list()

AIC_cures <- data.frame("Model" = character(N_m),
                        "Cure_0" = numeric(N_m), "Cure_1" = numeric(N_m),
                        "Cure_2" = numeric(N_m), "Cure_3" = numeric(N_m),
                        "Cure_4" = numeric(N_m), "Cure_5" = numeric(N_m),
                        "Cure_10" = numeric(N_m), "Cure_15" = numeric(N_m),
                        "Cure_20" = numeric(N_m), stringsAsFactors = FALSE)

area_cures <- data.frame("Model" = character(N_m),
                         "Cure_0" = numeric(N_m),
                         "Cure_1" = numeric(N_m),
                         "Cure_2" = numeric(N_m),
                         "Cure_3" = numeric(N_m),
                         "Cure_4" = numeric(N_m),
                         "Cure_5" = numeric(N_m),
                         "Cure_10" = numeric(N_m),
                         "Cure_15" = numeric(N_m),
                         "Cure_20" = numeric(N_m),
                         "Cure_0_res" = numeric(N_m),
                         "Cure_1_res" = numeric(N_m),
                         "Cure_2_res" = numeric(N_m),
                         "Cure_3_res" = numeric(N_m),
                         "Cure_4_res" = numeric(N_m),
                         "Cure_5_res" = numeric(N_m),
                         "Cure_10_res" = numeric(N_m),
                         "Cure_15_res" = numeric(N_m),
                         "Cure_20_res" = numeric(N_m),
                         "No_cure" = numeric(N_m),
                         "No_cure_res" = numeric(N_m),
                         stringsAsFactors = FALSE)
area_surv <- function(x, y) {
  step <- x[2] - x[1]
  len  <- length(y)
  val  <- sum(y[2 : (len - 1)])
  val  <- val + 0.5 * (y[1] + y[len])
  val  <- step * val
  return(val)
}

for (j in 1:length(cure_vals)) {
  parms_fit_frame[[j]] <- data.frame("Function" = character(0),
                                     "Param_Name" = character(0),
                                     "Value_mix" = numeric(0),
                                     "Cov_mix_p1" = numeric(0),
                                     "Cov_mix_p2" = numeric(0),
                                     "Cov_mix_p3" = numeric(0),
                                     "Cov_mix_PI" = numeric(0),
                                     stringsAsFactors = FALSE)
  init_pos[[j]] <- 1
}


## Fit and plot the various models defined above to the data
for (i in 1:length(models)) {
  # Open graphics device
  jpeg(paste("reports/plot_models_", models[i], ".jpg", sep = ""),
       width = 9, height = 6, units = "in", res = 500)

  # Fit the models, using the "external" cure fractions,
  #  using maximum likelihood
  fsr_fits_ <- flexsurvreg(Surv(as.numeric(TIME), CNSR) ~ 1, data = trial_data,
                           dist = models[i])
  surv_fits_ <- survfit(Surv(as.numeric(TIME), CNSR) ~ 1, data = trial_data)
  
  opt_obj_f_ <- list()
  for (j in 1:length(cure_vals)) {
    fsr_use <- fsr_fits_$res[, "est"]
    
    opt_obj_f_[[j]] <- try(optim(
      par = c(fsr_use), # dummy starting parameter for the cure-fraction
      ll.mix_f,
      pi_ = cure_vals[j],
      table = trial_data,
      time_col = "TIME",
      cen_col = "CNSR",
      rate_col = "rate_mod",
      method = "BFGS",
      obj = fsr_fits_,
      hessian = TRUE
    ),
    silent = FALSE)
    
    if (1) {
      print(init_pos[[j]])

      invHuse_ <- ginv(opt_obj_f_[[j]]$hessian)
      parm_names <- rownames(fsr_fits_$res)
      print(parm_names)
      row_util <- length(parm_names)
      miu <- models[i]
      
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]] + row_util),
                           "Function"] <- miu
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]] + row_util - 1),
                           "Param_Name"] <- parm_names
      parms_fit_frame[[j]][(init_pos[[j]] + row_util),
                           "Param_Name"] <- "PI"
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]] + row_util - 1),
                           "Value_mix"] <- opt_obj_f_[[j]]$par
      parms_fit_frame[[j]][(init_pos[[j]] + row_util),
                           "Value_mix"] <- 0
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]] + row_util - 1),
                           4:(4 + row_util - 1)] <- invHuse_
      init_pos[[j]] <- init_pos[[j]] + row_util + 1
    }
  }


  lx <- length(xty)
  St_ <- list()
  St_base_ <- list()
  y_ <- list()

  # Plot survival estimates
  colors__ <- rainbow(length(cure_vals))
  plot(fsr_fits_, xlim = c(0, 25), main = "   ",
       xlab = "Time (years)", ylab = "Survival", est = FALSE, ci = FALSE)

  fits_f1 <- data.frame("years" = numeric(lx), "weeks" = numeric(lx),
    "Cure_0" = numeric(lx), "Cure_1" = numeric(lx),
    "Cure_2" = numeric(lx), "Cure_3" = numeric(lx),
    "Cure_4" = numeric(lx), "Cure_5" = numeric(lx),
    "Cure_10" = numeric(lx), "Cure_15" = numeric(lx),
    "Cure_20" = numeric(lx), "background" = numeric(lx),
    stringsAsFactors = FALSE)

  fits_f1[, "weeks"] <- seq(0, lx - 1, by = 1)
  fits_f1[, "years"] <- xty
  fits_f1[, "background"] <- surv_mean

  cure_strs <- c("Cure_0", "Cure_1", "Cure_2", "Cure_3", "Cure_4", "Cure_5",
                 "Cure_10", "Cure_15", "Cure_20")
  v_leg <- c()

  # Test goodness-of-fit for each model type and cure fraction estimate,
  #  using AIC and area under the curve (AUC)
  for (j in 1:length(cure_vals)) {
    obj_ <- opt_obj_f_[[j]]
    len <- length(obj_$par)

    St_[[j]] <- do.call(fsr_fits_$dfns$p,
                        args = c(as.list(obj_$par[1:(len)]),
                                 list(q = xty, lower.tail = FALSE)))
    St_base_[[j]] <- do.call(fsr_fits_$dfns$p,
                             args = c(as.list(fsr_fits_$res[, "est"]),
                                      list(q = xty, lower.tail = FALSE)))
    y_[[j]] <- surv_mean * (cure_vals[j] + (1 - cure_vals[j]) * St_[[j]])

    fits_f1[, cure_strs[j]] <- St_[[j]]
    if ((1 == 1) & (j == 1)) {
      v_leg <- c(v_leg, paste("Parametric; AIC = ",
                              signif(AIC(fsr_fits_), 6), "; Mean (AUC) = ",
                              signif(area_surv(xty, St_base_[[j]]), 3)))
      pl_exp <- St_base_[[j]]
    }
    
    lines(xty, pl_exp, col = "blue", lwd = 2.5, lty = 2)
    lines(xty, y_[[j]], col = colors__[j], lwd = 2.5)

    AIC_cures$Model[i] <- models[i]
    area_cures$Model[i] <- models[i]
    AIC_cures[i, cure_strs[j]] <- 2 * (length(len)) + 2 * obj_$value
    
    i10 <- which(xty <= 10)
    area_cures[i, cure_strs[j]] <- area_surv(xty, y_[[j]])
    area_cures[i, "No_cure"]    <- area_surv(xty, St_base_[[j]])
    area_cures[i, paste(cure_strs[j],
                        "_res", sep = "")] <- area_surv(xty[i10], y_[[j]][i10])
    area_cures[i, "No_cure_res"] <- area_surv(xty[i10], St_base_[[j]][i10])
    v_leg <- c(v_leg, paste("Cure = ", signif(cure_vals[j] * 100, 2),
                            "%; ", "AIC = ",
                            signif(AIC_cures[i, cure_strs[[j]]], 6),
                            "; Mean (AUC) = ",
                            signif(area_cures[i, cure_strs[j]], 3))
    )
  }
  
  lines(fsr_fits_, est = FALSE, ci = FALSE)
  lines(surv_fits_, lwd = 2)

  legend("topright", c("KM", v_leg), col = c("black", "blue", colors__),
         lwd = c(1, rep(2.5, length(cure_vals) + 1)),
         lty = c(1, 2, rep(1, length(cure_vals) + 1)), cex = 1.05, ncol = 1)

  # Write to external file
  write.csv(fits_f1, file = paste("reports/cures_AIC.csv", sep = ""))
  
  # Close graphics device
  dev.off()
}
