##########################################################################
# Worked example of an analysis where the cure fraction is estimable     #
# from the trial data, using simulated BRIM3 data                        #
##########################################################################

###############################################################################
# In this script, anonymized and masked data from the BRIM3 trial, which
#  compared dacarbazine and vemurafenib in the treatment of
#  BRAFV600 mutation-positive metastatic melanoma in (Chapman et al., 2011,
#  N Engl J Med), were used to provide a
#  worked example of how to estimate cure fractions from a clinical trial
###############################################################################

# Load packages and custom functions --------------------------------------
require(flexsurv)
require(MASS)

source("functions/funs_hazard.R")
source("functions/funs_likelihood.R")

# Load trial data ---------------------------------------------------------
brim3 <- read.csv("data/trials/brim3_simulated.csv")


# Define models and data containers ---------------------------------------

# Parametric models to be considered in estimation ====
models <- c(
  "exponential", "weibull", "llogis", "lognormal", "gompertz",
  "gamma", "gengamma"
)
N_m <- length(models)


# Initialize containers for modelling outputs ====
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

est_cure <- data.frame(
  function_type = character(N_m), cure_rate = numeric(N_m),
  stderr = numeric(N_m), lower.ci = numeric(N_m),
  upper.ci = numeric(N_m), AIC = numeric(N_m),
  BIC = numeric(N_m), stringsAsFactors = FALSE
)

# Cure fraction values for later comparison ===
cure_vals <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2)
N_c <- length(cure_vals)

# Initialize containers for comparison of cure fraction values
parms_fit_frame <- list()
init_pos <- list()

AIC_cures <- data.frame(matrix(ncol = N_c + 1, nrow = N_c))
colnames(AIC_cures) <- c("Model", paste0("Cure_", cure_vals * 100))

AREA_cures <- data.frame(matrix(ncol = N_c * 2 + 3, nrow = N_c))
colnames(AREA_cures) <- c("Model",
                          paste0("Cure_", cure_vals * 100),
                          paste0("Cure_", cure_vals * 100, "_res"),
                          "No_cure", "No_cure_res")

# Function to calculate area under survival curve
area_surv <- function(x, y) {
  step <- x[2] - x[1]
  len <- length(y)
  val <- sum(y[2:(len - 1)])
  val <- val + 0.5 * (y[1] + y[len])
  val <- step * val
  
  return(val)
}

# Start analysis ----------------------------------------------------------

# Estimate the Kaplan-Meier curve ====
KM <- survfit(Surv(as.numeric(TIME), CNSR) ~ 1, data = brim3)

# Estimate the various models specified above
for (i in 1:N_m) {# CHECK THIS
  # Calculate model fits
  print(models[i])
  print("pre fit")
  fsr_fits_[[i]] <- flexsurvreg(Surv(as.numeric(TIME), CNSR) ~ 1,
                                data = brim3, dist = models[i])
  print("post fit")
  # fsr_use <- fsr_fits_[[i]]$res[, "est"] #CHECK THIS
  # len <- length(fsr_use)
  # 
  # Run optimization
  opt_obj_[[i]] <- try(optim(
    par = c(fsr_use, 0.23), # dummy starting parameter for the cure fraction
    ll.mix, # 
    table = brim3,
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
  
  # Store results in predefined containers
  est_cure$function_type[i] <- models[i]
  est_cure$cure_rate[i]     <- opt_obj_[[i]]$par[length(fsr_use) + 1]
  est_cure$stderr[i]        <- sqrt(invHuse[length(fsr_use) + 1,
                                            length(fsr_use) + 1])
  est_cure$lower.ci[i]      <- est_cure$cure_rate[i] - 
    qt(0.975, fsr_fits_[[i]]$N) * est_cure$stderr[i]
  est_cure$upper.ci[i]      <- est_cure$cure_rate[i] + 
    qt(0.975, fsr_fits_[[i]]$N) * est_cure$stderr[i]
  est_cure$AIC[i]           <- 2 * (length(fsr_use) + 1) + 
    2 * opt_obj_[[i]]$value
  est_cure$BIC[i]           <- (length(fsr_use) + 1) *
    log(fsr_fits_[[i]]$N) + 2 * opt_obj_[[i]]$value
  
  St_     <- do.call(fsr_fits_[[i]]$dfns$p,
                 args = c(as.list(opt_obj_[[i]]$par[1:(len)]),
                          list(q = xty, lower.tail = FALSE)))
  St_base <- do.call(fsr_fits_[[i]]$dfns$p,
                     args = c(as.list(fsr_fits_[[i]]$res[, "est"]),
                              list(q = xty, lower.tail = FALSE)))
  y_      <- surv_mean * (est_cure$cure_rate[i] +
                            (1 - est_cure$cure_rate[i]) * St_)
  
  St_models_cures[[i]]      <- St_
  y_models_cures[[i]]       <- y_
  St_base_models_cures[[i]] <- St_base
}

# Calculate the hazard rates and mean survival ====
hazard_out <- hazard_time(
  table_ = brim3,
  evttme = "TIME",
  sex = "SEX",
  age = "AGE",
  year = "YEAR",
  country_trial = "COUNTRY",
  country_output = NULL
)

brim3$rate_mod <- hazard_out$rate_vec
surv_mean <- colSums(hazard_out$surv_mat, na.rm = T) /
  length(which(!is.na(hazard_out$surv_mat[, 1])))
xty <- hazard_out$xty # Modeling time horizon


# Print outputs ====
write.csv(est_cure, file = "reports/est_cures_BRIM3.csv")

png("reports/cure_parametric_fits_BRIM3.png")

# Kaplan-Meier curve
plot(KM, conf.int = F, lwd = 5, xlab = "Time (years)",
     ylab = "Survival proportion", xlim = c(0, 10), lty = 2,
     cex.axis = 1.5, cex.lab = 1.5)

# Parametric models
for (i in 1:7) {
  lines(xty, y_models_cures[[i]], col = palette()[i + 1], lwd = 2)
}
lines(KM, conf.int = F, lwd = 5, xlab = "Time (years)", lty = 2)
legend("topright", c("KM", models[1:7]), col = c("black", palette()[2:8]),
       lwd = 4, lty = c(2, rep(1, 7)), cex = 1.5)
dev.off()


# Compare parametric models with cure fractions ---------------------------

## For each cure fraction value, initialise data frame ====
for (j in 1:N_c) {
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

## Calculate fits for each parametric model and cure fraction value
for (i in 1:N_m) {
  png(paste("reports/plot_models_", models[i], "_BRIM3.png", sep = ""),
      width = 9, height = 6, units = "in", res = 500)

  fsr_fits_ <- flexsurvreg(Surv(as.numeric(TIME), CNSR) ~ 1,
                           data = brim3, dist = models[i])
  surv_fits_ <- survfit(Surv(as.numeric(TIME), CNSR) ~ 1,
                        data = brim3)
  opt_obj_f_ <- list()

  for (j in 1:N_c) {
    fsr_use <- fsr_fits_$res[, "est"]
    opt_obj_f_[[j]] <- try(optim(
      par = c(fsr_use), # dummy staring parameter for the cure-fraction
      ll.mix_f,
      pi_ = cure_vals[j],
      table = brim3,
      time_col = "TIME",
      cen_col = "CNSR",
      rate_col = "rate_mod",
      method = "BFGS",
      obj = fsr_fits_,
      hessian = TRUE
    ),
    silent = FALSE
    )
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
  cure_strs <- paste0("Cure_", cure_vals * 100)

  colors__ <- rainbow(N_c)

  plot(fsr_fits_, xlim = c(0, 25), main = "   ",
       xlab = "Time (years)", ylab = "Survival proportion",
       est = FALSE, ci = FALSE)


  fits_f1 <- data.frame(matrix(ncol = N_c + 3, nrow = lx))
  colnames(fits_f1) <- c("years", "weeks", cure_strs, "background")

  fits_f1[, "weeks"] <- seq(0, lx - 1, by = 1)
  fits_f1[, "years"] <- xty
  fits_f1[, "background"] <- surv_mean

  v_leg <- c()

  for (j in 1:N_c) {
    obj_ <- opt_obj_f_[[j]]
    len <- length(obj_$par)

    St_[[j]] <- do.call(fsr_fits_$dfns$p,
                        args = c(as.list(obj_$par[1:(len)]),
                                 list(q = xty, lower.tail = FALSE)))
    St_base_[[j]] <- do.call(fsr_fits_$dfns$p,
                             args = c(as.list(fsr_fits_$res[, "est"]),
                                      list(q = xty, lower.tail = FALSE)))
    # See Equation 1 in the paper
    y_[[j]] <- surv_mean * (cure_vals[j] + (1 - cure_vals[j]) * St_[[j]])

    fits_f1[, cure_strs[j]] <- St_[[j]]
    if ((1 == 1) & (j == 1)) {
      v_leg <- c(v_leg, paste("Parametric"))
      pl_exp <- St_base_[[j]]
    }
    lines(xty, pl_exp, col = "blue", lwd = 2.5, lty = 2)
    lines(xty, y_[[j]], col = colors__[j], lwd = 2.5)

    AIC_cures$Model[i] <- models[i]
    AREA_cures$Model[i] <- models[i]
    AIC_cures[i, cure_strs[j]] <- 2 * (length(len)) + 2 * obj_$value
    i10 <- which(xty <= 10)
    AREA_cures[i, cure_strs[j]] <- area_surv(xty, y_[[j]])
    AREA_cures[i, "No_cure"] <- area_surv(xty, St_base_[[j]])
    AREA_cures[i, paste(cure_strs[j],
                        "_res", sep = "")] <- area_surv(xty[i10], y_[[j]][i10])
    AREA_cures[i, "No_cure_res"] <- area_surv(xty[i10], St_base_[[j]][i10])
    v_leg <- c(v_leg, paste("Cure = ", signif(cure_vals[j] * 100, 2), "%"))
  }
  
  lines(fsr_fits_, est = FALSE, ci = FALSE)
  lines(surv_fits_, lwd = 2)

  legend("topright",
         c("KM and confidence interval", v_leg),
         col = c("black", "blue", colors__),
         lwd = c(1, rep(2.5, N_c + 1)),
         lty = c(1, 2, rep(1, N_c + 1)), cex = 1.05, ncol = 1)

  dev.off()
}
