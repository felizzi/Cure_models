###########################################################################
# Function to build hazard and residual survival curves for each subject  #
###########################################################################

###############################################################################
# Inputs:
#  Data frame with clinical and demographic info
#  Event time column (in years)
#  Year in which the subject entered the trial
###############################################################################


# Load custom functions ---------------------------------------------------
source("functions/funs_long_term_survival.R")


# Read in and prepare background mortality for countries ------------------
country_map <- new.env(hash = TRUE, parent = emptyenv())

c_mort_list <- read.csv("data/CountryList.csv")
for (i in 1:length(c_mort_list$Code)) {
  country_use <- as.character(c_mort_list$Code3[i])
  country_map[[country_use]] <- as.character(c_mort_list$Country[i])
}
## Add additional keys to the map for NED SUI GER
country_map[["SUI"]] <- "Switzerland"
country_map[["NED"]] <- "Netherlands"
country_map[["GER"]] <- "Germany"

hazard_time <- function(table_, evttme, sex, age, year,
                        country_trial, country_output) {
  # Load in tables
  ## Loop over the countries in the trial
  N_p <- length(table_[, sex]) ## number of rows in the datatable

  rate_vec <- matrix(NA, nrow = N_p, ncol = 1)
  surv_vec <- matrix(NA, nrow = N_p, ncol = 1)

  xty <- seq(0, 30, by = 0.2) ## timespan -- 30 years
  surv_mat <- matrix(NA, nrow = N_p, ncol = length(xty))

  Nxty <- length(xty)
  surv_mean_table <- data.frame(week = numeric(Nxty),
                                residual_surv = numeric(Nxty))

  ## Loop over subjects in the trial
  for (i in 1:N_p) {
    country_in_use <- as.character(table_[i, country_trial])
    year_in_use <- table_[i, year]
    gender_in_use <- table_[i, sex]

    ## Load mortality tables
    ### Men
    if (gender_in_use == "MALE" | gender_in_use == "M") {
      table_in_use <- read.csv(paste0("data/", country_map[[country_in_use]],
                                      "/Male/Mortality.csv", sep = ""))
      max_year <- max(table_in_use$Year) ## select max year if trial year missing
      if (year_in_use > max_year) {
        year_in_use <- max_year
      }
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i] <- rate_c_years(as.numeric(table_[i, age]) +
                                    as.numeric(table_[i, evttme]), table_in_use)
      scvy_all <- Surv_cv_years(as.numeric(table_[i, age]) +
                                  xty, table_in_use)
    }
    
    ### Women
    if (gender_in_use == "FEMALE" | gender_in_use == "F") {
      table_in_use <- read.csv(paste0("data/", country_map[[country_in_use]],
                                      "/Female/Mortality.csv", sep = ""))
      max_year <- max(table_in_use$Year)
      if (year_in_use > max_year) {
        year_in_use <- max_year
      }
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i] <- rate_c_years(as.numeric(table_[i, age]) +
                                    as.numeric(table_[i, evttme]), table_in_use)
      scvy_all <- Surv_cv_years(as.numeric(table_[i, age]) + xty, table_in_use)
    }
    
    surv_mat[i, ] <- scvy_all / scvy_all[1]
  
  ## Create return object with hazard and survival rates
  ret_obj <- list()
  
  ret_obj$rate_vec <- rate_vec
  ret_obj$surv_mat <- surv_mat

  surv_mean <- colSums(surv_mat) / length(table_[, sex])
  ret_obj$surv_mean <- surv_mean
  ret_obj$xty <- xty

  return(ret_obj)
  }
}
