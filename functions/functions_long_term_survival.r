################################
# Long-term survival functions #
################################

################################
# Functions to calculate the mortality hazard rate and survival functions,
#  using a lifetable-style input table by age (expressed in years).
# Both for rates and survival, there is
#  - a function to calculate the rate/survival value
#  - a function to provide calculated rate/survival values as a vector
################################

# Hazard functions

## Calculate mortality hazard rate 
rate_c_years <- function(years, table) {
  # Prepare qx column from table
  table$qx <- as.numeric(table$qx)
  # Interpolate rates where necessary
  yr_l <- floor(years)
  yr_h <- ceiling(years)
  
  if (yr_l != yr_h) {
    s_h <- table$qx[yr_h + 2] # starting year is 0
    s_l <- table$qx[yr_h + 1]
    value <- (s_h - s_l) * (years - yr_l) + s_l
  } else {
    value <- table$qx[yr_h + 1]
  }
  # For very high ages, set mortality hazard rate to 0
  if (yr_h > 109) {
    value <- 0
  }
  return(value)
}

## Create vector of mortality hazard rates
rate_cv_years <- function(years, table) {
  n <- length(years)
  v_ret <- rep(0, length(years))
  
  for (i in 1:n) {
    v_ret[i] <- rate_c_years(years[i], table = table)
  }
  return(v_ret)
}

# Survival functions

## Calculate survival
Surv_c_years <- function(years, table) {
  # Interpolate survival where necessary
  yr_l <- floor(years)
  yr_h <- ceiling(years)
  
  if (yr_l != yr_h) {
    s_h <- table$lx[yr_h + 1] / 1e5
    s_l <- table$lx[yr_h + 0] / 1e5
    value <- (s_h - s_l) * (years - yr_l) + s_l
  } else {
    value <- table$lx[yr_h + 1] / 1e5
  }
  # For very high ages, set mortality hazard rate to 0
  if (yr_h > 109) {
    value <- 0
  }
  # Print result to R console
  print(paste0("value = ", value))
  return(value)
}

## Create vector of survival probabilities
Surv_cv_years <- function(years, table) {
  n <- length(years)
  v_ret <- rep(0, length(years))
  # Print initial message to R console
  print("Surv cv invoked")
  
  for (i in 1:n) {
    v_ret[i] <- Surv_c_years(years[i], table = table)
  }
  return(v_ret)
}
