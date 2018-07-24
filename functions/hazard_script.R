########################################
# Hazards and residual survival curves #
########################################

################################################################################
# Function to build hazard and residual survival curves for each subject, using
#  - Data frame with demographic and other required information
#  - Event times (in years)
#  - Gender (female [F] or male [M])
#  - Year during which the subject entered the study/follow-up period
################################################################################

# Clean and set up countries (background mortality is country-specific)

## Preparing the environment for country data
country_map <- new.env(hash = TRUE, parent = emptyenv())
c_mort_list <- read.csv("libraries/CountryList.csv")

### Store country names in environment prepared above
for (i in 1:length(c_mort_list$Code)) {
  country_use <- as.character(c_mort_list$Code3[i])
  country_map[[country_use]] <- as.character(c_mort_list$Country[i])
}
## For Germany, the Netherlands and Switzerland, country codes must be changed
country_map[["SUI"]] <- "Switzerland"
country_map[["NED"]] <- "Netherlands"
country_map[["GER"]] <- "Germany"

## Calculating the hazard, using country-specific life tables
source("functions/functions_long_term_survival.r")

hazard_time <- function(table_, evttme, sex, age, year, country_trial,
                        country_output) {
  # Load tables and loop over countries relevant for the trial
  Np <- length(table_[, sex]) # number of rows in the datatable
  rate_vec <- matrix(data = NA, nrow = Np, ncol = 1)
  surv_vec <- matrix(data = NA, nrow = Np, ncol = 1)

  xty <- seq(from = 0, to = 30, by = 0.2) # timespan -- 30 years
  surv_mat <- matrix(data = NA, nrow = Np, ncol = length(xty))

  Nxty <- length(xty)
  surv_mean_table <- data.frame(week = numeric(Nxty),
                                residual_surv = numeric(Nxty))

  # For each subject in the trial  
  for (i in 1:Np) {
    country_in_use <- as.character(table_[i, country_trial])
    year_in_use    <- table_[i, year]
    gender_in_use  <- table_[i, sex]

    ## Load background mortality tables
    ### For men
    if (gender_in_use == "MALE" | gender_in_use == "M") {
      table_in_use <- read.csv(paste0("libraries/",
                                      country_map[[country_in_use]],
                                      "/Male/Mortality.csv", sep = ""))
      #Select the latest year in case the year in the trial is missing
      max_year <- max(table_in_use$Year) 
      if (year_in_use > max_year) {
        year_in_use <- max_year
      }
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i]  <- rate_c_years(as.numeric(table_[i, age]) + 
                                     as.numeric(table_[i, evttme]),
                                   table_in_use)
      scvy_all <- Surv_cv_years(as.numeric(table_[i, age]) + xty, table_in_use)
      
      # Print messages if a year is not available
      if (length(is.na(scvy_all) > 0)) {
        print(paste(country_in_use, ": ", country_map[[country_in_use]], "; ",
                    gender_in_use, ": ", table_[i, evttme], "; ",
                    year_in_use, ": ", table_[i, age], " | ",
                    "year missing in table"))
      }

      if (is.na(rate_vec[i])) {
        if (length(table_in_use$Year) == 0) {
          print(paste(country_in_use, ": ", country_map[[country_in_use]], "; ",
                      gender_in_use, ": ", table_[i, evttme], "; ",
                      year_in_use, ": ", table_[i, age], " | ",
                      "year missing in table"))
        } else {
          print(paste(country_in_use, ": ", country_map[[country_in_use]], "; ",
                      gender_in_use, ": ", table_[i, evttme], "; ",
                      year_in_use, ": ", table_[i, age], " | ",
                      "year present in table"))
        }
      }
    }
    
    ### For women
    if (gender_in_use == "FEMALE" | gender_in_use == "F") {
      table_in_use <- read.csv(paste0("libraries/",
                                      country_map[[country_in_use]],
                                      "/Female/Mortality.csv", sep = ""))
      #Select the latest year in case the year in the trial is missing
      max_year <- max(table_in_use$Year) 
      if (year_in_use > max_year) {
        year_in_use <- max_year
      }
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i]  <- rate_c_years(as.numeric(table_[i, age]) + 
                                     as.numeric(table_[i, evttme]),
                                   table_in_use)
      scvy_all <- Surv_cv_years(as.numeric(table_[i, age]) + xty, table_in_use)
      
      # Print messages if a year is not available
      if (length(is.na(scvy_all) > 0)) {
        print(paste(country_in_use, ": ", country_map[[country_in_use]], "; ",
                    gender_in_use, ": ", table_[i, evttme], "; ",
                    year_in_use, ": ", table_[i, age], " | ",
                    "year missing in table"))
      }
      
      if (is.na(rate_vec[i])) {
        if (length(table_in_use$Year) == 0) {
          print(paste(country_in_use, ": ", country_map[[country_in_use]], "; ",
                      gender_in_use, ": ", table_[i, evttme], "; ",
                      year_in_use, ": ", table_[i, age], " | ",
                      "year missing in table"))
        } else {
          print(paste(country_in_use, ": ", country_map[[country_in_use]], "; ",
                      gender_in_use, ": ", table_[i, evttme], "; ",
                      year_in_use, ": ", table_[i, age], " | ",
                      "year present in table"))
        }
      }
    }
    
  # Calculate survival matrix  
  surv_mat[i, ] <- scvy_all / scvy_all[1]
  }

  # Object to return from function
  ret_obj <- list()
  
  ret_obj$rate_vec <- rate_vec
  ret_obj$surv_mat <- surv_mat

  surv_mean <- colSums(surv_mat) / length(table_[, sex]) # mean survival cure
  ret_obj$surv_mean <- surv_mean
  ret_obj$xty <- xty

  return(ret_obj)
}
