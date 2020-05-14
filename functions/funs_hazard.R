#### function to build the hazard for the subject and the residual survival curves for each subject
#### it takes as input 
#### the data frame with the info
#### the event time colum --- in years 
#### the colum with the gender of the subject M/F 
#### the year  in which the subject enetered the 
country_map <- new.env(hash=T, parent=emptyenv())
c_mort_list <- read.csv("data/CountryList.csv")
for (i in 1:length(c_mort_list$Code)){
  country_use <- as.character(c_mort_list$Code3[i])
  country_map[[country_use]] <- as.character(c_mort_list$Country[i])
}
## add additional keys to the map NED SUI GER
country_map[["SUI"]] <- "Switzerland"
country_map[["NED"]] <- "Netherlands"
country_map[["GER"]] <- "Germany"
source("functions/funs_long_term_survival.R")
hazard_time <- function(table_, evttme, sex, age, year, country_trial, country_output){
  # load in the file the tables 
  ## loop over the countries in the trial 
  Np <- length(table_[,sex]) ## number of rows in the datatable
  
  rate_vec <- matrix(,nrow = Np, ncol = 1);
  surv_vec <- matrix(,nrow = Np, ncol = 1);
  
  xty <-seq(0,30, by = 0.2) ## timespan -- 30 years 
  surv_mat <- matrix(,nrow = Np, ncol = length(xty));
  
  Nxty <- length(xty)
  #surv_mean_table_all <- data.frame(week = numeric(Nxty), residual_surv = numeric(Nxty))
  surv_mean_table <- data.frame(week = numeric(Nxty), residual_surv = numeric(Nxty))
  
  for (i in 1:Np){ ### loop over the subjects in the trial 
    ## manual adaptations of the country-codes
    country_in_use <- as.character(table_[i,country_trial])
    year_in_use   <- table_[i, year]
    gender_in_use <- table_[i,sex]
    
    
    ## load the tables 
    if (gender_in_use == "MALE" | gender_in_use == "M"){
      table_in_use <- read.csv(paste0("data/",country_map[[country_in_use]],"/Male/Mortality.csv", sep = ""))
      max_year <- max(table_in_use$Year) ## selection of the max year in case the year in the trial is missing
      if (year_in_use > max_year){ year_in_use <- max_year}
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme]), table_in_use);
      scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use);

    }
    if (gender_in_use == "FEMALE" | gender_in_use == "F"){
      table_in_use <- read.csv(paste0("data/",country_map[[country_in_use]],"/Female/Mortality.csv", sep = ""))
      max_year <- max(table_in_use$Year) ## selection of the max year in case the year in the trial is missing
      if (year_in_use > max_year){ year_in_use <- max_year}
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme]), table_in_use);
      
      scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use);
          }
    surv_mat[i,] <- scvy_all/scvy_all[1];
    #  print("MALE")
      #rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/365.25, table_in_use_M);
  #    rate_vec_USA[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_male_2013)
   #   surv_vec_USA[i] <- Surv_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_male_2013)/Surv_c_years(as.numeric(table_[i,age]), USA_table_male_2013)
      #scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_M);
  #    scvy_usa <- Surv_cv_years( as.numeric(table_[i,age]) + xty, USA_table_male_2013)
      #scvy_uk <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_UK_M);
      
   # }else{
  #    print("FEMALE")
      #rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/365.25, table_in_use_F)
  #    rate_vec_USA[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_female_2013)
   #   surv_vec_USA[i] <- Surv_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_female_2013)/Surv_c_years(as.numeric(table_[i,age]), USA_table_female_2013)
      #scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_F)
  #    scvy_usa <- Surv_cv_years( as.numeric(table_[i,age]) + xty, USA_table_female_2013)
      #scvy_uk <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_UK_F);
    }
    
    #surv_mat_all[i,] <- scvy_all/scvy_all[1];
   # surv_mat_usa[i,] <- scvy_usa/scvy_usa[1];
    #surv_mat_uk[i,] <- scvy_uk/scvy_uk[1];
  #}
  
  
  
  ###########################################################
  #mean survival curve 
  #surv_mean_all <- colSums(surv_mat_all)/Np
  #surv_mean_usa <- colSums(surv_mat_usa)/Np
  #surv_mean_uk <- colSums(surv_mat_uk)/Np
  #print("dimension of SURV_MAT UK")
  #print(dim(surv_mat_uk))
  
  ret_obj <- list()
  ret_obj$rate_vec <- rate_vec
  ret_obj$surv_mat <- surv_mat
  
  ###########################################################
  #mean survival curve 
  surv_mean <- colSums(surv_mat)/length(table_[,sex])
  ret_obj$surv_mean <- surv_mean
  ret_obj$xty <- xty
  #ret_obj$surv_mean_all <- surv_mean_all
  #ret_obj$surv_mean_uk <- surv_mean_uk
  #ret_obj$surv_mean_usa <- surv_mean_usa
  #ret_obj$time <- xty
  return(ret_obj)
}