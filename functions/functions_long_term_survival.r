#/*******************************************************************************
#* XXXXXXXX - (short program description)
#********************************************************************************
#* Project : cd66600a.pbe
#* Study : m66601g.pbe
#* Author : Felizzi, Federico {MGPM~Basel} 
#* SAS version : 9.2
#* $Id: sasheader.sas,v 1.6 2014-01-31 16:39:43+01 ravenhas Exp $
#* $Source: /opt/BIOSTAT/razordb/RAZOR_UNIVERSE/DOMAIN_01/acp/Archive/RZ_VCS/acp/sas/sasheader.sas,v $
#********************************************************************************
#* (full program description)
#*
#*
#*
#*
#********************************************************************************
#* Inputs :
#*
#*
#* Outputs :
#*
#*
#* Called macros :
#*
#*
#----- for macros only -------
#********************************************************************************
#* Positional parameters :
#*
#* Named parameters :
#*
#----- end for macros only -----
#*******************************************************************************/

  #series of functions that will be used by the fitting program 
  
  #hazard and survival functions
  
  rate_c_years <- function(years, table){
    
    yr_l <- floor(years)
    yr_h <- ceiling(years)
    if ( yr_l != yr_h){
      s_h <- table$qx[yr_h + 2]
      s_l <- table$qx[yr_h+1]
      value <- (s_h - s_l)*(years - yr_l) + s_l
    }else{
      value <- table$qx[yr_h+1]
    }
    if (yr_h > 109){ value <- 0 }
    return(value)
  }



rate_cv_years <- function(years, table){
  n <- length(years);
  v_ret <- rep(0, length(years));
  for (i in 1:n){
    v_ret[i] <- rate_c_years(years[i], table = table);
  }
  return(v_ret)
}


Surv_c_years <- function(years, table){
  
  yr_l <- floor(years)
  yr_h <- ceiling(years)
  if ( yr_l != yr_h){
    s_h <- table$lx[yr_h+1]/1e5
    s_l <- table$lx[yr_h+0]/1e5
    value <- (s_h - s_l)*(years - yr_l) + s_l
  }else{
    value <- table$lx[yr_h+1]/1e5
  }
  if (yr_h > 109){ value <- 0 }
  return(value)
}



Surv_cv_years <- function(years, table){
  n <- length(years);
  v_ret <- rep(0, length(years));
  for (i in 1:n){
    v_ret[i] <- Surv_c_years(years[i], table = table);
  }
  return(v_ret)
}
