### script for the second example, namely with the addition of an external-input cure 

trial_data <- read.csv(paste0("libraries/example_1.csv"))
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

#surv_mean <- colSums(hazard_out$surv_mat, na.rm = T)/
surv_mean <- colSums(hazard_out$surv_mat, na.rm = T)/length(which(!is.na(hazard_out$surv_mat[,1])))

models <- c("exponential", "weibull", "llogis", "lognormal", "gompertz", "gamma", "gengamma")
cure_vals <- c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2)

parms_fit_frame <- list()

init_pos <- list()

AIC_cures <- data.frame("Model" = character(7), "Cure_0"= numeric(7) , "Cure_1" = numeric(7),  "Cure_2" = numeric(7), "Cure_3" = numeric(7), "Cure_4" = numeric(7), "Cure_5" = numeric(7), "Cure_10" = numeric(7), "Cure_15" = numeric(7), "Cure_20" = numeric(7),  stringsAsFactors=FALSE)
AREA_cures <- data.frame("Model" = character(7), "Cure_0"= numeric(7) , "Cure_1" = numeric(7),  "Cure_2" = numeric(7), "Cure_3" = numeric(7), "Cure_4" = numeric(7), "Cure_5" = numeric(7), "Cure_10" = numeric(7), "Cure_15" = numeric(7), "Cure_20" = numeric(7), "Cure_0_res"= numeric(7) , "Cure_1_res" = numeric(7),  "Cure_2_res" = numeric(7), "Cure_3_res" = numeric(7), "Cure_4_res" = numeric(7), "Cure_5_res" = numeric(7), "Cure_10_res" = numeric(7), "Cure_15_res" = numeric(7), "Cure_20_res" = numeric(7), "No_cure" = numeric(7),  "No_cure_res" = numeric(7),stringsAsFactors=FALSE)
area_surv <- function(x, y){
  step <- x[2] - x[1]
  len <- length(y)
  val <- sum(y[2:(len-1)])
  val <- val + 0.5*(y[1] + y[len])
  val <- step*val
  return(val)
}


for (j in 1:length(cure_vals)){
  parms_fit_frame[[j]] <- data.frame("Function" = character(0), "Param_Name" = character(0), "Value_mix" = numeric(0), "Cov_mix_p1" = numeric(0), "Cov_mix_p2" = numeric(0), "Cov_mix_p3" = numeric(0), "Cov_mix_PI" = numeric(0), stringsAsFactors=FALSE)
  init_pos[[j]] <- 1
}

for (i in 1:length(models)){
  jpeg(paste("reports/plot_models_",models[i],".jpg",sep = ""), width = 9, height = 6, units="in", res=500)
  
  fsr_fits_ <- flexsurvreg(Surv(as.numeric(TIME), CNSR)~1, data = trial_data, dist = models[i])
  surv_fits_ <- survfit(Surv(as.numeric(TIME), CNSR)~1, data = trial_data )
  opt_obj_f_ <- list()
  
  for (j in 1:length(cure_vals)){
    
    fsr_use <- fsr_fits_$res[,'est'];
    opt_obj_f_[[j]] <- try(optim(par=c(fsr_use), #dummy staring parameter for the cure-fraction 
                                          ll.mix_f, 
                                          pi_ = cure_vals[j],
                                          table = trial_data,
                                          time_col = "TIME",
                                          cen_col = "CNSR",
                                          rate_col = "rate_mod",
                                          method = "BFGS", 
                                          obj = fsr_fits_,
                                          hessian = TRUE) 
                                    , silent = FALSE)
    if (1){ 
      print(init_pos[[j]])
      
      invHuse_ <- ginv(opt_obj_f_[[j]]$hessian)
      parm_names <- rownames(fsr_fits_$res)
      print(parm_names)
      row_util <- length(parm_names)
      miu <- models[i]
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]]+row_util),"Function"] <- miu
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]]+row_util-1),"Param_Name"] <- parm_names
      parms_fit_frame[[j]][(init_pos[[j]]+row_util),"Param_Name"] <- "PI"
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]]+row_util-1),"Value_mix"] <- opt_obj_f_[[j]]$par
      parms_fit_frame[[j]][(init_pos[[j]]+row_util),"Value_mix"] <- 0
      parms_fit_frame[[j]][init_pos[[j]]:(init_pos[[j]]+row_util-1),4:(4+row_util - 1)] <- invHuse_;
      init_pos[[j]] <- init_pos[[j]] + row_util + 1
    }
    
    
    
  }
  
  
  lx <- length(xty)
  St_ <- list()
  St_base_ <- list()
  y_ <- list()
  
  colors__ <- rainbow(length(cure_vals))
  
  plot(fsr_fits_, xlim = c(0,25), main = "   ", xlab = "Time (years)", ylab = "Survival", est = FALSE, ci = FALSE)
  
  
  fits_f1 <- data.frame("years" = numeric(lx), "weeks" = numeric(lx), 
                              "Cure_0" = numeric(lx),
                              "Cure_1" = numeric(lx),
                              "Cure_2" = numeric(lx),
                              "Cure_3" = numeric(lx),
                              "Cure_4" = numeric(lx),
                              "Cure_5" = numeric(lx), 
                              "Cure_10"  = numeric(lx),
                              "Cure_15" = numeric(lx),
                              "Cure_20" = numeric(lx),
                              "background" = numeric(lx),
                              stringsAsFactors=FALSE)
  
  fits_f1[,"weeks"] <- seq(0,lx-1, by = 1)
  fits_f1[,"years"] <- xty
  fits_f1[,"background"] <- surv_mean
  
  cure_strs <- c("Cure_0", "Cure_1", "Cure_2", "Cure_3","Cure_4", "Cure_5", "Cure_10", "Cure_15", "Cure_20")
  v_leg <- c()
  
  for (j in 1:length(cure_vals)){
    obj_ <- opt_obj_f_[[j]]
    len <- length(obj_$par)
    
    St_[[j]] <- do.call(fsr_fits_$dfns$p, args = c(as.list(obj_$par[1:(len)]),list(q = xty, lower.tail = FALSE)))
    St_base_[[j]] <- do.call(fsr_fits_$dfns$p, args = c(as.list(fsr_fits_$res[,'est']),   list(q = xty, lower.tail = FALSE)))
    y_[[j]] <- surv_mean*(cure_vals[j] + (1-cure_vals[j])*St_[[j]])
    
    fits_f1[,cure_strs[j]] <- St_[[j]]
    if ((1 == 1) & (j ==1) ){
      v_leg <- c(v_leg,paste("Parametric; AIC = ", signif(AIC(fsr_fits_),6),"; Mean (AUC) = ", signif(area_surv(xty, St_base_[[j]]),3)))
      pl_exp <- St_base_[[j]]
    }
    lines(xty, pl_exp, col = "blue", lwd = 2.5, lty = 2 )
    lines(xty, y_[[j]], col = colors__[j], lwd = 2.5)
    
    AIC_cures$Model[i] <- models[i]
    AREA_cures$Model[i] <- models[i]
    AIC_cures[i,cure_strs[j]]<- 2*(length(len)) + 2*obj_$value
    i10 <- which(xty<=10)
    AREA_cures[i,cure_strs[j]] <- area_surv(xty, y_[[j]])
    AREA_cures[i,"No_cure"] <- area_surv(xty, St_base_[[j]])
    AREA_cures[i,paste(cure_strs[j],"_res",sep="")] <- area_surv(xty[i10], y_[[j]][i10])
    AREA_cures[i,"No_cure_res"] <- area_surv(xty[i10], St_base_[[j]][i10])
    v_leg <- c(v_leg,
               paste("Cure = ",signif(cure_vals[j]*100,2),"%; ","AIC = ",signif(AIC_cures[i, cure_strs[[j]]],6), "; Mean (AUC) = ", signif(AREA_cures[i,cure_strs[j]],3))
    )
    
  }
  lines(fsr_fits_,  est = FALSE, ci = FALSE)
  lines(surv_fits_, lwd = 2)
  
  legend("topright",c("KM", v_leg), col = c("black","blue",colors__), lwd = c(1,rep(2.5,length(cure_vals) +1 )), lty = c(1,2,rep(1,length(cure_vals) + 1)),cex = 1.05,ncol = 1)
  
  
  
  #createSheet(wb, name = models[i])
  #writing into sheets within an Excel workbook : 
  #writing ChickWeight data frame into chickSheet
  #write.csv(Atezo_fits_f1, file= paste(wd,"/",out_fold,"/","Fit_cure_atezo_",models[i],"_",subg,"_norway.csv",sep = ""))
  #writeWorksheet(wb, Atezo_fits_f, sheet = models[i], startRow = 1, startCol = 1)
  #saveWorkbook(wb)
  
  #dev.off()
}
