###################################################
###             MSM functions                   ###
###################################################


###################################################
###         Function: state.model               ###
###       Fits then returns models, AIC,        ###
###       BIC, and plot for each dist           ###
###################################################
###    * Modify to change distributions fit     ###
###   Note: remember AIC.rec/BIC.rec indexing   ###
###################################################

state.model <- function(data, covariates, title){
  # data :      dataset
  # covariates: the string of covariates to be used
  # title :     the title of the plots
  
  fmla <- as.formula(paste("Surv(time, status)~ ", covariates))
  # fmla <- as.formula(paste("Surv(Tstart, Tstop, status)~ ", covariates)) 
  
  AIC.rec <- array(NA, dim = length(models))
  BIC.rec <- array(NA, dim = length(models))
  
  # Exponential
  tryCatch({
    state.model.exp <- NA
    state.model.exp <- flexsurvreg(fmla, dist="exp", data = data) 
    print("exp") # state.model.exp$events to get number of events
    
    AIC.rec[1] <- AIC(state.model.exp)
    BIC.rec[1] <- BIC(state.model.exp)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".exp", ".png"), width = 850, height = 700, res = 150)
    plot.model.exp <- plot(state.model.exp, col = 1, main = title)
    legend("topright", cex = 0.7, "exp" , col = 1, bty="n")
  },
  error = function(e){
    print("exp")
    print(e$message)
    return(state.model.exp)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Weibull
  tryCatch({
    state.model.weib <- NA
    state.model.weib <- flexsurvreg(fmla, dist="weibull", data = data)
    print("weib")
    
    AIC.rec[2] <- AIC(state.model.weib)
    BIC.rec[2] <- BIC(state.model.weib)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".weib", ".png"), width = 850, height = 700, res = 150)
    plot.model.weib <- plot(state.model.weib, col = 1, main = title)
    legend("topright", cex = 0.7, "weib" , col = 1, bty="n")
  },
  error = function(e){
    print("weib")
    print(e$message)
    return(state.model.weib)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Gamma
  tryCatch({
    state.model.gamma <- NA
    state.model.gamma <- flexsurvreg(fmla, dist="gamma", data = data) 
    print("gamma")
    
    AIC.rec[3] <- AIC(state.model.gamma)
    BIC.rec[3] <- BIC(state.model.gamma)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".gamma", ".png"), width = 850, height = 700, res = 150)
    plot.model.gamma <- plot(state.model.gamma, col = 4, main = title)
    legend("topright", cex = 0.7, "gamma" , col = 4, bty="n")
  },
  error=function(e){
    print("gamma")
    print(e$message)
    return(state.model.gamma) 
  })

  if (length(dev.list())!=0) {dev.off()}
    
  # LogNormal
  tryCatch({
    state.model.lnorm <- NA
    state.model.lnorm <- flexsurvreg(fmla, dist="lnorm", data = data)
    print("lnorm")
    
    AIC.rec[4] <- AIC(state.model.lnorm)
    BIC.rec[4] <- BIC(state.model.lnorm)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".lnorm", ".png"), width = 850, height = 700, res = 150)
    plot.model.lnorm <- plot(state.model.lnorm, col = 4, main = title)
    legend("topright", cex = 0.7, "lnorm" , col = 4, bty="n")
  },
  error=function(e){
    print("lnorm")
    print(e$message)
    return(state.model.lnorm)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Gompertz
  tryCatch({
    state.model.gompertz <- NA
    state.model.gompertz <- flexsurvreg(fmla, dist="gompertz", data = data)
    print("gompertz")
    
    AIC.rec[5] <- AIC(state.model.gompertz)
    BIC.rec[5] <- BIC(state.model.gompertz)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".gompertz", ".png"), width = 850, height = 700, res = 150)
    plot.model.gompertz <- plot(state.model.gompertz, col = 4, main = title)
    legend("topright", cex = 0.7, "gompertz" , col = 4, bty="n")
  },
  error=function(e){
    print("gompertz")
    print(e$message)
    return(state.model.gompertz)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # LogLogistic
  tryCatch({
    state.model.llogis <- NA
    state.model.llogis <- flexsurvreg(fmla, dist="llogis", data = data)   
    print("llogis")
    
    AIC.rec[6] <- AIC(state.model.llogis)
    BIC.rec[6] <- BIC(state.model.llogis)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".llogis", ".png"), width = 850, height = 700, res = 150)
    plot.model.llogis <- plot(state.model.llogis, col = 4, main = title)
    legend("topright", cex = 0.7, "llogis" , col = 4, bty="n")
  },
  error=function(e){
    print("llogis")
    print(e$message)
    return(state.model.llogis)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Spline
  tryCatch({
    state.model.spline <- NA
    state.model.spline <- flexsurvspline(fmla, k = 3, data = data)
    print("spline")
    
    AIC.rec[7] <- AIC(state.model.spline)
    BIC.rec[7] <- BIC(state.model.spline)
    
    knot <- state.model.spline$knots
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".spline", ".png"), width = 850, height = 700, res = 150)
    plot.model.spline <- plot(state.model.spline, col = 6, main = title)
    legend("topright", cex = 0.7, "spline" , col = 6, bty="n")
  },
  error=function(e){
    print("spline")
    print(e$message)
    return(state.model.spline)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Exponential Cure
  tryCatch({
    state.model.exp.c <- NA
    state.model.exp.c <- flexsurvcure(fmla, dist="exp", data = data, mixture = T) 
    print("exp cure")
    
    
    AIC.rec[8] <- AIC(state.model.exp.c)
    BIC.rec[8] <- BIC(state.model.exp.c)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".exp.c", ".png"), width = 850, height = 700, res = 150)
    plot.model.exp <- plot(state.model.exp.c, col = 1, main = title)
    legend("topright", cex = 0.7, "exp.c" , col = 1, bty="n")
  },
  error = function(e){
    print("exp cure")
    print(e$message)
    return(state.model.exp.c)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Weibull Cure
  tryCatch({
    state.model.weib.c <- NA
    state.model.weib.c <- flexsurvcure(fmla, dist="weibull", data = data, mixture = T) 
    print("weib cure")
    
    AIC.rec[9] <- AIC(state.model.weib.c)
    BIC.rec[9] <- BIC(state.model.weib.c)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".weib.c", ".png"), width = 850, height = 700, res = 150)
    plot.model.weib <- plot(state.model.weib.c, col = 1, main = title)
    legend("topright", cex = 0.7, "weib.c" , col = 1, bty="n")
  },
  error = function(e){
    print("weib cure")
    print(e$message)
    return(state.model.weib.c)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Gamma Cure
  tryCatch({
    state.model.gamma.c <- NA
    state.model.gamma.c <- flexsurvcure(fmla, dist="gamma", data = data, mixture = T) 
    print("gamma cure")
    
    AIC.rec[10] <- AIC(state.model.gamma.c)
    BIC.rec[10] <- BIC(state.model.gamma.c)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".gamma.c", ".png"), width = 850, height = 700, res = 150)
    plot.model.gamma <- plot(state.model.gamma.c, col = 1, main = title)
    legend("topright", cex = 0.7, "gamma.c" , col = 1, bty="n")
    if (length(dev.list())!=0) {dev.off()}
  },
  error = function(e){
    print("gamma cure")
    print(e$message)
    return(state.model.gamma.c)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # LogNormal Cure
  tryCatch({
    state.model.lnorm.c <- NA
    state.model.lnorm.c <- flexsurvcure(fmla, dist="lnorm", data = data, mixture = T) 
    print("lnorm cure")
    
    AIC.rec[11] <- AIC(state.model.lnorm.c)
    BIC.rec[11] <- BIC(state.model.lnorm.c)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".lnorm.c", ".png"), width = 850, height = 700, res = 150)
    plot.model.lnorm <- plot(state.model.lnorm.c, col = 1, main = title)
    legend("topright", cex = 0.7, "lnorm.c" , col = 1, bty="n")
  },
  error = function(e){
    print("lnorm cure")
    print(e$message)
    return(state.model.lnorm.c)
  })

  if (length(dev.list())!=0) {dev.off()}
    
  # Gompertz Cure
  tryCatch({
    state.model.gompertz.c <- NA
    state.model.gompertz.c <- flexsurvcure(fmla, dist="gompertz", data = data, mixture = T) 
    print("gompertz cure")
    
    AIC.rec[12] <- AIC(state.model.gompertz.c)
    BIC.rec[12] <- BIC(state.model.gompertz.c)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".gompertz.c", ".png"), width = 850, height = 700, res = 150)
    plot.model.gompertz <- plot(state.model.gompertz.c, col = 1, main = title)
    legend("topright", cex = 0.7, "gompertz.c" , col = 1, bty="n")
  },
  error = function(e){
    print("gompertz cure")
    print(e$message)
    return(state.model.gompertz.c)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # LogLogistic Cure
  tryCatch({
    state.model.llogis.c <- NA
    state.model.llogis.c <- flexsurvcure(fmla, dist="llogis", data = data, mixture = T) 
    print("llogis cure")
    
    AIC.rec[13] <- AIC(state.model.llogis.c)
    BIC.rec[13] <- BIC(state.model.llogis.c)
    
    png(file = paste0(output_file_directory, "/figs_msm/", title,".llogis.c", ".png"), width = 850, height = 700, res = 150)
    plot.model.llogis <- plot(state.model.llogis.c, col = 1, main = title)
    legend("topright", cex = 0.7, "llogis.c" , col = 1, bty="n")
  },
  error = function(e){
    print("llogis cure")
    print(e$message)
    return(state.model.lnorm.c)
  })
  
  if (length(dev.list())!=0) {dev.off()}
  
  # Store all results in res
  # if model is NA, store NA as the value otherwise take model output
  res = list(Exponential   = ifelse(is.na(state.model.exp),         NA, state.model.exp), 
             Weibull       = ifelse(is.na(state.model.weib),        NA, state.model.weib),
             Gamma         = ifelse(is.na(state.model.gamma),       NA, state.model.gamma),
             LogNormal     = ifelse(is.na(state.model.lnorm),       NA, state.model.lnorm),
             Gompertz      = ifelse(is.na(state.model.gompertz),    NA, state.model.gompertz),
             LogLogistic   = ifelse(is.na(state.model.llogis),      NA, state.model.llogis),
             Spline        = ifelse(is.na(state.model.spline),      NA, state.model.spline),
             Exponential.c = ifelse(is.na(state.model.exp.c),       NA, state.model.exp.c),
             Weibull.c     = ifelse(is.na(state.model.weib.c),      NA, state.model.weib.c),
             Gamma.c       = ifelse(is.na(state.model.gamma.c),     NA, state.model.gamma.c),
             LogNormal.c   = ifelse(is.na(state.model.lnorm.c),     NA, state.model.lnorm.c),
             Gompertz.c    = ifelse(is.na(state.model.gompertz.c),  NA, state.model.gompertz.c),
             LogLogistic.c = ifelse(is.na(state.model.llogis.c),    NA, state.model.llogis.c),
             AIC           = ifelse(is.na(AIC.rec), NA, AIC.rec),
             BIC           = ifelse(is.na(BIC.rec), NA, BIC.rec),
             knots         = ifelse(exists("knot"), list(knot), NA))
  return(res)
  print("end")
}

####################################################
###        Function: fun.model.list              ###
### uses transition matrix, and long format data ###
### to fit models specified in state.fit for all ###
### transitions and stores them in a list        ###
#################################################### 

fun.model.list <- function(v_n, m_t, df, cov.lists){
  # v_n = vector of state names
  # m_t = transition matrix
  # cov.lists = list of covariates for each transition
  # df = long format data frame
  if(length(cov.lists) < sum(!is.na(m_t))){message("cov.lists too short")
  } else if(length(cov.lists) > sum(!is.na(m_t))){message("cov.lists too long")}
  
  model.list <- list()
  
  for (i in 1:sum(!is.na(m_t))){  # for 1 to the number of transitions
    transition <- which(m_t == i) # index in the transition matrix where # transition = i
    
    # i is the transition number
    # transition is the index for transition i in the transition matrix
    s.from <- transition %% length(v_n)     # from state = row of the transition matrix, or remainder when divided by # states 
    s.to <- ceiling(transition/length(v_n)) # to state = column of the transition matrix, or rounded up division by # states
    
    print(paste("transition", i, "-", v_n[s.from], "to", v_n[s.to], sep = " "))
    
    assign(paste0("model_", v_n[s.from], "_", v_n[s.to]),              # assign "model_" s.from "_" s.to with state.model function output
           state.model(data = subset(df, to == s.to & from == s.from), # only include rows that transition from s.from and to s.to
                       covariates = cov.lists[i],                      # index list of covariate lists by transition number
                       title      = paste(v_n[s.from], "to", v_n[s.to]))) 
    model.list[[i]] <- eval(parse(text = paste0("model_", v_n[s.from], "_", v_n[s.to]))) # append output to model.list 
  }
  return(model.list)
}

########################################### 
#                                         #  
#        Function: fit.model              #                    
#         return AIC and BIC              #
#                                         # 
########################################### 

fit.model <- function (data, row.names){
  
  AIC <- as.data.frame(data$AIC, row.names) 
  BIC <- as.data.frame(data$BIC, row.names)
  
  AIC.min   <- sapply(AIC, function(x) if(length(min <- which.min(x))) min else NA)
  BIC.min   <- sapply(BIC, function(x) if(length(min <- which.min(x))) min else NA)
  
  data.name <- deparse(substitute(data))
  
  cov <- sapply(data[1:length(row.names)], function(x) ifelse(is.na(x["cov"]), 0, 1)) # check if cov matrix is NA 
  
  fit.res <- list(data.name, AIC, BIC, AIC.min, BIC.min, cov)
  
  return(fit.res)
}

####################################################
###          Function: fun.fit.list              ###
### uses list output by fun.model.list to create ###
###  a list of AIC and BIC values for all models ###
###         fitted for all transitions           ###
#################################################### 

fun.fit.list <- function(m_t, model.list, models){
  # m_t = transition matrix
  # model.list = model list fitted by fit.model.list
  # models = vector of model names 
  fit.list <- list()
  
  for (i in 1:sum(!is.na(m_t))){ # sum !is.na(m_t) is the number of transitions, from 1 to number of transitions 
    transition <- which(m_t == i) # index in the transition matrix where # transition = i
    s.from <- transition %% length(v_n)     # from state = row of the transition matrix, or remainder when divided by # states 
    s.to <- ceiling(transition/length(v_n)) # to state = column of the transition matrix, or rounded up division by # states
    
    fit.list[[i]] <- fit.model(data      = model.list[[i]],
                               row.names = models)
    # name the transition 
    fit.list[[i]][[1]] <- paste0("model_", v_n[s.from], "_", v_n[s.to])
  }
  return(fit.list)
}

###################################################
###         Function: models.AICBIC             ###
###       Creates a dataframe with the          ###
###        AIC and BIC of all models            ###
###           for all transitions               ###
###################################################

# Function to create a dataframe containing the AIC and BIC of all models

models.AICBIC <- function(models, tran.names, fit.list){
  # tran.names = transition names
  # fit.list = list of fitted models AIC/BIC from fit.model function for all transitions (fun.fit.list output)
  v_col_names <- c("transition", "model", "AIC", "BIC", "cov")
  a_model_sel <- array(data = NA,  # empty array to store transition, model, AIC, BIC
                       dim = c(length(models), length(v_col_names), length(fit.list)),
                       dimnames = list(1:length(models), # rows = models
                                       v_col_names,      # columns = transition, model, AIC, BIC, cov
                                       tran.names))      # 3rd dimension = transitions
  
  for (i in 1:length(fit.list)){
    a_model_sel[ , 1, i] <- fit.list[[i]][[1]]               #transition
    a_model_sel[ , 2, i] <- models                           #model
    a_model_sel[ , 3, i] <- c(fit.list[[i]][[2]]$`data$AIC`) #AIC
    a_model_sel[ , 4, i] <- c(fit.list[[i]][[3]]$`data$BIC`) #BIC
    a_model_sel[ , 5, i] <- c(fit.list[[i]][[6]])            #cov
  }
  
  # Store in a dataframe
  df_all_models <- plyr::adply(.data = a_model_sel, .margins = c(1,1,3)) %>% 
    mutate(AIC        = round(as.numeric(AIC), digits = 2),
           BIC        = round(as.numeric(BIC), digits = 2),
           transition = as.factor(transition)) %>% 
    dplyr::select(all_of(v_col_names)) 
  return(df_all_models)
}

########################################### 
#                                         #
#          Function: output.f             #  
#      Extract and organize output        # 
#         for each transition             #
#                                         #
###########################################
#      Modify if distributions change     #
###########################################

 output.f <- function(model){
   
   model.name <- deparse(substitute(model))
   
   tryCatch({
     model.exp <- model.exp.cov <- NA
    #model.exp <- model[["Exponential"]]$res
     model.exp.cov <- model[["Exponential"]]$cov
   },
   error=function(e){
     e
     print(paste("Exp NA"))
   })
   
   tryCatch({
     model.weib <- model.weib.cov <- NA
     #model.weib <- model[["Weibull"]]$res
     model.weib.cov <- model[["Weibull"]]$cov
   },
   error=function(e){
     e
     print(paste("Weib NA"))
   })
   
   tryCatch({
     model.gam <- model.gam.cov <- NA
     #model.gam <- model[["Gamma"]]$res
     model.gam.cov <- model[["Gamma"]]$cov
   },
   error=function(e){
     e
     print(paste("Gamma NA"))
   })
   
   tryCatch({
     model.lnorm <- model.lnorm.cov <- NA
     #model.lnorm <- model[["LogNormal"]]$res 
     model.lnorm.cov <- model[["LogNormal"]]$cov
   },
   error=function(e){
     e
     print(paste("lnorm NA"))
   })
   
   tryCatch({
     model.gompertz <- model.gompertz.cov <- NA
     #model.gompertz <- model[["Gompertz"]]$res 
     model.gompertz.cov <- model[["Gompertz"]]$cov
   },
   error=function(e){
     e
     print(paste("gompertz NA"))
   })
   
   tryCatch({
     model.llog <- model.llog.cov <- NA
     #model.llog <- model[["LogLogistic"]]$res 
     model.llog.cov <- model[["LogLogistic"]]$cov
   },
   error=function(e){
     e
     print(paste("llog NA"))
   })
   
   op <- list( model.name,
               "exp",
               #model.exp,
               model.exp.cov,
               "weib",
               #model.weib,
               model.weib.cov,
               "gam",
              # model.gam,
               model.gam.cov,
               "llog",
               #model.llog,
               model.llog.cov,
               "lnorm",
               #model.lnorm,
               model.lnorm.cov, 
               "gompertz",
               #model.gompertz,
               model.gompertz.cov)
   return(op)
 } 



