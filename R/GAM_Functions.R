
  
  # 2024 Institute for Clinical Evaluative Sciences. All rights reserved.
  
  # TERMS OF USE:
  
  ##Not for distribution.## This code and data is provided to the user solely for its own non-commerical use by individuals and/or not-for-profit corporations. Users shall not distribute without express written permission from the Institute for Clinical Evaluative Sciences.
  
  ##Not-for-profit.## This code and data may not be used in connection with profit generating activities.
  
  ##No liability.## The Institute for Clinical Evaluative Sciences makes no warranty or representation regarding the fitness, quality or reliability of this code and data.
  
##No Support.## The Institute for Clinical Evaluative Sciences will not provide any technological, educational or informational support in connection with the use of this code and data.

##Warning.## By receiving this code and data, user accepts these terms, and uses the code and data, solely at its own risk.


################################################################################# #


rates_f <- function(data){
  
  ti1 <-1
  firstyr <- 1985
  lastyr <- 2019
  years <- (lastyr-firstyr)+1
  label <- 15*years
 
  observedrates <- with(data,
                    tapply(N, list(floor( A /ti1)*ti1+ti1/2,
                                   floor((P-firstyr)/ti1)*ti1+ti1/2+firstyr), sum ) /
                      tapply(Y,list(floor( A /ti1)*ti1+ti1/2,
                                    floor((P-firstyr)/ti1)*ti1+ti1/2+firstyr), sum ) *10^6)
  observedrates <- as.data.frame(observedrates)
  
  observedrates <- gather(observedrates, P, rate, "1985.5":"2019.5")
  observedrates$P <- rep(firstyr:lastyr, each = 15, times = 1)
  observedrates$A <- rep(0:14, each = 1, times = years)
  observedrates$status <- rep("observed", times = label)
  return(as.data.frame(observedrates))
}


Pfr_export_f <- function(data){
  
  a.pt <- seq(0, 14.9, 1/10)
  p.predict <- c(2020:2040)
  
  res.list <- vector(mode = "list", length = length(p.predict))
  for (i in 1:length(p.predict)) {
    res.list[[i]] <- data.frame(A=a.pt, P=p.predict[i] , Y=10^6)
  }
  return(do.call(rbind,res.list))
  
}

ci.pred2 <- function(obj, newdata, Exp = NULL){
  
  zz <- predict(obj, newdata = newdata, se.fit = TRUE, 
                type = "link")
  zz <- cbind(zz$fit, zz$se.fit)
  colnames(zz) <- c("lograte","logSE")
  return(zz)
}


gammodel_f <- function(data, modelno){
#Nov2021 results uses:
  # firstyr_ocr <- 1970
  # firstknot <- firstyr_ocr-5
  # firstyr <- 1985 #first year of observed data
  # firstyr_knot <- firstyr+1 #first year of observed data minus 1
  # lastyr <- firstyr+((nrow(data)/15)/2-1)-1 #last year of observed data minus 1
  # lastyr_pred <- 2040+5 #last knot is 5 years after the last predicted year
  # knots <- list(x = c(firstknot, firstyr_knot, lastyr, lastyr_pred))
  
#for predicting 2020-2040:
#knots <- list(x = c(1965, 1986, 2018, 2044)) #Nov 26 2021- predict 1970-2040 (knots at inner points of observed data and extreme end points (+5 years after first observed yr and 5+ years after last predicted year))

  #March 2022 results uses:
  firstknot <- 1950 #first knot is first predicted year
  firstyr <- 1985 #first year of observed data
  firstyr_knot <- firstyr-1 #first year of observed data minus 1
  lastyr <- firstyr+((nrow(data)/15)/2-1)+1 #last year of observed data plus 1 (to cover full range of observed data, plus a little)
  lastyr_pred <- 2040 #last knot is last predicted year
  knots <- list(x = c(firstknot, firstyr_knot, lastyr, lastyr_pred))
  
  
if(modelno == 1){
  m1 <- gam(N ~ s(A,k=10, bs = "bs", m = c(3,1)) + s(P,k=10, bs = "bs", m = c(3,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
} else if(modelno == 2){
  m1 <- gam( N ~ s(A,k=10, bs = "bs", m = c(3,2)) + s(P,k=10, bs = "bs", m = c(3,2)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
} else if(modelno == 3){
  m1 <- gam( N ~ s(A,k=10, bs = "bs", m = c(3,2,1)) + s(P,k=10, bs = "bs", m = c(3,2,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
} else if(modelno == 4){
  m1 <- gam(N ~ s(A,k=5, bs = "bs", m = c(3,1)) + s(P,k=5, bs = "bs", m = c(3,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
} else if(modelno == 5){
  m1 <- gam( N ~ s(A,k=5, bs = "bs", m = c(3,2)) + s(P,k=5, bs = "bs", m = c(3,2)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)  
} else if(modelno == 6){
  m1 <- gam( N ~ s(A,k=5, bs = "bs", m = c(3,2,1)) + s(P,k=5, bs = "bs", m = c(3,2,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
}

#gam.check(m1)
#summary(m1)
#ci.exp(m1)
#m1$coefficients
#vcov.gam(m1)
  
return(m1)

}

##### - use if running one gam model (with gammodel_f function, not *gammodels_f*)
predrates_f <- function(data){
  
  rates <- as.data.frame(ci.pred2(m1, Pfr_export))
  rates2 <- bind_cols(Pfr_export, rates)
  rates2$age <- rep(seq(0,14), each = 10, times = 20)
  rates2$rate <- exp(rates2$lograte)
  rates2$SE <- exp(rates2$logSE)
  rates2$lower <- exp(rates2$lograte-1.96*rates2$logSE)
  rates2$upper <- exp(rates2$lograte+1.96*rates2$logSE)
  
  rates_median <- rates2 %>% group_by(age, P) %>% arrange(desc(rate)) %>% slice(5) %>% ungroup()
 
  predrates <- rates_median
  
  predrates <- predrates[c(-1, -3)]
  predrates <- predrates %>% dplyr::rename(A = age)
  predrates$status <- rep("predicted", times = 300)
  return(as.data.frame(predrates))
}

plotrates_f <- function(data){
  plotrates <- dplyr::bind_rows(observedrates, predrates)
  plotrates$A <- as.factor(plotrates$A)
  return(as.data.frame(plotrates))
}


plot_f <- function(data){
  #old - from older gam pred functions file #plot <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(subtitle = "Projected cancer incidence rates for 2020-2040", x = "Year of diagnosis", y = "Rate per million person-years") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A) + theme(legend.position = "bottom") #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 1.0, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  plot <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 2020-2040") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A) + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 1.0, decimal.mark = '.'))
  return(plot)
}

plot_freey_f <- function(data){
  plot_freey <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 2020-2040, free y axes, tsCV using multi-year forecasting") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A, scales = "free_y") + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + theme(axis.text.y = element_blank())
  #plot_freey <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(subtitle = "Projected cancer incidence rates for 2020-2040", x = "Year of diagnosis", y = "Rate per million person-years") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A, scales = "free_y") + theme(legend.position = "bottom") #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
return(plot_freey)
}


plot_CI_f <- function(data){
  #plot <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(subtitle = "Projected cancer incidence rates for 2020-2040", x = "Year of diagnosis", y = "Rate per million person-years") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A) + theme(legend.position = "bottom") #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 1.0, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  #plot_CI <- plot + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)
  plot <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 2020-2040, 95% CIs") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A) + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 1.0, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  plot_CI <- plot + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)
  return(plot_CI)
}


plot_CI_freey_f <- function(data){
  #plot_freey <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(subtitle = "Projected cancer incidence rates for 2020-2040", x = "Year of diagnosis", y = "Rate per million person-years") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A, scales = "free_y") + theme(legend.position = "bottom") #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  #plot_CI_freey <- plot_freey + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)
  plot_freey <- ggplot(plotrates, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 2020-2040, free y axes and 95% CIs, tsCV using multi-year forecasting") + ggtitle(names(res.list[i])) + facet_wrap(plotrates$A, scales = "free_y") + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates$rate), max(plotrates$rate)), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  plot_CI_freey <- plot_freey + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)
  return(plot_CI_freey)
}


Pfr_export_all_f <- function(data){
  
  a.pt <- seq(0, 14.9, 1/10)
  #p.predict_all <- c(1970:2040)
  p.predict_all <- c(1950:2040)
  
  res.list <- vector(mode = "list", length = length(p.predict_all))
  for (i in 1:length(p.predict_all)) {
    res.list[[i]] <- data.frame(A=a.pt, P=p.predict_all[i] , Y=10^6)
  }
  return(do.call(rbind,res.list))
}

predrates_all_f <- function(data){

  rates_all <- as.data.frame(ci.pred2(m1, Pfr_export_all))
  rates2 <- bind_cols(Pfr_export_all, rates_all)
  rates2$age <- rep(seq(0,14), each = 10, times = length(p.predict_all))
  rates2$rate <- exp(rates2$lograte)
  rates2$SE <- exp(rates2$logSE)
  rates2$lower <- exp(rates2$lograte-1.96*rates2$logSE)
  rates2$upper <- exp(rates2$lograte+1.96*rates2$logSE)
  
  rates_median <- rates2 %>% group_by(age, P) %>% arrange(desc(rate)) %>% slice(5) %>% ungroup()
  
  predrates_all <- rates_median
  
  predrates_all <- predrates_all[c(-1, -3)]
  predrates_all <- predrates_all %>% dplyr::rename(A = age)
  predrates_all$status <- rep("predicted", times = 15*length(p.predict_all))
  return(as.data.frame(predrates_all))
}

plotrates_all_f <- function(data){
  plotrates_all <- dplyr::bind_rows(observedrates, predrates_all)
  plotrates_all$A <- as.factor(plotrates_all$A)
  return(as.data.frame(plotrates_all))
}


plot_all_f <- function(data){
  
  plot_all <- ggplot(plotrates_all, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 1950-2040") + ggtitle(names(res.list[i])) + facet_wrap(plotrates_all$A) + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates_all$rate), max(plotrates_all$rate)), labels = scales::number_format(accuracy = 1.0, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  return(print(plot_all))
}

plot_all_freey_f <- function(data){
  
  plot_all_freey <- ggplot(plotrates_all, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 1950-2040, free y axes, tsCV using multi-year forecasting") + ggtitle(names(res.list[i])) + facet_wrap(plotrates_all$A, scales = "free_y") + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates_all$rate), max(plotrates_all$rate)), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  return(print(plot_all_freey))
}


plot_CI_all_f <- function(data){
  
  plot_all <- ggplot(plotrates_all, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 1950-2040, 95% CIs") + ggtitle(names(res.list[i])) + facet_wrap(plotrates_all$A) + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates_all$rate), max(plotrates_all$rate)), labels = scales::number_format(accuracy = 1.0, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  plot_CI_all <- plot_all + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)
  return(print(plot_CI_all))
}


plot_CI_all_freey_f <- function(data){
  
  plot_all_freey <- ggplot(plotrates_all, aes(x=P, y=rate)) + geom_point(aes(color=status)) + labs(x = "Year of diagnosis", y = "Rate per million person-years", subtitle = "Projected cancer incidence rates for 1950-2040, free y axes and 95% CIs, tsCV using multi-year forecasting") + ggtitle(names(res.list[i])) + facet_wrap(plotrates_all$A, scales = "free_y") + theme(legend.position = "bottom") + theme(axis.text.y = element_blank()) #+ scale_y_continuous(breaks = c(min(plotrates_all$rate), max(plotrates_all$rate)), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) #+ theme(axis.text.y = element_blank())
  plot_CI_all_freey <- plot_all_freey + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1)
  return(print(plot_CI_all_freey))
}

#Time-series cross-validation - 1 year ahead predictions

#create training datasets with ids 
data_tr_f <- function(data){
  
  #ids  <- c(1:33) #get all unique combinations for .id. . when predicting 1988 forward
  ids <- c(1:21) #when predicting 2000 forward
  ids <- as.data.frame(ids)
  
  data_new <- tsibble(data, key = A, index = P)
  
  data_tr <- data_new %>%
    #stretch_tsibble(.init = 3, .step = 1) #when predicting 1988 forward
    stretch_tsibble(.init = 15, .step = 1) %>% #stretch training datasets to predict P+1 yrs forward (P = 1999)
    relocate(P, A, .id)
  
  return(data_tr)
}


##### for each of the 21 split training sets, calculate observed rates for the period covered
rates_tscv_f <- function(data){
  
  ti1 <-1
  firstyr <- 1985
  lastyr <- 1985+((nrow(data)/15)/2-1)
  years <- (lastyr-firstyr)+1
  
  rates <- with(data,
                tapply(N, list(floor( A /ti1)*ti1+ti1/2,
                               floor((P-firstyr)/ti1)*ti1+ti1/2+firstyr), sum ) /
                  tapply(Y,list(floor( A /ti1)*ti1+ti1/2,
                                floor((P-firstyr)/ti1)*ti1+ti1/2+firstyr), sum ) *10^6)
  rates <- as.data.frame(rates)
  
  rates_tscv <- rates %>% pivot_longer(everything())
  rates_tscv <- rates_tscv %>% dplyr::rename(P = name, rate = value)
  rates_tscv$A <- rep(0:14, each = years, times = 1)
  rates_tscv$P <- rep(firstyr:lastyr, each = 1, times = 15)
  rates_tscv$status <- rep("observed", times = nrow(rates_tscv))
  return(as.data.frame(rates_tscv))

}


#Time-series cross-validation - 5 years ahead predictions

Pfr_tscv5_f <- function(data){
  
  firstyr_tscv5 <- 1985
  a.pt_tscv5 <- seq(0, 14.9, 1/10)
  #p.predict_tscv5 <- ( (1985+((nrow(data)/15)/2)) : (1985+((nrow(data)/15)/2)+4)  )
  
  #res.list_tscv5 <- vector(mode = "list", length = length(p.predict_tscv5))
  #for (i in 1:length(p.predict_tscv5)) {
  #  res.list_tscv5[[i]] <- data.frame(A=a.pt_tscv5, P=p.predict_tscv5[i] , Y=10^6)
  #}
  #return(do.call(rbind,res.list_tscv5))
  
  p.predict_tscv5 <- firstyr_tscv5+((nrow(data)/15)/2)+4
  
  Pfr_tscv5 <- data.frame(A=a.pt_tscv5, P=p.predict_tscv5, Y=10^6)
  return(Pfr_tscv5)
  
}

#run MULTIPLE fitted gams for TSCV per each training dataset - predicting 5 years forward
gammodels_tscv5_f<- function(data){
  #used for Nov2021 results - predicting 1970-2040 using 1985-2019 data
  # firstyr_ocr <- 1970
  # firstknot <- firstyr_ocr-5
  # firstyr <- 1985
  # firstyr_knot <- firstyr+1 #first year of observed data (1985) - 1
  # lastyr <- firstyr+((nrow(data)/15)/2-1)-1 #last year of observed data - 1
  # lastyr_pred_tscv5 <- (firstyr+((nrow(data)/15)/2)+4)+5 #last knot is 5th year after the last predicted year
  # knots_tscv5 <- list(x = c(firstknot, firstyr_knot, lastyr, lastyr_pred_tscv5))

  #used for March2022 results - predicting 1950-2040 using 1985-2019 data, altering knot placements to cover all observed data, plus a little, and knots at extreme end points of prediction
  firstknot <- 1950
  firstyr <- 1985
  firstyr_knot <- firstyr-1 #first year of observed data (1985) - 1 = 1984
  lastyr <- firstyr+((nrow(data)/15)/2-1)+1 #last year of observed data - 1
  lastyr_pred_tscv5 <- (firstyr+((nrow(data)/15)/2)+4) #last knot is the last predicted year
  knots_tscv5 <- list(x = c(firstknot, firstyr_knot, lastyr, lastyr_pred_tscv5))
  
  m1_tscv5 <- gam(N ~ s(A,k=10, bs = "bs", m = c(3,1)) + s(P,k=10, bs = "bs", m = c(3,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots_tscv5)
  #b-splines, cubic splines with 1st derivative penalty (penalises deviations from a flat function)
  
  m2_tscv5 <- gam( N ~ s(A,k=10, bs = "bs", m = c(3,2)) + s(P,k=10, bs = "bs", m = c(3,2)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots_tscv5)
  #b-splines, cubic splines with 2nd derivative penalty (penalising the curvature of the spline) (conventional cubic spline)
  
  m3_tscv5 <- gam( N ~ s(A,k=10, bs = "bs", m = c(3,2,1)) + s(P,k=10, bs = "bs", m = c(3,2,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots_tscv5)
  #b-splines, cubic splines with 1st and 2nd derivative penalties
  
  m4_tscv5 <- gam(N ~ s(A,k=5, bs = "bs", m = c(3,1)) + s(P,k=5, bs = "bs", m = c(3,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots_tscv5)
  #repeating same spline specification as 3 gams above, but reducing number of knots to 5 and 10
  
  m5_tscv5 <- gam( N ~ s(A,k=5, bs = "bs", m = c(3,2)) + s(P,k=5, bs = "bs", m = c(3,2)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots_tscv5)
  
  m6_tscv5 <- gam( N ~ s(A,k=5, bs = "bs", m = c(3,2,1)) + s(P,k=5, bs = "bs", m = c(3,2,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots_tscv5)
  
  p_m1 <- as_tibble(ci.pred2(m1_tscv5, Pfr_tscv5)) %>%
    dplyr::rename(fit_k10_31 = lograte, se_k10_31 = logSE)
  p_m2 <- as_tibble(ci.pred2(m2_tscv5, Pfr_tscv5)) %>%
    dplyr::rename(fit_k10_32 = lograte, se_k10_32 = logSE)
  p_m3 <- as_tibble(ci.pred2(m3_tscv5, Pfr_tscv5)) %>%
    dplyr::rename(fit_k10_321 = lograte, se_k10_321 = logSE)
  p_m4 <- as_tibble(ci.pred2(m4_tscv5, Pfr_tscv5)) %>%
    dplyr::rename(fit_k5_31 = lograte, se_k5_31 = logSE)
  p_m5 <- as_tibble(ci.pred2(m5_tscv5, Pfr_tscv5)) %>%
    dplyr::rename(fit_k5_32 = lograte, se_k5_32 = logSE)
  p_m6 <- as_tibble(ci.pred2(m6_tscv5, Pfr_tscv5)) %>%
    dplyr::rename(fit_k5_321 = lograte, se_k5_321 = logSE)
  
  crit <- 1.96
  
  models <- bind_cols(Pfr_tscv5, p_m1, p_m2, p_m3, p_m4, p_m5, p_m6) %>% pivot_longer(fit_k10_31:se_k5_321, names_sep = '_', names_to = c('variable', 'spline', 'penalty')) %>% 
    pivot_wider(names_from = variable, values_from = value) %>% dplyr::rename(lograte = fit, logSE = se) %>% 
    mutate(rate = exp(lograte), SE = exp(logSE), lwr_ci = exp(lograte - (crit * logSE)), upr_ci = exp(lograte + (crit * logSE)))
  
  #models$age <- rep(seq(0,14), each = 60, times = 5) #times 5 because we're predicting 5 years, not 1
  models$age <- rep(seq(0,14), each = 60, times = 1) #omly need the 5th year of the 5 year predictions
  
  models <- models %>% group_by(age, P, spline, penalty) %>% arrange(desc(rate)) %>% slice(5) %>% ungroup()
  models <- models[c(-3)]
  models$status2 <- rep("tsCV predicted", times = nrow(models))
  
  converge <- c(m1_tscv5$converged, m2_tscv5$converged, m3_tscv5$converged, m4_tscv5$converged, m5_tscv5$converged, m6_tscv5$converged)
  converge <- as.data.frame(converge)
  converge$modelno <- c(1,2,3,4,5,6)
  
  AIC <- c(m1_tscv5$aic, m2_tscv5$aic, m3_tscv5$aic, m4_tscv5$aic, m5_tscv5$aic, m6_tscv5$aic)
  AIC <- round(AIC,3)
  AIC <- as.data.frame(AIC)
  AIC$modelno <- c(1,2,3,4,5,6)
  
  df_m1 <- filter(models, spline == "k10" & penalty == "31")
  df_m2 <- filter(models, spline == "k10" & penalty == "32")
  df_m3 <- filter(models, spline == "k10" & penalty == "321")
  df_m4 <- filter(models, spline == "k5" & penalty == "31")
  df_m5 <- filter(models, spline == "k5" & penalty == "32")
  df_m6 <- filter(models, spline == "k5" & penalty == "321")
  
  df_models_tscv5 <- list(df_m1, df_m2, df_m3, df_m4, df_m5, df_m6, converge, AIC)
  return(df_models_tscv5)
}

#run all 6 gam models and predict - to plot and compare all 6 models' predictions

gammodels_all_f<- function(data){
  #Nov2021
  # firstyr_ocr <- 1970 #hide when predicting 1985 onward
  # firstknot <- firstyr_ocr-5 #hide when predicting 1985 onward
  # firstyr <- 1985 #first year of observed data is either 1985 or 1970
  # #firstknot <- firstyr-5 #first knot is 5 years before first year of observed data #unhide when predicting 1985 forward
  # firstyr_knot <- firstyr+1 #first year of observed data minus 1
  # lastyr <- firstyr+((nrow(data)/15)/2-1)-1 #last year of observed data minus 1
  # lastyr_pred <- 2040+5 #last knot is 5 years after the last predicted year
  # knots <- list(x = c(firstknot, firstyr_knot, lastyr, lastyr_pred))
  
  #Mar2022
  firstknot <- 1950 
  firstyr <- 1985 #first year of observed data is either 1985 or 1970
  firstyr_knot <- firstyr-1 #first year of observed data minus 1
  lastyr <- firstyr+((nrow(data)/15)/2-1)+1 #last year of observed data minus 1
  lastyr_pred <- 2040 #last knot is 5 years after the last predicted year
  knots <- list(x = c(firstknot, firstyr_knot, lastyr, lastyr_pred))
  
  #wiggliness penalty covers the range of x in observed data. to extend the penalty beyond range of x, pass in set of end points over which knots will be defined by specifying the two extreme end points that enclose region for prediction and two interior knots that cover range of data, plus a little"
  
  m1_tscv5 <- gam(N ~ s(A,k=10, bs = "bs", m = c(3,1)) + s(P,k=10, bs = "bs", m = c(3,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
  #b-splines, cubic splines with 1st derivative penalty (penalises deviations from a flat function)
  
  m2_tscv5 <- gam( N ~ s(A,k=10, bs = "bs", m = c(3,2)) + s(P,k=10, bs = "bs", m = c(3,2)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
  #b-splines, cubic splines with 2nd derivative penalty (penalising the curvature of the spline) (conventional cubic spline)
  
  m3_tscv5 <- gam( N ~ s(A,k=10, bs = "bs", m = c(3,2,1)) + s(P,k=10, bs = "bs", m = c(3,2,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
  #b-splines, cubic splines with 1st and 2nd derivative penalties
  
  m4_tscv5 <- gam(N ~ s(A,k=5, bs = "bs", m = c(3,1)) + s(P,k=5, bs = "bs", m = c(3,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
  #repeating same spline specification as 3 gams above, but reducing number of knots to 5 and 10
  
  m5_tscv5 <- gam( N ~ s(A,k=5, bs = "bs", m = c(3,2)) + s(P,k=5, bs = "bs", m = c(3,2)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
  
  m6_tscv5 <- gam( N ~ s(A,k=5, bs = "bs", m = c(3,2,1)) + s(P,k=5, bs = "bs", m = c(3,2,1)) + offset(log(Y)), family=poisson, data=data, method = "REML", knots = knots)
  
  p_m1 <- as_tibble(ci.pred2(m1_tscv5, Pfr_export_all)) %>%
    dplyr::rename(fit_k10_31 = lograte, se_k10_31 = logSE)
  p_m2 <- as_tibble(ci.pred2(m2_tscv5, Pfr_export_all)) %>%
    dplyr::rename(fit_k10_32 = lograte, se_k10_32 = logSE)
  p_m3 <- as_tibble(ci.pred2(m3_tscv5, Pfr_export_all)) %>%
    dplyr::rename(fit_k10_321 = lograte, se_k10_321 = logSE)
  p_m4 <- as_tibble(ci.pred2(m4_tscv5, Pfr_export_all)) %>%
    dplyr::rename(fit_k5_31 = lograte, se_k5_31 = logSE)
  p_m5 <- as_tibble(ci.pred2(m5_tscv5, Pfr_export_all)) %>%
    dplyr::rename(fit_k5_32 = lograte, se_k5_32 = logSE)
  p_m6 <- as_tibble(ci.pred2(m6_tscv5, Pfr_export_all)) %>%
    dplyr::rename(fit_k5_321 = lograte, se_k5_321 = logSE)
  
  crit <- 1.96
  
  models_all <- bind_cols(Pfr_export_all, p_m1, p_m2, p_m3, p_m4, p_m5, p_m6) %>% pivot_longer(fit_k10_31:se_k5_321, names_sep = '_', names_to = c('variable', 'spline', 'penalty')) %>% 
    pivot_wider(names_from = variable, values_from = value) %>% dplyr::rename(lograte = fit, logSE = se) %>% 
    mutate(rate = exp(lograte), SE = exp(logSE), lwr_ci = exp(lograte - (crit * logSE)), upr_ci = exp(lograte + (crit * logSE)))
  
  models_all$age <- rep(seq(0,14), each = 60, times = length(unique(Pfr_export_all$P)))
  models_all <- models_all %>% group_by(age, P, spline, penalty) %>% arrange(desc(rate)) %>% slice(5) %>% ungroup()
  models_all <- models_all[c(-3)]
  models_all$status <- case_when(models_all$spline == "k10" & models_all$penalty == "31" ~ "predicted, k=10, penalty=1",
                                 models_all$spline == "k10" & models_all$penalty == "32" ~ "predicted, k=10, penalty=2",
                                 models_all$spline == "k10" & models_all$penalty == "321" ~ "predicted, k=10, penalty=2,1",
                                 models_all$spline == "k5" & models_all$penalty == "31" ~ "predicted, k=5, penalty=1",
                                 models_all$spline == "k5" & models_all$penalty == "32" ~ "predicted, k=5, penalty=2",
                                 models_all$spline == "k5" & models_all$penalty == "321" ~ "predicted, k=5, penalty=2,1")
  return(models_all)
}

