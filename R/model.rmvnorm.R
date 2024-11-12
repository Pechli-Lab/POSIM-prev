#' Adding uncertainty with rmvnorm
#'
#' @param MSM_est list with selected model, covariates, estimates for each transition
#' @param MSM_cov list with covariance matrices for each transition
#' @param index index of the transition in MSM lists
#' @param n_sim number of simulations being run
#'
#' @return matrix of parameter values



model.rmvnorm.f <- function(MSM_est, MSM_cov, index, n_sim = n_sim){
  #browser()
  transition <- names(MSM_cov)[index] # transition name
  rows <- which(MSM_est$transition==transition) # row index of MSM_est for this transition
  dist.v <- paste(MSM_est$model[rows[1]]) # chosen distribution for this transition model
  
  # d.data: vector of model estimates <- flexsurvreg(Surv(t, event) ~ x1 + x2 + x3, dist = "...")$res[,1]
  # if multiple models are chosen, need to select estimates for the first model - dist.v model
  rows <- rows[which(MSM_est$model[rows] == dist.v)]
  d.data <- MSM_est$estimate[rows]
  
  # vc.data: variance covariance matrix <- flexsurvreg(Surv(t, event) ~ x1 + x2 + x3, dist = "...")$cov
  vc.data <- as.matrix(MSM_cov[[index]][-1])
  vc.data <- round(vc.data, 6)
  
  # function returns transition probability
  ### GAMMA, WEIBULL, LOGLOGISTIC ###
  if(dist.v == "Gamma" | dist.v == "Weibull" | dist.v == "LogLogistic"){
    
    # 1) Transform MSM model est
    logPar1 <- log(d.data[1]) # log estimate paramter 1
    logPar2 <- log(d.data[2]) # log estimate paramter 2
    if(length(d.data)>2){ # covariate coefficients
      B <- d.data[3:length(d.data)]
      op.par <- c(logPar1, logPar2, B)
    }else {
      op.par <- c(logPar1, logPar2)
    }
    # 2) Run MSM output in rmvnorm - multivariate normal distribution with mean = op.par and sigma = vc.data (variance covariance matrix)
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F) # Multivariate normal distribution sample # of simulations
    
    #3) Transform rmvnorm MSM model est
    t.par1 <- exp(mvnorm.res[,1]) # transform multivariate normal estimates
    t.par2 <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }
    
    ### CURE GAMMA, WEIBULL, LOGLOGISTIC ###
  } else if(dist.v == "Gamma.c" | dist.v == "Weibull.c" | dist.v == "LogLogistic.c"){
    Par1 <- log(d.data[1]/(1 - d.data[1])) # theta estimate
    logPar2 <- log(d.data[2]) # log estimate paramter 2
    logPar3 <- log(d.data[3]) # log estimate paramter 3
    if(length(d.data)>3){
      B <- d.data[4:length(d.data)]
      op.par <- c(Par1, logPar2, logPar3, B)
    }else {
      op.par <- c(Par1, logPar2, logPar3)
    }
    
    
    
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F) # Multivariate normal distribution
    
    t.par1 <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2 <- exp(mvnorm.res[,2])
    t.par3 <- exp(mvnorm.res[,3])
    if(length(d.data)>3){
      Beta <- mvnorm.res[,4:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, t.par3, Beta)
    }else {
      results <- cbind(t.par1, t.par2, t.par3)
    }
    
    ### EXPONENTIAL ###
  } else if (dist.v == "Exponential"){
    
    logPar1 <- log(d.data[1]) # log rate estimate
    if(length(d.data)>1){
      B <- d.data[2:length(d.data)]
      op.par <- c(logPar1, B)
    }else {
      op.par <- c(logPar1)
    }
    
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)
    
    t.par1 <- exp(mvnorm.res[,1])
    if(length(d.data)>1){
      Beta <- mvnorm.res[,2:ncol(mvnorm.res)]
      results <- cbind(t.par1, Beta)
    }else {
      results <- cbind(t.par1)
    }
    
    ### CURE EXPONENTIAL ###
  } else if (dist.v == "Exponential.c"){
    Par1 <- log(d.data[1]/(1 - d.data[1])) # theta estimate
    logPar2 <- log(d.data[2]) # log rate estimate
    if(length(d.data)>2){
      B <- d.data[3:length(d.data)]
      op.par <- c(Par1, logPar2, B)
    }else {
      op.par <- c(Par1, logPar2)
    }
    
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)
    t.par1 <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2 <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }
    
    ### LOGNORMAL, GOMPERTZ ###
  } else if (dist.v == "LogNormal"| dist.v == "Gompertz"){
    Par1 <- d.data[1] # meanlog(LogNormal). shape(Gompertz)
    logPar2 <- log(d.data[2]) # transform sdlog(LogNormal) or rate(Gompertz) estimate
    if(length(d.data)>2){
      B <- d.data[3:length(d.data)]
      op.par <- c(Par1, logPar2, B)
    }else {
      op.par <- c(Par1, logPar2)
    }
    
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)
    
    t.par1 <- mvnorm.res[,1]
    t.par2 <- exp(mvnorm.res[,2])
    if(length(d.data)>2){
      Beta <- mvnorm.res[,3:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, Beta)
    }else {
      results <- cbind(t.par1, t.par2)
    }
    
    ### CURE LOGNORMAL, GOMPERTZ ###
  } else if (dist.v == "LogNormal.c"| dist.v == "Gompertz.c"){
    Par1 <- log(d.data[1]/(1 - d.data[1])) # theta
    Par2 <- d.data[2] # meanlog(LogNormal). shape(Gompertz)
    logPar3 <- log(d.data[3]) # transform sdlog(LogNormal) or rate(Gompertz) estimate
    if(length(d.data)>3){
      B <- d.data[4:length(d.data)]
      op.par <- c(Par1, Par2, logPar3, B)
    }else {
      op.par <- c(Par1, Par2, logPar3)
    }
    
    
    
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)
    
    t.par1 <- exp(mvnorm.res[,1])/(1 + exp(mvnorm.res[,1]))
    t.par2 <- mvnorm.res[,2]
    t.par3 <- exp(mvnorm.res[,3])
    if(length(d.data)>3){
      Beta <- mvnorm.res[,4:ncol(mvnorm.res)]
      results <- cbind(t.par1, t.par2, t.par3, Beta)
    }else {
      results <- cbind(t.par1, t.par2, t.par3)
    }
    
    ### SPLINE ###
  } else if (dist.v == "Spline"){
    op.par <- d.data # model estimates
    mvnorm.res <- mvtnorm::rmvnorm(n_sim, op.par, vc.data, checkSymmetry = F)
    results <- mvnorm.res
    
  } else {
    print("no distribution found")
    results <- NA
  }
  return(results)
}
