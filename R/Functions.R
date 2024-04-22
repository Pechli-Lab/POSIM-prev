#----------------------------------------------------------------------------#
####   Function to check if transition probability array/matrix  is valid ####
#----------------------------------------------------------------------------#
#' Check if transition array is valid
#'
#' \code{check_transition_probability} checks if transition probabilities are in \[0, 1\].
#'
#' @param a_P A transition probability array or matrix.
#' @param err_stop Logical variable to stop model run if set up as TRUE. Default = FALSE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = FALSE
#'
#' @return
#' This function stops if transition probability array is not valid and shows 
#' what are the entries that are not valid
#' @import utils
#' @export
check_transition_probability <- function(a_P,
                                         err_stop = FALSE, 
                                         verbose = FALSE) {
  
  a_P <- as.array(a_P)
  m_indices_notvalid <- arrayInd(which(a_P < 0 | a_P > 1), 
                                 dim(a_P))
  
  if(dim(m_indices_notvalid)[1] != 0){
    v_rows_notval   <- rownames(a_P)[m_indices_notvalid[, 1]]
    v_cols_notval   <- colnames(a_P)[m_indices_notvalid[, 2]]
    v_cycles_notval <- dimnames(a_P)[[3]][m_indices_notvalid[, 3]]
    
    df_notvalid <- data.frame(`Transition probabilities not valid:` = 
                                matrix(paste0(paste(v_rows_notval, v_cols_notval, sep = "->"),
                                              "; at cycle ",
                                              v_cycles_notval), ncol = 1), 
                              check.names = FALSE)
    
    if(err_stop) {
      stop("Not valid transition probabilities\n",
           paste(capture.output(df_notvalid), collapse = "\n"))
    }
    
    if(verbose){
      warning("Not valid transition probabilities\n",
              paste(capture.output(df_notvalid), collapse = "\n"))
    } 
  }
}

#----------------------------------------------------------------------------#
####   Function to check if sum of transition probabilities equal to one  ####
#----------------------------------------------------------------------------#
#' Check if the sum of transition probabilities equal to one. 
#'
#' \code{check_sum_of_transition_array} checks if each of the rows of the 
#' transition matrices sum to one. 
#' 
#' @param a_P A transition probability array.
#' @param n_states Number of health states.
#' @param n_t Number of cycles.
#' @param err_stop Logical variable to stop model run if set up as TRUE. 
#' Default = TRUE.
#' @param verbose Logical variable to indicate print out of messages. 
#' Default = TRUE
#' @return 
#' The transition probability array and the cohort trace matrix.
#' @import dplyr
#' @export
check_sum_of_transition_array <- function(a_P,
                                          n_states,
                                          n_t,  
                                          err_stop = TRUE, 
                                          verbose  = TRUE) {
  
  a_P <- as.array(a_P)
  d <- length(dim(a_P))
  # For matrix
  if (d == 2) {
    valid <- sum(rowSums(a_P))
    if (abs(valid - n_states)> 1e-04 ) {
      if(err_stop) {
        browser()
        stop("This is not a valid transition Matrix")
      }
      
      if(verbose){
        warning("This is not a valid transition Matrix")
      } 
    }
  } else {
    # For array
    valid <- (apply(a_P, d, function(x) sum(rowSums(x))) == n_states)
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n_t)))) {
      if(err_stop) {
        stop("This is not a valid transition Matrix")
      }
      
      if(verbose){
        warning("This is not a valid transition Matrix")
      } 
    }
  }
}


#---------------------------------------------------------------------------#
#### R function to sample states for multiple individuals simultaneously ####
#---------------------------------------------------------------------------#

# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. https://www.ncbi.nlm.nih.gov/pubmed/29587047

################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property 
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating 
# institutions and may not be used, reproduced, modified, distributed or adapted 
# in any way without appropriate citation.

################################################################################
# Developed by Petros Pechlivanoglou

samplev <- function(m.Probs) {
  # Arguments:
   # m.Probs: matrix with probabilities (n.i * n.s)
  # Returns:
   # ran: n.i x m matrix filled with sampled health state(s) per individual
  d <- dim(m.Probs) # dimensions of the matrix filled with the multinomical probabilities for the health states
  n <- d[1] # first dimension - number of rows (number of individuals to sample for)
  k <- d[2] # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]] # extract the names of the health states considered for sampling
  if (!length(lev)) # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k # create a sequence from 1:k (number of health states considered)
  # create a matrix
  ran <- rep(lev[1], n) # create the matrix ran, filled with the first health state of the levels
  U <- t(m.Probs) # transposed m.Probs matrix n.i x n.s --> n.s x n.i
  for(i in 2:k) { # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual
  }
  if (any((U[k, ] - 1) > 1e-05)) # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  un <- rep(runif(n), rep(k, n)) # sample from a uniform distribution of length n*k
  ran <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  ran # return the new health state per individual n.i x m
} # close the function


#---------------------------------------------------------------------------#
####                    R functions for visualization                    ####
#---------------------------------------------------------------------------#


# plot health state trace
plot_m_TR <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_states,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_states,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}

