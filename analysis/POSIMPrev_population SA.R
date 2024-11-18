######### Scenario analysis with low / high population growth scenarios

ON.sim_demog <- readRDS("~/GitHub/POSIM-prev/data/ON.sim_demog.Rds") #load sample paths of the projected population by sex and single year of age (generated in pop_projection.Rmd file)

no <- 500
no_slice <- no/4

ON.sim_demog[["male"]] <- ON.sim_demog[["male"]][,,1:no] #Take first 500/1,000 sample paths
ON.sim_demog[["female"]] <- ON.sim_demog[["female"]][,,1:no] #Take first 500/1,000 sample paths


#for each sample path - sum pop across all ages/years into a final total
maxpop_m <- apply(ON.sim_demog[["male"]][,,],3,sum)
maxpop_m <- as.data.frame(maxpop_m)

maxpop_f <- apply(ON.sim_demog[["female"]][,,],3,sum)
maxpop_f <- as.data.frame(maxpop_f)

maxpop <- cbind(maxpop_m,maxpop_f)

maxpop <- maxpop %>% dplyr::mutate(maxpop_sum = maxpop_m + maxpop_f, sample_path = rep(1:no))

maxpop_low <- maxpop %>% dplyr::slice_min(maxpop_sum, n = no_slice) #take the 125 sample paths with lowest maxpop_sum value

maxpop_high <- maxpop %>% dplyr::slice_max(maxpop_sum, n = no_slice) #take the 125 sample paths with highest maxpop_sum value

maxpop_low_v <- maxpop_low$sample_path #save the string of lowest sample paths as a vector

maxpop_high_v <- maxpop_high$sample_path #save the string of highest sample paths as a vector

# setwd("~/GitHub/POSIM-prev/analysis")
# saveRDS(maxpop_low_v, "maxpop_low_v.rds")
# saveRDS(maxpop_high_v, "maxpop_high_v.rds")


######## child pop

#for each sample path - sum pop across ages 0-14/all years into a final total
rm(list=ls())

ON.sim_demog <- readRDS("~/GitHub/POSIM-prev/data/ON.sim_demog.Rds")

ON.sim_demog[["male"]] <- ON.sim_demog[["male"]][1:15,,1:no] #Take first 500/1,000 sample paths, children 0-14
ON.sim_demog[["female"]] <- ON.sim_demog[["female"]][1:15,,1:no] #Take first 500/1,000 sample paths, children 0-14

maxpop_m <- apply(ON.sim_demog[["male"]][,,],3,sum)
maxpop_m <- as.data.frame(maxpop_m)

maxpop_f <- apply(ON.sim_demog[["female"]][,,],3,sum)
maxpop_f <- as.data.frame(maxpop_f)

maxpop <- cbind(maxpop_m,maxpop_f)

maxpop <- maxpop %>% dplyr::mutate(maxpop_sum = maxpop_m + maxpop_f, sample_path = rep(1:no))

maxpop_low <- maxpop %>% dplyr::slice_min(maxpop_sum, n = no_slice) #take the 125 sample paths with lowest maxpop_sum value

maxpop_high <- maxpop %>% dplyr::slice_max(maxpop_sum, n = no_slice) #take the 125 sample paths with highest maxpop_sum value

maxpop_low_child_v <- maxpop_low$sample_path #save the string of lowest sample paths as a vector

maxpop_high_child_v <- maxpop_high$sample_path #save the string of highest sample paths as a vector

#Save
# setwd("~/GitHub/POSIM-prev/analysis")
# saveRDS(maxpop_low_child_v, "maxpop_low_child_v.rds")
# saveRDS(maxpop_high_child_v, "maxpop_high_child_v.rds")


## Generate results ----------------------------------------------------------------------------------------------------

# *** Note - the below code is for illustration purposes. To run the code below, user is first required to run the model in POSIMPrev_microsimulation.R

# rm(list=ls())
# setwd("~/GitHub/POSIM-prev/analysis")
# 
# maxpop_low_v <- readRDS("maxpop_low_v.rds")
# maxpop_high_v <- readRDS("maxpop_high_v.rds")
# maxpop_low_child_v <- readRDS("maxpop_low_child_v.rds")
# maxpop_high_child_v <- readRDS("maxpop_high_child_v.rds")
# 
# if (!require('pacman')) install.packages('pacman'); library(pacman)
# 
# p_load("tidyverse", "reshape2", "knitr", "markdown", "data.table", "zoo", "stats", "here")
# 
# results <- c("prev_sim_kk", "prev_rate_kk",  "dxyear_sim_imm_kk")
# 
# init_year <- 1970               # starting year for the microsimulation model
# max_year  <- 2040               # end year
# 
# 
# #Read in results from microsimulation iterations
# 
# event.path <- here::here("analysis", "microsim_output")
# batchno <- 1:56 #number of batches run
# 
# 
# batchresults <- vector(mode = "list", length = length(results)) #for each result - catch the results from each batch/iteration
# names(batchresults) <- c(results)
# 
# for(i in 1:length(batchresults)){ #i = # of different results
#   
#   batchresults_kk <- list() #for each result - store results from all iterations in a sub-list
#   # for(b in 1:length(batchno)){
#   for(b in batchno){
#     
#     batchresults_kk[[b]] <-  read.csv(paste0(event.path, "/batch", b,"rawdataoutput/", names(batchresults)[i], ".csv"))
#   }
#   batchresults[[i]] <- batchresults_kk
#   
# }
# 
# #### for each result, bind cols from all batches
# results_kk <- vector(mode = "list", length = length(results))
# names(results_kk) <- c(results)
# 
# for(i in 1:length(results_kk)){
#   
#   results_kk[[i]] <- do.call(cbind, batchresults[[i]])
#   
#   results_kk[[i]] <- results_kk[[i]][,1:500]
#   
# }
# 
# 
# 
# #1970-2040 prevalence, all cancers combined - LOW growth scenario
# res_prev_low <- apply(results_kk[["prev_sim_kk"]][, maxpop_low_v], 1, quantile, c(0.025,0.5,0.975), na.rm = TRUE)
# res_prev_low <- t(res_prev_low)
# res_prev_low <- data.table(res_prev_low) %>% dplyr::rename(prev_sim = "50%", lower = "2.5%", upper = "97.5%") %>% dplyr::mutate(year =  rep((init_year):(max_year)))
# res_prev_low <- data.table(sapply(res_prev_low, round))
# res_prev_low$scenario <- rep("low pop")
# res_prev_low <- res_prev_low %>% filter(year == 2025 | year == 2030 | year == 2035 | year == 2040)
# 
# #1970-2040 prevalence rate, all cancers combined - LOW growth scenario
# res_prev_rate_low <- apply(results_kk[["prev_rate_kk"]][, maxpop_low_v], 1, quantile, c(0.025,0.5,0.975), na.rm = TRUE)
# res_prev_rate_low <- t(res_prev_rate_low)
# res_prev_rate_low <- data.table(res_prev_rate_low) %>% dplyr::rename(rate = "50%", lower_rate = "2.5%", upper_rate = "97.5%") %>% mutate(year =  rep((init_year):(max_year)))
# res_prev_rate_low <- data.table(sapply(res_prev_rate_low, round)) #round
# 
# res_prev_low <- left_join(res_prev_low, res_prev_rate_low, by = c("year"))
# 
# 
# #1970-2040 prevalence, all cancers combined - HIGH growth scenario
# res_prev_high<- apply(results_kk[["prev_sim_kk"]][, maxpop_high_v], 1, quantile, c(0.025,0.5,0.975), na.rm = TRUE)
# res_prev_high<- t(res_prev_high)
# res_prev_high<- data.table(res_prev_high) %>% dplyr::rename(prev_sim = "50%", lower = "2.5%", upper = "97.5%") %>% dplyr::mutate(year =  rep((init_year):(max_year)))
# res_prev_high <- data.table(sapply(res_prev_high, round))
# res_prev_high$scenario <- rep("high pop")
# res_prev_high <- res_prev_high %>% filter(year == 2025 | year == 2030  | year == 2035 | year == 2040)
# 
# #1970-2040 prevalence rate, all cancers combined - HIGH growth scenario
# res_prev_rate_high <- apply(results_kk[["prev_rate_kk"]][, maxpop_high_v], 1, quantile, c(0.025,0.5,0.975), na.rm = TRUE)
# res_prev_rate_high <- t(res_prev_rate_high)
# res_prev_rate_high <- data.table(res_prev_rate_high) %>% dplyr::rename(rate = "50%", lower_rate = "2.5%", upper_rate = "97.5%") %>% mutate(year =  rep((init_year):(max_year)))
# res_prev_rate_high <- data.table(sapply(res_prev_rate_high, round)) #round
# 
# res_prev_high <- left_join(res_prev_high, res_prev_rate_high, by = c("year"))
# 
# SApop_prev <- rbind(res_prev_low, res_prev_high)
# 
# ############ #
# 
# #overall incidence counts by year - LOW growth scenario
# res_dxyear_sim_imm <- apply(results_kk[["dxyear_sim_imm_kk"]][, maxpop_low_child_v], 1, quantile, c(0.025,0.5,0.975), na.rm = TRUE)
# res_dxyear_sim_imm <- t(res_dxyear_sim_imm)
# dxyear_sim_imm_df <- data.table(res_dxyear_sim_imm) %>% dplyr::rename(n_adj = "50%", lower = "2.5%", upper = "97.5%") %>% mutate(year =  rep((init_year):(max_year-1)))
# dxyear_sim_imm_df_low <- data.table(sapply(dxyear_sim_imm_df, round))
# dxyear_sim_imm_df_low$scenario <- rep("low pop")
# dxyear_sim_imm_df_low <- dxyear_sim_imm_df_low %>% filter(year == 2025 | year == 2030 | year == 2039)
# 
# #overall incidence counts by year - HIGH growth scenario
# res_dxyear_sim_imm <- apply(results_kk[["dxyear_sim_imm_kk"]][, maxpop_high_child_v], 1, quantile, c(0.025,0.5,0.975), na.rm = TRUE)
# res_dxyear_sim_imm <- t(res_dxyear_sim_imm)
# dxyear_sim_imm_df <- data.table(res_dxyear_sim_imm) %>% dplyr::rename(n_adj = "50%", lower = "2.5%", upper = "97.5%") %>% mutate(year =  rep((init_year):(max_year-1)))
# dxyear_sim_imm_df_high <- data.table(sapply(dxyear_sim_imm_df, round))
# dxyear_sim_imm_df_high$scenario <- rep("high pop")
# dxyear_sim_imm_df_high <- dxyear_sim_imm_df_high %>% filter(year == 2025 | year == 2030 | year == 2039)
# 
# SApop_inc <- rbind(dxyear_sim_imm_df_low, dxyear_sim_imm_df_high)





  