# Parts of this material were developed using code from the Decision Analysis in R for Technologies in Health (DARTH) workgroup. 

# Please cite the following publication when using or adapting this code:
# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 
# https://journals.sagepub.com/doi/abs/10.1177/0272989X18754513

# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide. Copyright, 
# trademarks, trade names and any and all associated intellectual property are 
# exclusively owned by THE HOSPITAL FOR Sick CHILDREN and the collaborating 
# institutions. These materials may be used, reproduced, modified, distributed 
# and adapted with proper attribution.

## Setup ----------------------------------------------------------------------------------------------------
options(scipen=999)

## Load packages --------------------------------------------------------------------------------------------

rm(list=ls())
if (!require('pacman')) install.packages('pacman'); library(pacman)

p_load("here", "devtools", "scales", "ellipse", "ggplot2", "lazyeval", "igraph", "truncnorm", "ggraph", "reshape2", "knitr", "markdown", "stringr", "demography", "readr", "data.table", "zoo", "mvtnorm", "flexsurv", "flexsurvcure", "HMDHFDplus", "tidyr", "stats", "survival", "diagram", "ggforce", "survminer", "plotly", "plyr", "dplyr", "ggpubr", "gridExtra", "reshape")

#install_github("DARTH-git/dampack", force = TRUE)
p_load_gh("DARTH-git/dampack")

#install_github("DARTH-git/darthtools", force = TRUE)
p_load_gh("DARTH-git/darthtools")

setwd("~/GitHub/POSIM-prev/analysis")

## Create folder structure to save results from each iteration -----------------------------------------------

##We ran batches of 9 iterations at a time due to computational demand

##First create a holding folder called 'SA_inc_output'
# dir.create(paste0("~/GitHub/POSIM-prev/analysis/","SA_surv_output"))

##250/9 - need 27.7 runs of 9 sims to reach 250 iterations. Create 28 folders
# batch <- c(1:28)
# for (b in batch){
#   folder<-dir.create(paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/","batch",b,"rawdataoutput"))
# }

## Load functions ---------------------------------------------------------------------------------------------

source("~/GitHub/POSIM-prev/R/Functions.R") #microsim functions
source("~/GitHub/POSIM-prev/R/trans prob function.R") #transition probability calculation for individuals with cancer history


## Load model inputs derived from parametric multi-state survival modeling ------------------------------------

#MSM results for patients diagnosed 1970-1989
load("m_t_70to89.Rdata") # MSM matrix (diagnosis --> death)
l_MSM_est_70to89 <- readRDS("l_MSM_est_70to89_nocure.rds") #model estimates
norm.mat.all_70to89 <- readRDS("norm.mat.all_70to89_nocure_500.rds") #uncertainty
MSM_knots_70to89    <- read.csv("df_model_knots_70to89_nocure.csv") #knot locations for models that required it
MSM_knots_70to89 <- MSM_knots_70to89 %>% select(transition, knots, ICCC_regroup)
MSM_knots_70to89$transition <- factor(MSM_knots_70to89$transition, levels = unique(MSM_knots_70to89$transition))
MSM_knots_70to89$ICCC_regroup[MSM_knots_70to89$ICCC_regroup =="Non-Hodgkin lymphomas and other lymphomas\n"] = "Non-Hodgkin lymphomas and other lymphomas"

#MSM results for patients diagnosed 1990-2019
load("m_t.Rdata") # MSM matrix (diagnosis --> CRE --> death)
l_MSM_est    <- readRDS("l_MSM_est_90to19_nocure.rds") #model estimates
norm.mat.all <- readRDS("norm.mat.all_90to19_nocure_500.rds") #uncertainty
MSM_knots    <- read.csv("df_model_knots_90to19_nocure.csv") #knot locations for models that required it
MSM_knots <- MSM_knots %>% select(transition, knots, ICCC_regroup)
MSM_knots$transition <- factor(MSM_knots$transition, levels = unique(MSM_knots$transition))
MSM_knots$ICCC_regroup[MSM_knots$ICCC_regroup =="Non-Hodgkin lymphomas and other lymphomas\n"] = "Non-Hodgkin lymphomas and other lymphomas"


## Model specifications -----------------------------------------------------------------------------------------

batchno <- c(1:28) #for this scenario analysis, we ran 250 iterations

for (b in batchno) { #start batch loop

print(paste("Starting batch =", b)) #print which batch has started running
  
set.seed(b)

n_sim <- 9
max_n_sim <- 9*56

# Model structure
v_n        <- c("tobeborn", "noC", "C", "CRE", "adult",  "dead", "candead")  # vector with state names

n_states   <- length(v_n)     # number of states
n_i_init   <- 7000000         # number of individuals

country   <- "CAN"    # country to extract HMD data for
con.name  <- "Canada"
StatCan_start_year <- 1971      # StatCan data starts at 1971, not 1970
StatCan_end_year <- 2020

init_year <- 1970               # starting year for the microsimulation model

start_year_pred <- 2020         # starting year for predictions (births, population, incident cases, etc.)
max_year  <- 2040               # end year
max_age   <- 100                # Max age for the cohort to be followed
max_dx_age <- 14                # Max age at cancer diagnosis

dxages <- 0:14
v_years   <- init_year:max_year # the years that the model will allow people to be born.

sex_var <- c(0, 1)  # 1 = male, 0 = female
sex_var_f <- c("Male", "Female")

n_t       <- max_year - init_year # number of total cycles for the model

ICCCgroups <- c("ALL", "AML and other leukemias", "Astrocytoma", "Bone tumours", "Germ cell tumours", "Hepatic tumours", "Hodgkin lymphomas", "Neuroblastoma", "Non-Hodgkin lymphomas and other lymphomas", "Other CNS neoplasms", "Other epithelial and unspecified neoplasms", "Renal tumours", "Retinoblastoma", "Soft tissue sarcomas") #in alphabetical order - used in ICES survival models

attainedage_group <- c("0-19", "20-39", "40-59", "60+")
fiveyrperiods <- c("2020-2024", "2025-2029", "2030-2034", "2035-2039")

attainedage_group2 <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
                        "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+")

decade_dx <- c("1970-79", "1980-89", "1990-99", "2000-09", "2010-19", "2020-29", "2030-39")
decade_dx_long <- c("decade=1970-79", "decade=1980-89", "decade=1990-99", "decade=2000-09", "decade=2010-19", "decade=2020-29", "decade=2030-39")
decade_time <- rep(c(0:70, 0:61, 0:51, 0:41, 0:31, 0:21, 0:11), times = 1)
decade_strata_long <- rep(decade_dx_long, times = c(length(0:70), length(0:61), length(0:51), length(0:41), length(0:31), length(0:21),length(0:11)))
template_surv <- data.frame(time = decade_time, strata = as.factor(decade_strata_long))

period_dx_long <- c("period=1970-74", "period=1975-79", "period=1980-84", "period=1985-89", "period=1990-94", "period=1995-99", "period=2000-04", "period=2005-09", "period=2010-14", "period=2015-19", "period=2020-24", "period=2025-29", "period=2030-34", "period=2035-39")
period_time <- rep(c(0:70, 0:66, 0:61, 0:56, 0:51, 0:46, 0:41, 0:36, 0:31, 0:26, 0:21, 0:16, 0:11, 0:6), times = 1)
period_strata_long <- rep(period_dx_long, times = c(length(0:70), length(0:66), length(0:61), length(0:56), length(0:51), length(0:46), length(0:41), length(0:36), length(0:31), length(0:26), length(0:21), length(0:16), length(0:11), length(0:6)))
template_surv5 <- data.frame(time = period_time, strata = as.factor(period_strata_long))

survstatus_group <- c("<5", "5-9", "10-14", "15-19", "20+")
survstatus_group2 <- c("0-1 years", "2 years", "3 years", "4 years", "5+ years")

## Generate projections of mortality rates from lee-carter modeling ----------------------------------------------

# source("~/GitHub/POSIM-prev/R/hmd_function.R") #human mortality database function
# mort_con  <- hmd.mx2(country, "petros.pechlivanoglou@sickkids.ca","DARTH", Ontario = T)
# saveRDS(mort_con, file = "mort_con_backup.Rds")

mort_con <- readRDS("~/GitHub/POSIM-prev/data/mort_con_backup.Rds") #mortality rates from the Human Mortality Database
small  <- extract.years(mort_con, years <- c(init_year:max_year))
small  <- extract.ages(small, 0:max_age, FALSE)

# Lee-Carter mortality projections
start_year_lca <- start_year_pred
end_year_lca   <- max_year

lca_m      <- lca(small, series = "male",   interpolate = T, max.age = max_age)
lca_fort_m <- forecast(lca_m, h = end_year_lca - start_year_lca + 1)
lca_f      <- lca(small, series = "female",   interpolate = T, max.age = max_age)
lca_fort_f <- forecast(lca_f, h = end_year_lca - start_year_lca + 1)

#mortality projections for all ages - 0-110+: for use in no cancer and cancer states
p_mort_female_all <- cbind(small$rate$female, lca_fort_f$rate$female)
p_mort_male_all   <- cbind(small$rate$male, lca_fort_m$rate$male)
colnames(p_mort_female_all) <- colnames(p_mort_female_all) <- init_year: (init_year + ncol(p_mort_male_all)-1)
colnames(p_mort_male_all) <- colnames(p_mort_male_all) <- init_year: (init_year + ncol(p_mort_male_all)-1)

dim_all <- length(p_mort_female_all)

df_mort <- data.frame(
  sexMale  = rep(c(1, 0), each = dim_all), #1 = male, 0 = female
  
  year_current = c(rep(colnames(p_mort_female_all), each = nrow(p_mort_female_all)),
                   rep(colnames(p_mort_male_all)  , each = nrow(p_mort_male_all))) ,
  
  age_current = c(rep(rownames(p_mort_female_all), times = ncol(p_mort_female_all)),
                  rep(rownames(p_mort_male_all)  , times = ncol(p_mort_male_all))),
  
  Rate   = c(matrix(p_mort_female_all, ncol = 1),
             matrix(p_mort_male_all,   ncol = 1))
)
df_mort$year_current <- as.numeric(df_mort$year_current)
df_mort$year_current <- rep(init_year:max_year, each = length(unique(df_mort$age_current)), times = length(unique(df_mort$sexMale)))
df_mort$age_current  <- as.numeric(df_mort$age_current)


## Load cancer incidence rate projections from generalized additive modeling -------------------------------------

predrates <- read.csv("~/GitHub/POSIM-prev/data/GAM_incidence_results.csv")

predrates <- predrates %>%
  select(!X) %>%
  dplyr::rename(cancer_type = group, year = P, age = A) %>%
  dplyr::filter(year >= init_year)

ocr_pogo_years <- init_year:max_year #predicted years for incidence rates in GAM models with OCR/POGO/ICES data

combinations <-length(ICCCgroups)*length(ocr_pogo_years)*length(unique(predrates$age))
predrates$sex <- rep(0:1, each = combinations)
predrates$sex <- factor(predrates$sex,
                        levels = c("0", "1"),
                        labels = c("female", "male"))

y <- c(predrates$cancer_type) #here cancer type variable is fixed to remove _Female and _Male from cancer type group names
y2 <- unlist(strsplit(y, split = "_"))
predrates$cancer_type <- y2[seq(1,length(y)*2, by=2)]

predrates$cancer_type <- as.factor(predrates$cancer_type) #prints in alphabetical order, which was used in the surv. models but ICCC groups were ordered by their roman numerals for ICES incidence analyses.
ICCCgroups_incidence <- c("ALL", "AML and other leukemias", "Hodgkin lymphomas", "Non-Hodgkin lymphomas and other lymphomas", "Astrocytoma", "Other CNS neoplasms", "Neuroblastoma", "Retinoblastoma", "Renal tumours", "Hepatic tumours", "Bone tumours", "Soft tissue sarcomas", "Germ cell tumours", "Other epithelial and unspecified neoplasms")
predrates$cancer_type <- factor(predrates$cancer_type, levels = unique(ICCCgroups_incidence))

## Vary estimates of late mortality relative risk for childhood cancer survivors ---------------------------------

#load eTable5 information from Yeh et al. JAMA Oncology 2020 - 'Life Expectancy of Adult Survivors of Childhood Cancer Over 3 Decades'; doi: 10.1001/jamaoncol.2019.5582.
Yeh <- read.csv("~/GitHub/POSIM-prev/data/Yeh_JAMA2020_eTable5.csv") 
colnames(Yeh)[1] <- "Cohort"
Yeh <- Yeh %>% dplyr::filter(Cohort == "Overall")
Yeh <- rbind(Yeh, Yeh [3,]) #apply 1990-99 RRs to 2000-09
Yeh$Decade[4] = "   2000-09"

Yeh$Decade_start <- as.numeric(substr(Yeh$Decade,start = 3, stop = 7))
Yeh$Decade_stop <- Yeh$Decade_start+9

#age 30 RRs/SE
Yeh$age_30_SE <- (log(Yeh$age_30) - log(Yeh$LI_30))/1.96
Yeh$age_30_SE_2 <- (log(Yeh$UI_30) - log(Yeh$age_30))/1.96
Yeh$age_30_SE_3 <- colMeans(rbind(Yeh$age_30_SE, Yeh$age_30_SE_2))

#age 40 RRs/SE
Yeh$age_40_SE <- (log(Yeh$age_40) - log(Yeh$LI_40))/1.96
Yeh$age_40_SE_2 <- (log(Yeh$UI_40) - log(Yeh$age_40))/1.96
Yeh$age_40_SE_3 <- colMeans(rbind(Yeh$age_40_SE, Yeh$age_40_SE_2))

#age 50 RRs/SE
Yeh$age_50_SE <- (log(Yeh$age_50) - log(Yeh$LI_50))/1.96
Yeh$age_50_SE_2 <- (log(Yeh$UI_50) - log(Yeh$age_50))/1.96
Yeh$age_50_SE_3 <- colMeans(rbind(Yeh$age_50_SE, Yeh$age_50_SE_2))

#age 60 RRs/SE
Yeh$age_60_SE <- (log(Yeh$age_60) - log(Yeh$LI_60))/1.96
Yeh$age_60_SE_2 <- (log(Yeh$UI_60) - log(Yeh$age_60))/1.96
Yeh$age_60_SE_3 <- colMeans(rbind(Yeh$age_60_SE, Yeh$age_60_SE_2))

Yeh_df <- data.frame(Cohort = Yeh$Cohort, Decade_start = Yeh$Decade_start, Decade_stop = Yeh$Decade_stop, age_30_log = log(Yeh$age_30), age_30_SE = Yeh$age_30_SE_3, age_40_log = log(Yeh$age_40), age_40_SE = Yeh$age_40_SE_3, age_50_log = log(Yeh$age_50), age_50_SE = Yeh$age_50_SE_3, age_60_log = log(Yeh$age_60), age_60_SE = Yeh$age_60_SE_3)

Yeh_df_melt <- reshape2::melt(Yeh_df, id.vars=c("Decade_start", "Decade_stop", "Cohort"))

Yeh_mean <- Yeh_df_melt %>% filter(substring(variable, first = 8) == "log")
Yeh_mean <- Yeh_mean %>% slice(rep(1:n(), each = 10))
Yeh_mean <- Yeh_mean[order(Yeh_mean$Decade_start, Yeh_mean$variable),]
ages <- length(unique(Yeh_mean$variable))
decades <- length(unique(Yeh_mean$Decade_start))
Yeh_mean$Age <- rep(30:69, decades)

Yeh_se <- Yeh_df_melt %>% filter(substring(variable, first = 8) == "SE")
Yeh_se <- Yeh_se %>% slice(rep(1:n(), each = 10))
Yeh_se <- Yeh_se[order(Yeh_se$Decade_start, Yeh_se$variable),]
Yeh_se$Age <- rep(30:69, decades)

#gen vals from each dist.
k = 0
RRmort <- array(dim = c(length(min(Yeh_mean$Age):max(Yeh_mean$Age)), decades, max_n_sim), dimnames = list(min(Yeh_mean$Age):max(Yeh_mean$Age), unique(Yeh_mean$Decade_start), 1:max_n_sim))
for (j in 1:decades){
  for (i in 1:length(min(Yeh_mean$Age):max(Yeh_mean$Age))){
    k = k+1
    RRmort[i,j,] <- rlnorm(max_n_sim, Yeh_mean$value[[k]], Yeh_se$value[[k]])
    
  }
}
RRmort_final <- array(dim = c(length(min(Yeh_mean$Age):max(Yeh_mean$Age)), length(init_year:max(Yeh_mean$Decade_stop)), max_n_sim), dimnames = list(min(Yeh_mean$Age):max(Yeh_mean$Age), init_year:max(Yeh_mean$Decade_stop), 1:max_n_sim))

RR_long <- list()

for (k in 1:max_n_sim){
  RRmort_final[,,k] <- matrix(rep(t(RRmort[,,k]), each=10), ncol = length(init_year:max(Yeh_mean$Decade_stop)), byrow=T)
  RR_long[[k]]<- reshape2::melt(RRmort_final[,,k])
  colnames(RR_long[[k]]) <- c("age_current", "dxyear", "RR")
}

## Organize historical population data -------------------------------------------------------------------------------

# HMD_pop_all <- readCHMDweb(provID = "ont", item = "Population", fixup = TRUE)
# write.csv(HMD_pop_all, "../Population prediction/HMD_pop_backup.csv", row.names=FALSE)

HMD_pop_all <- read.csv("~/GitHub/POSIM-prev/data/HMD_pop_backup.csv") #Population estimates from the Human Mortality Database

HMD_pop <- HMD_pop_all %>%
  dplyr:::filter(Year >= init_year & Year < StatCan_start_year) %>%
  dplyr:::rename(pop_female = Female, pop_male = Male) %>%
  dplyr:::filter(Age < max_age) #pull up to 99. 100 - 110+ is summed separately below to get '100 and over'

#1970 only, 0-99
HMD_pop_m <- HMD_pop %>% 
  dplyr:::select(Year, Age, pop_male) %>%
  mutate(pop_male = round(pop_male), Sex = rep("male")) %>%
  dplyr:::rename(Pop = pop_male) %>%
  dplyr:::select(Year, Sex, Age, Pop)

#1970 only, 0-99
HMD_pop_f <- HMD_pop %>% 
  dplyr:::select(Year, Age, pop_female) %>%
  mutate(pop_female = round(pop_female), Sex = rep("female")) %>%
  dplyr:::rename(Pop = pop_female) %>%
  dplyr:::select(Year, Sex, Age, Pop)

HMDpop <- rbind(HMD_pop_m, HMD_pop_f)

HMD_pop_100over <- HMD_pop_all %>%
  dplyr:::filter(Year >= init_year & Year <= 2001) %>% #pull HMD pop estimates for 100 - 110+ for 1970 - 2001
  dplyr:::rename(pop_female = Female, pop_male = Male) %>%
  dplyr:::filter(Age >= max_age)

#sum to get '100 and over' for 1970-2001
HMD_pop_100over_m <- HMD_pop_100over %>% group_by(Year) %>% dplyr::summarise(pop_male = sum(pop_male)) %>% dplyr::rename(Pop = pop_male) %>% mutate(Sex = rep("male"), Age = rep(max_age), Pop = round(Pop))

HMD_pop_100over_f <- HMD_pop_100over %>% group_by(Year) %>% dplyr::summarise(pop_female = sum(pop_female)) %>% dplyr::rename(Pop = pop_female) %>% mutate(Sex = rep("female"), Age = rep(max_age), Pop = round(Pop))

HMD_pop_m_90to100 <- HMD_pop_all %>%
  dplyr:::filter(Year >= (init_year+1) & Year <= 2001) %>% #grab pop estimates for people 90-99 between 1971-2001 (what StatCan doesn't report)
  dplyr:::rename(pop_female = Female, pop_male = Male) %>%
  dplyr:::filter(Age >= (max_age-10) & Age < max_age) %>%
  dplyr:::select(Year, Age, pop_male) %>%
  mutate(pop_male = round(pop_male), Sex = rep("male")) %>%
  dplyr:::rename(Pop = pop_male) %>%
  dplyr:::select(Year, Sex, Age, Pop)

HMD_pop_f_90to100 <- HMD_pop_all %>%
  dplyr:::filter(Year >= (init_year+1) & Year <= 2001) %>% #grab pop estimates for people 90-99 between 1971-1990 (what StatCan doesn't report)
  dplyr:::rename(pop_female = Female, pop_male = Male) %>%
  dplyr:::filter(Age >= (max_age-10) & Age < max_age) %>%
  dplyr:::select(Year, Age, pop_female) %>%
  mutate(pop_female = round(pop_female), Sex = rep("female")) %>%
  dplyr:::rename(Pop = pop_female) %>%
  dplyr:::select(Year, Sex, Age, Pop)

HMD_pop_90to100 <- rbind(HMD_pop_m_90to100, HMD_pop_f_90to100)
HMD_pop_100over <- rbind(HMD_pop_100over_m, HMD_pop_100over_f)

### StatCan
popON <- read.csv("~/GitHub/POSIM-prev/data/pop_ON_Oct2021.csv") #Pop estimates for each year on JULY 1st (*StatCan captures period from July 1 to June 30.)
colnames(popON)[1] <- "REF_DATE"

#change pop data from representing pop on July 1st to Jan 1st
popON$Pop[1] <- popON$VALUE[1]
for (i in 1:nrow(popON)) {
  if(popON$REF_DATE[i] == StatCan_start_year) {
    popON$Pop[i]  <-   popON$VALUE[i] 
  } else{
    popON$Pop[i] <-   (popON$VALUE[i-1] + popON$VALUE[i])/2
  }
}
popON$Pop <- round(popON$Pop)

popON <- popON %>% 
  select(!VALUE) %>%
  dplyr:::rename(Year = REF_DATE, Age = Age.group)
popON$Age <- as.numeric(sub(" y.*","",popON$Age))

popON <- popON %>% dplyr:::filter(Year >= StatCan_start_year & Year < StatCan_end_year) #now contains pop estimate on JAN 1st of each year
popON$Sex <- factor(popON$Sex,
                    levels = c("Females", "Males"),
                    labels = c("female", "male"))
popON$Sex <- as.character(popON$Sex)

popON <- rbind(HMDpop, popON)

popON <- dplyr::union(popON, HMD_pop_90to100)
popON <- dplyr::union(popON, HMD_pop_100over)
popON <- subset(popON, !is.na(Pop))

## Load results from stochastic population forecasting ------------------------------------------------------------------
ON.sim_demog <- readRDS("~/GitHub/POSIM-prev/data/ON.sim_demog.Rds") #load sample paths of the projected population by sex and single year of age (generated in pop_projection.Rmd file)

## Create matrices to store results from each iteration of the microsimulation model ------------------------------------

adj_fac_m <- matrix(NA, nrow = n_sim, ncol = 2) #store adjustment factor used for population size in each iteration

prev_sim <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #1970-2040 annual prevalence
prev_rate <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #1970-2040 annual prevalence rate per 100,000 pop
prev_imm <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #1970-2040 annual prevalence - identify proportion diagnosed out of country
prev_bornON <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #1970-2040 annual prevalence - identify those born and diagnosed in Ontario
prev_imm_dxON <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #1970-2040 annual prevalence - identify those born outside Ontario, diagnosed in Ontario
prev_sim_true_prev <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #1990-2019 annual prevalence to compare with ICES cohort (internal validation)
prev_sim_true_prev_70 <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #1970-2019 annual prevalence to compare with ICES cohort (internal validation)
prev_sim_true_prev_ICCC <- matrix(NA, nrow = length(init_year:max_year)*length(ICCCgroups), ncol = n_sim) #1990-2019 ICES prevalence comparison by cancer type
prev_sim_true_prev_70_ICCC <- matrix(NA, nrow = length(init_year:max_year)*length(ICCCgroups), ncol = n_sim) #1970-2019 ICES prevalence comparison by cancer type
prev_sim_char <- matrix(NA, nrow = length(init_year:max_year)*length(sex_var)*length(decade_dx), ncol = n_sim) #1970-2040 annual prevalence by sex and decade of diagnosis (all cancers combined)
prev_sim_ICCC <- matrix(NA, nrow = length(init_year:max_year)*length(ICCCgroups), ncol = n_sim) #1970-2040 annual prevalence by cancer type
prev_sim_ICCC_char <- matrix(NA, nrow = length(init_year:max_year)*length(ICCCgroups)*length(sex_var)*length(decade_dx), ncol = n_sim) #1970-2040 annual prevalence by cancer type, sex, decade of diagnosis
prev_rate_ICCC <- matrix(NA, nrow = length(init_year:max_year)*length(ICCCgroups), ncol = n_sim) #1970-2040 annual prevalence rate per 100,000 pop by cancer type
dxyear_sim_imm <- matrix(NA, nrow = length((init_year+1):max_year), ncol = n_sim) #1970-2039 annual overall incidence counts - ON-born individuals + immigrants diagnosed in ON
dxyear_sim_imm_rate <- matrix(NA, nrow = length((init_year+1):max_year), ncol = n_sim) #1970-2039 annual overall incidence rates - ON-born individuals + immigrants diagnosed in ON
dxyear_sim_imm5_df <- matrix(NA, nrow = length(fiveyrperiods), ncol = n_sim) #incidence counts overall by 5-year time periods
dxyear_sim_imm5_rate_df <- matrix(NA, nrow = length(fiveyrperiods), ncol = n_sim) #incidence rates overall by 5-year time periods
dxyeargroup_sim <- matrix(NA, nrow = length((init_year+1):max_year)*length(ICCCgroups), ncol = n_sim) #1970-2039 annual incidence counts by cancer type
dxyeargroup_sim_age <- matrix(NA, nrow = length((init_year+1):max_year)*length(ICCCgroups)*length(sex_var)*length(dxages), ncol = n_sim) #1970-2039 annual incidence counts by cancer type, sex, age at diagnosis
dxyeargroup_iccc_df <- matrix(NA, nrow = length(ICCCgroups)*length(fiveyrperiods), ncol = n_sim) #incidence counts by 5-year time periods and cancer type
dxyeargroup_iccc_rate_df <- matrix(NA, nrow = length(ICCCgroups)*length(fiveyrperiods), ncol = n_sim) #incidence rates by 5-year time periods and cancer type
attainedage_year_df_m <- matrix(NA, nrow = length(attainedage_group)*(length(c(start_year_pred,max_year))), ncol = n_sim) # prevalence counts by attained age groups for 2020 and 2040
attainedage_aspr_df_m <- matrix(NA, nrow = length(attainedage_group2)*(length(c(start_year_pred,max_year))), ncol = n_sim) # prevalence by attained age groups for 2020 and 2040
attainedage_ICCC_m <- matrix(NA, nrow = length(attainedage_group)*length(ICCCgroups)*(length(c(start_year_pred,max_year))), ncol = n_sim) #prevalence by attained age group and cancer type for 2020 and 2040
attainedage_ICCC_aspr_m  <- matrix(NA, nrow = length(attainedage_group2)*length(ICCCgroups)*(length(c(start_year_pred,max_year))), ncol = n_sim) #prevalence by ASPR attained age groups and cancer type for 2020 and 2040
survstatus_df <- matrix(NA, nrow = length(survstatus_group)*(length(c(start_year_pred,max_year))), ncol = n_sim) #prevalence by length of FU since diagnosis, in 5 yr time periods. years 2020, 2040
survstatus_df2 <- matrix(NA, nrow = length(survstatus_group2)*(length(c(start_year_pred,max_year))), ncol = n_sim) #prevalence by length of FU since diagnosis, in 1 yr time periods until 5 yrs from diagnosis. years 2020, 2040
survstatus_ICCC_m <- matrix(NA, nrow = length(survstatus_group)*length(ICCCgroups)*(length(c(start_year_pred,max_year))), ncol = n_sim) #prevalence by length of FU since diagnosis in 5yr periods and cancer type for 2020 and 2040
survstatus2_ICCC_m <- matrix(NA, nrow = length(survstatus_group2)*length(ICCCgroups)*(length(c(start_year_pred,max_year))), ncol = n_sim) #prevalence by length of FU since diagnosis in 1yr periods until 5yrs and cancer type for 2020 and 2040
surv_kk <- matrix(NA,nrow = length(decade_time), ncol = n_sim) #overall survival by decade of diagnosis, all cancers combined
surv_kk_ALL <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_AML <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_NHL <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_HL <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_AST <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_OCNS <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_NEU <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_RET <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_REN <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_HEP <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_BONE <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_SAR <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_GERM <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_OTHER <- matrix(NA,nrow = length(decade_time), ncol = n_sim)
surv_kk_ICCC_5 <- matrix(NA, nrow = length(period_time)*length(ICCCgroups), ncol = n_sim) #survival by 5yr time periods, all ICCC groups in same list
surv_kk_5 <- matrix(NA, nrow = length(period_time), ncol = n_sim) #overall survival by 5yr time periods of diagnosis, all cancers combined
surv_pfs_kk <- matrix(NA,nrow = length(decade_time), ncol = n_sim) #event-free survival by decade of diagnosis (all cancers combined)
surv_pfs_kk_ICCC <- matrix(NA, nrow = length(decade_time)*length(ICCCgroups), ncol = n_sim) #event-free survival by decade of diagnosis, all ICCC groups in same list
childpop_PSA <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #estimates of population size used in each iteration, children aged 0-14
childpop_agesex_PSA <- matrix(NA, nrow = length(init_year:max_year)*length(sex_var)*length(dxages), ncol = n_sim)
prevpop_PSA <- matrix(NA, nrow = length(init_year:max_year), ncol = n_sim) #estimates of population size used in each iteration, all ages

## Set up parallel processing ------------------------------------------------------------------------------------

library(foreach)
library(doParallel)
parallel::detectCores()
# n.cores <- parallel::detectCores() - 1
n.cores <- 3
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK",
  outfile="log_SA_surv.txt"
)

clusterCall(my.cluster, function() lapply(c("here", "devtools", "scales", "ellipse", "ggplot2", "lazyeval", "igraph", "truncnorm", "ggraph", "reshape2", "knitr", "markdown", "stringr", "demography", "readr", "data.table", "zoo", "darthtools", "mvtnorm", "flexsurv", "flexsurvcure", "HMDHFDplus", "tidyr", "stats", "survival", "diagram", "ggforce", "survminer", "plotly", "plyr", "dplyr", "ggpubr", "gridExtra", "reshape"), library, character.only = TRUE))

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

x <- foreach(
  i = 1:10, 
  .combine = 'c'
) %do% {
  sqrt(i)
}


## Run simulation loop--------------------------------------------------------------------------------------------

sim_results <- foreach (k = 1:n_sim) %dopar% {
  rm(list = "outcomes")
  
  n_i <- n_i_init # number of individuals
  
  iteration <- ((b - 1) * n_sim + k)
  
  #assemble annual pop estimates by sex
  #each iteration uses a different sample path of the projected population from the stochastic population forecasting
  popm.mean <- ON.sim_demog$male[,,iteration]
  popf.mean <- ON.sim_demog$female[,,iteration]
  colnames(popf.mean) <- colnames(popm.mean) 
  
  popm.mean <- as.data.frame(popm.mean)
  popf.mean <- as.data.frame(popf.mean)
  
  #Males
  popON_m_pred <- reshape::melt(popm.mean, na.rm = FALSE)
  popON_m_pred <- popON_m_pred %>% 
    dplyr:::rename(Year = variable, Pop = value) %>%
    mutate(Sex = rep("male"),
           Age = rep(0:max_age, times = length(levels(Year))),
           Pop = round(Pop),
           Year = as.numeric(as.character(Year))) %>%
    dplyr:::select(Year, Sex, Age, Pop) %>%
    arrange(Age) #arrange by age so rolling mean can be calculated to get pop estimates for Jan 1st instead of July 1st
  
  popON_m_pred$Pop_Jan1[1] <- popON_m_pred$Pop[1]
  for (i in 1:nrow(popON_m_pred)) {
    if(popON_m_pred$Year[i] == start_year_pred) {
      popON_m_pred$Pop_Jan1[i]  <-   popON_m_pred$Pop[i] 
    } else{
      popON_m_pred$Pop_Jan1[i] <-   (popON_m_pred$Pop[i-1] + popON_m_pred$Pop[i])/2
    }
  }
  popON_m_pred$Pop_Jan1 <- round(popON_m_pred$Pop_Jan1)
  popON_m_pred <- popON_m_pred %>% select(!Pop)
  colnames(popON_m_pred) <- c("Year", "Sex", "Age", "Pop")
  
  #Females
  popON_f_pred <- reshape::melt(popf.mean, na.rm = FALSE)
  popON_f_pred <- popON_f_pred %>% 
    dplyr:::rename(Year = variable, Pop = value) %>%
    mutate(Sex = rep("female"),
           Age = rep(0:max_age, times = length(levels(Year))),
           Pop = round(Pop),
           Year = as.numeric(as.character(Year))) %>%
    dplyr:::select(Year, Sex, Age, Pop) %>%
    arrange(Age) 
  
  popON_f_pred$Pop_Jan1[1] <- popON_f_pred$Pop[1]
  for (i in 1:nrow(popON_f_pred)) {
    if(popON_f_pred$Year[i] == start_year_pred) {
      popON_f_pred$Pop_Jan1[i]  <-   popON_f_pred$Pop[i] 
    } else{
      popON_f_pred$Pop_Jan1[i] <-   (popON_f_pred$Pop[i-1] + popON_f_pred$Pop[i])/2
    }
  }
  popON_f_pred$Pop_Jan1 <- round(popON_f_pred$Pop_Jan1)
  popON_f_pred <- popON_f_pred %>% select(!Pop)
  colnames(popON_f_pred) <- c("Year", "Sex", "Age", "Pop")
  
  #bind predicted pop for males and females
  popON_pred <- rbind(popON_m_pred,popON_f_pred)
  popON_pred <- popON_pred %>% filter(Year >= start_year_pred & Year <= max_year)
  
  #### merge historical and predicted population data by age and sex
  pop <- rbind(popON,popON_pred)
  pop <- pop[order(pop$Year,pop$Age),]
  colnames(pop) <- c("year", "sex", "age", "pop")
  # pop$pop[is.na(pop$pop)] <- 0
  
  #Calculate net migration in children to use inside the microsim's nocancer state
  
  #'pop' df contains pop for all historical and future years on JAN 1 of each year
  pop.ON_f <- pop %>% dplyr:::filter(sex == "female" & age <= max_age) %>% dplyr:::select(!sex)
  pop.ON_f <- as.matrix(reshape2::dcast(data = pop.ON_f, age ~ year)[,-1])
  row.names(pop.ON_f) <- c(0:max_age) #assign matrix rows to correspond with the age range
  
  extra <- matrix(0, nrow = nrow(pop.ON_f), ncol = (nrow(pop.ON_f)-ncol(pop.ON_f)))
  row.names(extra) <- c(0:max_age)
  colnames(extra) <- c(1:ncol(extra))
  
  #create matrix of same # of rows and cols
  pop.ON_f <- cbind(pop.ON_f, extra)
  pop.ON_f <- as.matrix(rbind(rep(0,nrow(pop.ON_f)), pop.ON_f))
  pop.ON_f <- as.matrix(cbind(rep(0,nrow(pop.ON_f)+1), pop.ON_f))
  pop.ON_f <- as.matrix(pop.ON_f)
  
  #create lists for the netpop, year and ages
  list2 <- list()
  year2 <- list()
  age2 <- list()
  
  for(i in 1:ncol(pop.ON_f)){
    help.dat <- pop.ON_f[, i:ncol(pop.ON_f)]
    
    list2[[i]] <- diff(diag(help.dat))
    year2[[i]] <- colnames(help.dat)[-1]
    age2[[i]] <- rownames(help.dat)[-1][1:length(list2[[i]])]
    
  }
  
  results2 <- cbind(unlist(year2),  unlist(age2), unlist(list2[-length(list2)]))
  results2 <- data.table(results2)
  results2 <- results2 %>%   dplyr:::rename(year = V1, age = V2, netpop = V3)
  results2 <- sapply(results2, as.numeric)
  results2 <- data.table(results2)
  results2 <- filter(results2, year >= init_year) #to remove the diag diff results from the extra empty cols
  # results2 <- filter(results2, age <= max_dx_age) #retain only ages 0-14
  results2 <- results2[order(results2$year,results2$age),]
  
  #Get the remaining diag results in lists - for all other ages that weren't calculated above
  list3 <- list()
  year3 <- list()
  age3 <- list()
  
  for(i in 1:nrow(pop.ON_f)){
    help.dat <- pop.ON_f[i:nrow(pop.ON_f), ] #this time by nrow, not ncol
    
    list3[[i]] <- diff(diag(help.dat))
    year3[[i]] <- colnames(help.dat)[-1][1:length(list3[[i]])]
    age3[[i]] <- rownames(help.dat)[-1]
    
  }
  
  results3 <- cbind(unlist(year3),  unlist(age3), unlist(list3[-length(list3)]))
  results3 <- data.table(results3)
  results3 <- results3 %>%   dplyr:::rename(year = V1, age = V2, netpop = V3) 
  results3 <- sapply(results3, as.numeric)
  results3 <- data.table(results3)
  results3 <- filter(results3, year >= init_year) #to remove the diag diff results from the extra empty cols
  # results3 <- filter(results3, age <= max_dx_age)
  results3 <- results3[order(results3$year,results3$age),]
  
  #this df now contains the pop to simulate at each year within the microsim (# who enters (+) or exits (-) by year, age and sex)
  popON_net_f <- dplyr::union(results2, results3)
  popON_net_f <- popON_net_f[order(popON_net_f$year,popON_net_f$age),]
  popON_net_f$sex <- rep("female")
  
  #Males
  
  #Read in pop size for historical and future years. 'pop' df contains pop for all historical and future years on JAN 1 of each year
  pop.ON_m <- pop %>% dplyr:::filter(sex == "male" & age <= max_age) %>% dplyr:::select(!sex)
  pop.ON_m <- as.matrix(reshape2::dcast(data = pop.ON_m, age ~ year)[,-1])
  row.names(pop.ON_m) <- c(0:max_age)
  
  extra <- matrix(0, nrow = nrow(pop.ON_m), ncol = (nrow(pop.ON_m)-ncol(pop.ON_m)))
  row.names(extra) <- c(0:max_age)
  colnames(extra) <- c(1:ncol(extra)) #give new cols a smaller numeric value than the calendar years being assessed
  
  # #create matrix of same # of rows and cols
  pop.ON_m <- cbind(pop.ON_m, extra)
  pop.ON_m <- as.matrix(rbind(rep(0,nrow(pop.ON_m)), pop.ON_m))
  pop.ON_m <- as.matrix(cbind(rep(0,nrow(pop.ON_m)+1), pop.ON_m))
  pop.ON_m <- as.matrix(pop.ON_m)
  
  #create lists for the netpop, year and ages
  list2 <- list()
  year2 <- list()
  age2 <- list()
  
  for(i in 1:ncol(pop.ON_m)){
    help.dat <- pop.ON_m[, i:ncol(pop.ON_m)]
    
    list2[[i]] <- diff(diag(help.dat))
    year2[[i]] <- colnames(help.dat)[-1]
    age2[[i]] <- rownames(help.dat)[-1][1:length(list2[[i]])]
    
  }
  
  results2 <- cbind(unlist(year2),  unlist(age2), unlist(list2[-length(list2)]))
  results2 <- data.table(results2)
  results2 <- results2 %>%   dplyr:::rename(year = V1, age = V2, netpop = V3)
  results2 <- sapply(results2, as.numeric)
  results2 <- data.table(results2)
  results2 <- filter(results2, year >= init_year) #to remove the diag diff results from the extra empty cols
  # results2 <- filter(results2, age <= max_dx_age) #retain only ages 0-14
  results2 <- results2[order(results2$year,results2$age),]
  
  #Get the remaining diag results in lists - for all other ages that weren't calculated above
  list3 <- list()
  year3 <- list()
  age3 <- list()
  
  for(i in 1:nrow(pop.ON_m)){
    help.dat <- pop.ON_m[i:nrow(pop.ON_m), ] #this time by nrow, not ncol
    
    list3[[i]] <- diff(diag(help.dat))
    year3[[i]] <- colnames(help.dat)[-1][1:length(list3[[i]])]
    age3[[i]] <- rownames(help.dat)[-1]
    
  }
  
  results3 <- cbind(unlist(year3),  unlist(age3), unlist(list3[-length(list3)]))
  results3 <- data.table(results3)
  results3 <- results3 %>%   dplyr:::rename(year = V1, age = V2, netpop = V3) 
  results3 <- sapply(results3, as.numeric)
  results3 <- data.table(results3)
  results3 <- filter(results3, year >= init_year) #to remove the diag diff results from the extra empty cols
  # results3 <- filter(results3, age <= max_dx_age)
  results3 <- results3[order(results3$year,results3$age),]
  
  #this df now contains the pop to simulate at each year within the microsim (# who enters (+) or exits (-) by year, age and sex)
  popON_net_m <- dplyr::union(results2, results3)
  popON_net_m <- popON_net_m[order(popON_net_m$year,popON_net_m$age),]
  popON_net_m$sex <- rep("male")
  
  #Combine net migration results to use in microsim - popON_net contains starting pop for 1970 and net migrants for years 1971+
  
  popON_net <- rbind(popON_net_m, popON_net_f)
  popON_net[is.na(popON_net)] <- 0 #change any rows with pop size = NA to 0
  popON_net <- popON_net[order(popON_net$year,popON_net$age),]
  
  popON_net <- popON_net %>%
    dplyr::rename(pop = netpop) %>%
    dplyr::mutate(year_birth = year - age)
  popON_net$sim_pop <- ceiling(popON_net$pop / sum(popON_net$pop)*n_i)
  
  #starting pop - 1970
  popON_init_year <- popON_net %>% 
    dplyr::filter(year == init_year)
  
  #contains 0-100 pop for 1970, age 0 (new births cohorts) for 1971-2040, and immigrants 1-100 for 1971-2040
  v_n_popsim_female <- popON_net %>%
    dplyr::filter(sex == "female" & year > init_year & pop > 0 & age < max_age)
  v_n_popsim_female <- rbind((popON_init_year%>%filter(sex == "female")), v_n_popsim_female)
  v_n_popsim_female <- v_n_popsim_female %>% mutate(age_inityear = (year - year_birth - (year - init_year)))
  v_n_popsim_female$year_entry <- ifelse(v_n_popsim_female$age_inityear < 0, v_n_popsim_female$year_birth, init_year)
  v_n_popsim_female$age_entry <- ifelse(v_n_popsim_female$age_inityear < 0, (v_n_popsim_female$year_entry - v_n_popsim_female$year_birth), v_n_popsim_female$age_inityear)
  v_n_popsim_female$year_imm <- ifelse(v_n_popsim_female$year != v_n_popsim_female$year_birth, v_n_popsim_female$year, NA)
  v_n_popsim_female$year_imm[v_n_popsim_female$year == init_year] <- NA
  
  v_n_popsim_male <- popON_net %>%
    dplyr::filter(sex == "male" & year > init_year & pop > 0 & age < max_age)
  v_n_popsim_male <- rbind((popON_init_year%>%filter(sex == "male")), v_n_popsim_male)
  v_n_popsim_male <- v_n_popsim_male %>% mutate(age_inityear = (year - year_birth - (year - init_year)))
  v_n_popsim_male$year_entry <- ifelse(v_n_popsim_male$age_inityear < 0, v_n_popsim_male$year_birth, init_year)
  v_n_popsim_male$age_entry <- ifelse(v_n_popsim_male$age_inityear < 0, (v_n_popsim_male$year_entry - v_n_popsim_male$year_birth), v_n_popsim_male$age_inityear)
  v_n_popsim_male$year_imm <- ifelse(v_n_popsim_male$year != v_n_popsim_male$year_birth, v_n_popsim_male$year, NA)
  v_n_popsim_male$year_imm[v_n_popsim_male$year == init_year] <- NA
  
  ############# Emigration
  
  popON_net_f_remove <- popON_net %>%
    dplyr::filter(sex =="female") %>%
    # dplyr::filter(age <= max_dx_age) %>%
    dplyr::filter(pop < 0) %>%
    mutate(sim_pop = ceiling(pop / sum(popON_net$pop)*n_i)) %>%
    mutate(pop = pop*(-1), sim_pop = sim_pop*(-1)) %>%
    dplyr::filter(sim_pop > 0)
  
  popON_net_m_remove <- popON_net %>%
    dplyr::filter(sex =="male") %>%
    # dplyr::filter(age <= max_dx_age) %>%
    dplyr::filter(pop < 0) %>%
    mutate(sim_pop = ceiling(pop / sum(popON_net$pop)*n_i))  %>%
    mutate(pop = pop*(-1), sim_pop = sim_pop*(-1)) %>%
    dplyr::filter(sim_pop > 0)
  
  # Calculate probability of being diagnosed with each cancer type
  
  predrates <- predrates[order(predrates$year, predrates$age, predrates$sex),]
  
  set.seed(iteration*100+1000)
  
  predrates$rate_sim <- exp(rnorm(nrow(predrates),predrates$lograte,predrates$logSE))
  
  matrates <- matrix(predrates$rate_sim, nrow = length(ICCCgroups_incidence), dimnames = list(ICCCgroups_incidence)) # matrix of incidence rates
  matrates <- t(matrates)
  matrates <- proportions(matrates,1)
  #df of characteristics that matches the corresponding rates:
  rate_chars <- data.frame(year = rep(ocr_pogo_years, each = length(unique(predrates$age))*length(unique(predrates$sex))), 
                           age = rep(unique(predrates$age), each = length(unique(predrates$sex)), times = length(ocr_pogo_years)), 
                           sexMale = rep(0:1, times = (nrow(predrates)/(length(ICCCgroups)*length(unique(predrates$sex)))))) #0 = female, 1 = male
  
  matrates <- cbind(rate_chars, matrates)
  
  #rearrange to match the ICCCgroups' variable order
  matrates <- matrates %>% select(year, age, sexMale, ALL, 'AML and other leukemias', Astrocytoma, 'Bone tumours', 'Germ cell tumours', 'Hepatic tumours', 'Hodgkin lymphomas', Neuroblastoma, 'Non-Hodgkin lymphomas and other lymphomas', 'Other CNS neoplasms', 'Other epithelial and unspecified neoplasms', 'Renal tumours', Retinoblastoma, 'Soft tissue sarcomas')
  
  ######### #
  predrates_sum <- predrates %>% group_by(age, sex, year) %>% dplyr::summarise(rate_sum = sum(rate)) #overall incidence rate per 10^6 *person-years* is summed for all cancer types combined
  predrates_sum <- left_join(predrates_sum, pop, by = c("age", "sex", "year"))
  
  predrates_sum <- predrates_sum %>% dplyr::rename(sexMale = sex)
  predrates_sum$sexMale <- as.numeric(ifelse(predrates_sum$sexMale == "female", 0, 1))
  predrates_sum <- left_join(predrates_sum, matrates, by = c("year", "age", "sexMale")) #predrates_sum now contains matrates but with the corresponding sum of predicted rates and pop size for each group
  
  predrates_sum <- predrates_sum %>%
    mutate(p_NCC = ((rate_sum/10^6) * pop) / pop) %>% #overall probability of being diagnosed
    dplyr:::rename(year_current = year, age_current = age)
  
  incidenceprobs <- predrates_sum %>% select(year_current, age_current, sexMale, p_NCC)
  
  
  # 04 Sample individual level characteristics
  v_sex <- c(rep ("male"  , sum(v_n_popsim_male$sim_pop)),
             rep ("female", sum(v_n_popsim_female$sim_pop)))
  
  v_year_birth <- c(rep(v_n_popsim_male$year_birth, times = v_n_popsim_male$sim_pop), #retains birth year for all
                    rep(v_n_popsim_female$year_birth, times = v_n_popsim_female$sim_pop))
  
  #year they enter the model
  v_year_entry <- c(rep(v_n_popsim_male$year_entry, times = v_n_popsim_male$sim_pop),
                    rep(v_n_popsim_female$year_entry, times = v_n_popsim_female$sim_pop))
  
  #age they enter the model
  v_age_entry <- c(rep(v_n_popsim_male$age_entry, times = v_n_popsim_male$sim_pop),
                   rep(v_n_popsim_female$age_entry, times = v_n_popsim_female$sim_pop))
  
  #record the year they actually immigrate into ON
  v_year_imm <- c(rep(v_n_popsim_male$year_imm, times = v_n_popsim_male$sim_pop),
                  rep(v_n_popsim_female$year_imm, times = v_n_popsim_female$sim_pop))
  
  v_age_sim <- init_year - v_year_birth
  v_age_sim[v_age_sim <0 | v_year_entry > init_year] = NA #anyone with entry year after 1970 will have negative v_age_sim. assign NA as they won't enter the microsim as part of initial pop.
  
  n_i <- length(v_sex)
  
  #dataframe which outlines characteristics for each synthetic individual
  df_X  <- data.frame(ID = 1:n_i,
                      sexMale = (v_sex == "male")*1, #male = 1, female = 0
                      year_birth = v_year_birth, #year of birth
                      year_entry = v_year_entry, #year they enter the model
                      age_current = v_age_sim, #age at the current microsim cycle
                      year_current = NA,
                      year_imm = v_year_imm, #retain true year of immigration (Will become 1900 for Ontario residents)
                      year_emm = NA,
                      age = NA, #age at cancer diagnosis
                      dxyear = NA, #year of cancer diagnosis
                      ICCCgroup = NA, #cancer type
                      v_dx0 = 0, #time (years) since diagnosis
                      v_cre0 = 0, #time (years) since cancer-related event
                      v_death0 = 0, #time (years) since death
                      age_death = NA, #age of death
                      year_death = NA, #year of death
                      year_cre = NA, #year of CRE
                      age_cre = NA) #age of CRE
  
  # adj_fac_m[k,1] <- sum(v_n_popsim_male$pop)/sum(v_n_popsim_male$sim_pop) #save the specific adj_fac needed for each iteration
  # adj_fac_m[k,2] <- n_i #save the specific n_i used for each iteration
  adj_fac_m <- numeric(2) 
  adj_fac_m[1] <- sum(v_n_popsim_male$pop)/sum(v_n_popsim_male$sim_pop) #save the specific adj_fac needed for each iteration
  adj_fac_m[2] <- n_i #save the specific n_i used for each iteration
  
  # df_X <- df_X %>% mutate(`I(log(dxyear))` = NA) #unhide if using log(dxyear) model
  
  ## Initialize dynamic characteristics
  
  v_M_init  <- rep("tobeborn", times = n_i) #specify initial health state for all individuals
  
  v_M_init[df_X$year_entry == init_year & df_X$age_current < max_dx_age] = "noC"
  
  v_M_init[df_X$year_entry == init_year & df_X$age_current >= 15] = "adult"
  
  df_X$age_current[df_X$year_birth == init_year & df_X$year_entry == init_year] = 0
  
  df_X$year_current[df_X$year_entry == init_year] = init_year #assigning year_current to be initial year for those in the starting population
  
  Probs <- function(M_t, df_X0, t, kk) {
    #Main INPUT
    #Df_X_allocate : the whole cohort to be allocated at time t
    #incidenceprobs :incidence probabilities overall
    #df_mort:       background mortality
    
    #OUTPUT:
    #p_NCD <- probability of non cancer death for all df_X_allocate
    #p_NCC <- probability of cancer diagnosis for all
    #browser()
    
    #assign cancer type
    v_Cancer_new <- M_t =="C" & df_X0$v_dx0 == 1
    
    #if(sum(df_X0$v_dx0 == 1)>0)browser()
    predrates_sumtemp <- predrates_sum
    predrates_sumtemp$age_current <- predrates_sumtemp$age_current +1  
    rates_iccc <- inner_join(df_X0[v_Cancer_new,], predrates_sumtemp, by = c("sexMale", "year_current", "age_current"))%>% dplyr::select(all_of(ICCCgroups))
    
    df_X0[v_Cancer_new,"ICCCgroup"] <- samplev(rates_iccc)
    
    #look up baseline probability and rate of dying based on individual characteristics
    p_NCD_all  <- inner_join(df_X0, df_mort, by = c("sexMale", "year_current", "age_current"))
    p_NCD_RR   <- left_join(df_X0, RR_long[[iteration]], by = c("dxyear", "age_current"))[,"RR"]
    p_NCD_RR[is.na(p_NCD_RR)] <- 1
    p_NCD2      <- rate_to_prob(prob_to_rate(p_NCD_all[, "Rate"])*p_NCD_RR)
    
    p_NCC  <- left_join(df_X0[,c("sexMale", "year_current", "age_current")], incidenceprobs, by = c("sexMale", "year_current", "age_current"))
    p_NCC2  <- p_NCC[, "p_NCC"]
    p_NCC2[is.na(p_NCC2)] = 0
    p_CD <- p_CCRE <- p_CRED <- c()
  
    
    ######### calculate transition probabilities for NCD, CCRE, CRED, CD using msm results
    
    #df_X [tobeallocated,] <<- df_X0
    # loops over transitions, 
    for(trans in 1:sum(!is.na(m_t))){
      if(trans ==3){
        p.trans    <- rep(NA, sum(M_t == "CRE"))
        df_X_trans <- df_X0[M_t == "CRE",]
        df_X_trans$dxyear[df_X_trans$dxyear>= 2020] <- 2019 #Scenario analysis <- assign survival probabilities during future time horizon as if year of dx was 2019
      }else{
        p.trans   <- rep(NA, sum(M_t == "C"))
        df_X_trans <- df_X0[M_t == "C",]
        df_X_trans$dxyear[df_X_trans$dxyear>= 2020] <- 2019
      }
      #if(t==42) View(df_X_trans)
      transition <- which(m_t == trans)
      v_n_surv <- c("Diagnosis", "CRE", "Death")  
      s.from <- transition %% length(v_n_surv)    
      s.to   <- ceiling(transition/length(v_n_surv))
      
      #if(t==43)browser()
      
      # loops over ICCCgoups
      for(j in seq_along(ICCCgroups)){
        # pulling individuals within the Jth group 
        X_j <- df_X_trans[df_X_trans$ICCCgroup == ICCCgroups[j],]
        #if(t==42)View(X_j)
        
        # If someone is in the group then
        if(nrow(X_j) > 0){
          n_i_j <- nrow(X_j)
          
          if(trans == 3){
            X_j$v_t <-  X_j$v_cre0 
          }else{
            X_j$v_t <-X_j$v_dx0 #CRE --> Death is transition #3
          }
          
          
          #pulling Jth estimates from MSM
          l_MSM_est_j <- l_MSM_est[[j]]
          norm.mat_j <- norm.mat.all[[j]]
          cov.trans <- as.character(l_MSM_est_j[[trans]]$covariate)
          
          
          X_temp <-as.matrix(X_j[, cov.trans[cov.trans %chin% names(X_j)]])
          X_temp[,"dxyear"] =  X_temp[,"dxyear"] - 1989 #unhide if using original MSM results (dxyear,age,sex)
          X_temp[,"age"] =  X_temp[,"age"] - 1 #needed to capture 0-14 instead of 1-15
          
          #print(paste(trans,j, sep ="-"))
          #  if(t==45 & trans ==2 & j == 9) browser()
          
          #if(sum(X_j$v_t == 2 & df_X_trans$ICCCgroup == "Other CNS neoplasms") > 0)browser()
          p.trans[df_X_trans$ICCCgroup == ICCCgroups[j]] <- try(t(model.dist.f(dist.v = l_MSM_est_j[[trans]]$model[[1]], # model distribution
                                                                               #d.data = norm.mat_j[[trans]][kk,],  # multivariate normal estimates #original - relying on kk means only rows 1-9 were utilized in each batch. norm.mat_j had 100 rows, not 500 
                                                                               d.data = norm.mat_j[[trans]][iteration,],  # multivariate normal estimates
                                                                               dat.x  = X_temp, # matrix of baseline covariates
                                                                               t      = X_j$v_t, # time spent in states #original
                                                                               model = paste0("model_", v_n_surv[s.from], "_", v_n_surv[s.to]),
                                                                               step  = 1,
                                                                               n_i_j = n_i_j, 
                                                                               j = j,
                                                                               MSM_knots = MSM_knots)),silent = T)
          
          
          ######## specifications for patients diagnosed before 1990 ####### #
          
          if(trans == 2){ #diagnosis to death
            
            s.from <- 1    #manually assign transition from diagnosis
            s.to   <- 3    #manually assign transition to death
            
            l_MSM_est_70to89_j <- l_MSM_est_70to89[[j]] #MSM results for patients diagnosed 1970-1989
            norm.mat_70to89_j <- norm.mat.all_70to89[[j]] #MSM results for patients diagnosed 1970-1989
            cov.trans <- as.character(l_MSM_est_70to89_j[[1]]$covariate)
            
            X_j_70to89 <- X_j[X_j$dxyear < 1990,]
            
            X_temp_70to89 <-as.matrix(X_j_70to89[, cov.trans[cov.trans %chin% names(X_j_70to89)]])
            X_temp_70to89[,"dxyear"] =  X_temp_70to89[,"dxyear"] - 1969 #original MSM results analyzed patients diagnosed 1990-2019 (covariates: dxyear,age,sex)
            X_temp_70to89[,"age"] =  X_temp_70to89[,"age"] - 1 #needed to capture 0-14 instead of 1-15
            
            # X_temp[,"dxyear"] =  X_temp[,"dxyear"] + 20  #dxyear is already transformed in this matrix. Add back 20 so that 1970 = 1, etc. (1990 = 21, 1991 = 22, etc.)
            # n_i_j_70to89 <- nrow(X_j[X_j$dxyear < 1990,])
            
            n_i_j_70to89 <- nrow(X_j_70to89)
            
            p.trans[df_X_trans$ICCCgroup == ICCCgroups[j] & df_X_trans$dxyear < 1990] <- try(t(model.dist.f(dist.v = l_MSM_est_70to89_j[[1]]$model[[1]], # model distribution
                                                                                                            #d.data = norm.mat_70to89_j[[1]][kk,],  # multivariate normal estimates
                                                                                                            d.data = norm.mat_70to89_j[[1]][iteration,],  # multivariate normal estimates
                                                                                                            dat.x  = X_temp_70to89, # matrix of baseline covariates
                                                                                                            t      = X_j_70to89$v_t, # time spent in states
                                                                                                            model = paste0("model_", v_n_surv[s.from], "_", v_n_surv[s.to]),
                                                                                                            step  = 1,
                                                                                                            n_i_j = n_i_j_70to89, 
                                                                                                            j = j,
                                                                                                            MSM_knots = MSM_knots_70to89)),silent = T)
            
            #assign(paste0("p.", v_n_surv[s.from], "_", v_n_surv[s.to]), p.trans)
          }
          
          if(trans == 1){ #diagnosis to CRE
            
            s.from <- 1    #diagnosis
            s.to   <- 2    #CRE
            p.trans[df_X_trans$ICCCgroup == ICCCgroups[j] & df_X_trans$dxyear < 1990] <- 0
            
            #assign(paste0("p.", v_n_surv[s.from], "_", v_n_surv[s.to]), p.trans)
          }
          
          ####################################################################### # 
          
          assign(paste0("p.", v_n_surv[s.from], "_", v_n_surv[s.to]), p.trans)
        }
        
        else{  
          p_CCRE <- 0
          p_CD <- 0
          p_CRED <- 0}
        
      }
      
      try(assign("p_CCRE" , p.Diagnosis_CRE),silent = T) 
      try(assign("p_CD" , p.Diagnosis_Death),silent = T) 
      try(assign("p_CRED" , p.CRE_Death),silent = T) 
      
    }
    
    p_CD[is.na(p_CD)] = 0
    p_CCRE[is.na(p_CCRE)] = 0
    p_CRED[is.na(p_CRED)] = 0
    
    p_CD <- as.numeric(p_CD)
    p_CCRE <- as.numeric(p_CCRE)
    p_CRED <- as.numeric(p_CRED)
    
    # create matrix of state transition probabilities
    m_p_t           <- matrix(0, nrow = n_states, ncol = length(M_t))  
    # give the state names to the rows
    rownames(m_p_t) <-  v_n
    
    #age under 15 where still at risk
    age_at_risk <- df_X0$age_current <= max_dx_age #ages at risk are below 14
    p_NCD_NC <- p_NCD2[M_t == "noC" & age_at_risk] #background mortality risk for children in nocancer state
    p_NCC <- p_NCC2[M_t == "noC" & age_at_risk] #prob of cancer diagnosis if children in nocancer state
    p_NCD_C  <- p_NCD2[M_t == "C"] #background mortality risk for those in cancer state 
    p_NCD_CRE  <- p_NCD2[M_t == "CRE"] #background mortality risk for those in CRE state
    
    # update m_p_t with the appropriate probabilities (all non-death probabilities are conditional on survival)
    m_p_t["tobeborn",   M_t == "noC" & age_at_risk ] <- 0 # tobeborn
    m_p_t["noC",   M_t == "noC" & age_at_risk] <- 1 - p_NCD_NC - p_NCC #nocancer, grabs gen. pop risk of death for corresponding age
    m_p_t["C",     M_t == "noC" & age_at_risk ] <- p_NCC #cancer
    m_p_t["CRE",        M_t == "noC" & age_at_risk]  <- 0 #CRE
    m_p_t["adult",      M_t == "noC" & age_at_risk ] <- 0      #adult
    m_p_t["dead",       M_t == "noC" & age_at_risk ] <- p_NCD_NC   #dead
    m_p_t["candead", M_t == "noC" & age_at_risk ] <- 0   #cancerdead
    
    #when a child reaches age 15, they transition automatically to the adult stage (probability = 100%)
    m_p_t["tobeborn", M_t == "noC" & df_X0$age_current >= (max_dx_age + 1) ] <- 0 # tobeborn
    m_p_t["noC", M_t == "noC" & df_X0$age_current >=  (max_dx_age + 1) ] <- 0 #nocancer
    m_p_t["C",   M_t == "noC" & df_X0$age_current >=  (max_dx_age + 1 )] <- 0 # cancer
    m_p_t["CRE",      M_t == "noC" & df_X0$age_current >=  (max_dx_age + 1 )] <- 0 #CRE
    m_p_t["adult",    M_t == "noC" & df_X0$age_current >= ( max_dx_age + 1) ] <- 1  # adult                  
    m_p_t["dead",       M_t == "noC" & df_X0$age_current >=  (max_dx_age + 1) ] <- 0  # dead 
    m_p_t["candead", M_t == "noC" & df_X0$age_current >= ( max_dx_age + 1) ] <- 0  # cancerdead 
    
    #transition probabilities for those with a cancer diagnosis
    m_p_t["tobeborn", M_t == "C" ]  <- 0 # tobeborn
    m_p_t["noC", M_t == "C" ]  <- 0 # nocancer
    m_p_t["C",   M_t == "C" ]  <- 1 - p_NCD_C - p_CD - p_CCRE # cancer
    m_p_t["CRE",      M_t == "C" ]  <- p_CCRE #CRE (time since cancer dx to CRE)
    m_p_t["adult",    M_t == "C" ]  <- 0  # adult                                                
    m_p_t["dead",     M_t == "C" ]  <- p_NCD_C # dead from other causes
    m_p_t["candead", M_t == "C" ]<- p_CD # cancerdead after cancer dx
    
    #transition probabilities for those with CRE
    m_p_t["tobeborn",   M_t == "CRE" ] <- 0 # tobeborn
    m_p_t["noC",   M_t == "CRE" ] <- 0 # nocancer
    m_p_t["C",     M_t == "CRE" ] <- 0 # cancer
    m_p_t["CRE",        M_t == "CRE" ] <- 1 - p_NCD_CRE - p_CRED #CRE (time since cancer dx to CRE)
    m_p_t["adult",      M_t == "CRE" ] <- 0  # adult                                                  
    m_p_t["dead",       M_t == "CRE" ] <- p_NCD_CRE # dead from other causes
    m_p_t["candead", M_t == "CRE" ] <- p_CRED # cancerdead after cancer dx
    
    return(list("mpt"  = t(m_p_t), "df_X0"= df_X0))
  }       
  
  MicroSim <- function(n_i, df_X, seed = 1, kk) {
    
    set.seed(seed)
    
    m_M <-  matrix(nrow = n_i, ncol = n_t + 1,   # m_M stores health state information over time for every individual
                   dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                   paste("cycle", 0:n_t, sep = " ")))
    
    m_M[, 1] <- v_M_init         # initial health state (col 1)
    toAllocate <- c("noC", "C", "CRE") #the non-absorbing states
    
    # open a loop for time running cycles 1 to n_t
    for (t in 1:n_t) {
      
      # define which people to allocate through the Probs function and whom manually, "by hand"
      tobeallocated <- m_M[,t] %chin% toAllocate
      
      v_Allocate    <- m_M[,t] [tobeallocated]
      df_X_Allocate <-    df_X [tobeallocated,]
      
      
      ################################################################################ # 

      # calculate the transition probabilities for the cycle based on health state t
      P_all <- Probs(v_Allocate, df_X0 = df_X_Allocate, t, kk)
      m_P <- P_all$mpt
      
      #check_transition_probability(m_P, verbose = TRUE)  #check if transition probabilities are between 0 and 1
      # check if each of the rows of the transition probabilities matrix sum to one
      check_sum_of_transition_array(m_P, n_states = sum(tobeallocated), n_t = n_t, verbose = TRUE)
      # sample the current health state and store that state in matrix m_M
      m_M[tobeallocated == T, t + 1]  <- samplev(m_P) # sample next year's state for those that need to be allocated    
      m_M[tobeallocated == F, t + 1]  <- m_M[tobeallocated == F, t ]    
      df_X [tobeallocated,] <- P_all$df_X0
      
      # this "gives entry" to people who are born or immigrate in the next cycle (1-14)
      m_M[df_X$year_entry == (init_year + t) & ((init_year + t) - df_X$year_birth) <= max_dx_age, t + 1 ] = "noC" #if entering at age 0-14, enter nocancer state
      m_M[df_X$year_entry == (init_year + t) & ((init_year + t) - df_X$year_birth) > max_dx_age, t + 1 ] = "adult" #if entering at age 15+, enter adult state
      
      id_imm <- df_X$year_entry == (init_year + t) & df_X$year_birth < (init_year + t)
      
      # if born in the next cycle - assign an age of 0 and year_current
      df_X$age_current[df_X$year_entry == (init_year + t) & df_X$year_birth == (init_year + t)] = 0 #if your year_birth is equal to the year of the current microsim cycle, age is assigned to be 0
      df_X$year_current[df_X$year_entry == (init_year + t) & df_X$year_birth == (init_year + t)] = init_year + t
      
      # if immigrating in the next cycle - assign an age at migration and year_current
      df_X$age_current[id_imm] = df_X$year_entry[id_imm] - df_X$year_birth[id_imm]
      df_X$year_current[id_imm] = df_X$year_entry[id_imm]
      
      # this ages people who immigrated in the previous years
      df_X$age_current[df_X$year_entry <  (init_year + t) &
                         m_M[, t + 1] != "dead" &
                         m_M[, t + 1] != "candead" ] = 
        df_X$age_current[df_X$year_entry <  (init_year + t) &
                           m_M[, t + 1] != "dead" &
                           m_M[, t + 1] != "candead" ] + 1
      
      m_M[df_X$age_current >= (max_age-1), t+1] = "dead" #Automatically move those 100 and over into the death state
      
      m_M[df_X$year_imm == (init_year + t) & is.na(df_X$year_death), 1:t] = "outON" #replace m_M to indicate out of country for those who haven't officially immigrated yet.
      m_M[df_X$year_imm == (init_year + t) & !is.na(df_X$year_death), ] = "outON" #if youve died before you immigrated, assign your entire trajectory to be out of country
      
      # updates year_current for people who have been born or immigrated in the previous years
      df_X$year_current[df_X$year_entry < (init_year + t)] = df_X$year_current[df_X$year_entry < (init_year + t)] + 1
      
      ##### For those with cancer:
      # update time since cancer onset for t + 1
      df_X$v_dx0[ m_M[, t + 1] == "C"]   <- df_X$v_dx0[m_M [, t + 1] == "C"]   + 1 #time in cancer state
      df_X$v_cre0[m_M[, t + 1] == "CRE"] <- df_X$v_cre0[m_M[, t + 1] == "CRE"] + 1 #time in CRE state
      
      #for those who get diagnosed this cycle, assign age var to == age_current and year var to == dxyear for that cycle to retain diagnosis characteristics
      df_X$age[df_X$v_dx0 == 1 & m_M[, t + 1] == "C"] = df_X$age_current[df_X$v_dx0 == 1 & m_M[, t + 1] == "C"] #diagnosis age
      df_X$dxyear[df_X$v_dx0 == 1 & m_M[, t + 1] == "C"] = df_X$year_current[df_X$v_dx0 == 1 & m_M[, t + 1] == "C"] #diagnosis year
      # df_X$`I(log(dxyear))`[df_X$v_dx0 == 1 & m_M[, t + 1] == "C"] = log(df_X$year_current[df_X$v_dx0 == 1 & m_M[, t + 1] == "C"] - 1949) #unhide if using logdxyear MSM results
      
      ##### For those who die:
      #these vars update for people who died who did not have a cancer dx
      # assign time since death
      df_X$v_death0[m_M[, t + 1] == "dead" | m_M[, t + 1] == "candead"] <- df_X$v_death0[m_M[, t + 1] == "dead" | m_M[, t + 1] == "candead"] + 1 #time in death state
      
      # assign age at death to those who die in the next cycle
      df_X$age_death[df_X$v_death0 == 1 & m_M[, t + 1] == "dead" | df_X$v_death0 == 1 & m_M[, t + 1] == "candead"] = df_X$age_current[df_X$v_death0 == 1 & m_M[, t + 1] == "dead" | df_X$v_death0 == 1 & m_M[, t + 1] == "candead"] + 1
      
      # assign year of death
      df_X$year_death[df_X$v_death0 == 1] = df_X$year_current[df_X$v_death0 == 1]
      
      # assign age at CRE to those who are develop a CRE in the next cycle
      df_X$age_cre[df_X$v_cre0 == 1 & m_M[, t + 1] == "CRE"] = df_X$age_current[df_X$v_cre0 == 1 & m_M[, t + 1] == "CRE"] + 1
      
      # assign year of CRE
      df_X$year_cre[df_X$v_cre0 == 1] = df_X$year_current[df_X$v_cre0 == 1]
      
      #Assigning adults an annual risk of NCD to exit into the death state
      # if (t >= 1){
      # age_adult <- data.frame(table(df_X$age_current[m_M[, t + 1] == "adult"], df_X$sexMale[m_M[, t + 1] == "adult"]))
      # age_adult$year_current = init_year + t
      df_X1 <- df_X[m_M[, t + 1] == "adult",] #pull df_X char. for all who occupy adult state in next cycle
      age_adult2 <- left_join(df_X1, df_mort, by = c("sexMale", "year_current", "age_current"))[,"Rate"]
      NCD_adult <- rate_to_prob(age_adult2)
      isdead_adult <- rbinom(n = length(NCD_adult), size = 1, prob = NCD_adult)
      ID_dead <- df_X1$ID[isdead_adult == 1]
      m_M[ID_dead, t + 1] = "dead"
      # }
      
      ################################ #
      #Remove emigrants - females
      # length_emm_old <- 0
      # age_emm_remove <- vector("numeric")
      
      # if(t==6)browser()
      
      emm_year_f <- popON_net_f_remove[popON_net_f_remove$year == (init_year + t),] #females only, emigrants to remove, by age, per year
      emm_year_f <- emm_year_f[emm_year_f$age < 99,]
      emm_ages <- c(emm_year_f$age) #grab only the ages for immigrants that need to be removed
      
      for (e in emm_ages){
        # df_X_age2 <- df_X$ID[(((init_year + t) - df_X$year_birth) == (e)) & df_X$sexMale == 0]
        df_X_age2 <- df_X$ID[(((init_year + t) - df_X$year_birth) == (e)) & df_X$sexMale == 0 & m_M[, t + 1] != "dead" & m_M[, t + 1] != "candead" & m_M[, t + 1] != "outON" & m_M[, t + 1] != "leftON"]
        size_emm <- emm_year_f$sim_pop[emm_year_f$age == e] #grab the total # of emigrants for that specific age for sampling
        size_emm[length(df_X_age2) < size_emm] = length(df_X_age2) #if there are not enough people in the model at X age, sample from the limited number of them
        age_emm <- sample(df_X_age2, size_emm)
        # length_emm <- length(age_emm) + length_emm_old
        # age_emm_remove <- c(age_emm_remove, age_emm)
        
        m_M[age_emm, t+1] = "leftON" #replace m_M to indicate out of country for those who emigrated
        # m_M[df_X$ID %in% age_emm, t+1] = "outON" #replace m_M to indicate out of country for those who emigrated
        # df_X$year_emm[df_X$ID %in% age_emm] = init_year + t
        df_X$year_emm[age_emm] = init_year + t
      }
      # age_emm_remove_f <- setDT(age_emm_remove) #save the IDs of those who should have this replacement done (females)
      
      ################################ #
      #Remove emigrants - males
      # length_emm_old <- 0
      # age_emm_remove <- vector("numeric")
      
      emm_year_m <- popON_net_m_remove[popON_net_m_remove$year == (init_year + t),] #males only, emigrants to add, by age, per year
      emm_year_m <- emm_year_m[emm_year_m$age < 99,]
      emm_ages <- c(emm_year_m$age)
      
      for (e in emm_ages){
        # df_X_age2 <- df_X$ID[(((init_year + t) - df_X$year_birth) == (e)) & df_X$sexMale == 1]
        df_X_age2 <- df_X$ID[(((init_year + t) - df_X$year_birth) == (e)) & df_X$sexMale == 1 & m_M[, t + 1] != "dead" & m_M[, t + 1] != "candead" & m_M[, t + 1] != "outON" & m_M[, t + 1] != "leftON"]
        size_emm <- emm_year_m$sim_pop[emm_year_m$age == e]
        size_emm[length(df_X_age2) < size_emm] = length(df_X_age2)
        age_emm <- sample(df_X_age2, size_emm)
        # length_emm <- length(age_emm) + length_emm_old
        # age_emm_remove <- c(age_emm_remove, age_emm)
        
        m_M[age_emm, t+1] = "leftON"
        # m_M[df_X$ID %in% age_emm, t+1] = "outON"
        # df_X$year_emm[df_X$ID %in% age_emm] = init_year + t
        df_X$year_emm[age_emm] = init_year + t
      }
      # age_emm_remove_m <- setDT(age_emm_remove)
      ################################ #
      
      ##Display simulation progress
      if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
        cat('\r', paste(t/n_t * 100, "% done", sep = " "))
      }
      print(t)
      
    } # close the loop for the time points 
    
    
    # store the results from the simulation in a list
    results <- list(m_M = m_M, df_X = df_X)
    
    return(results)
    
  }
  
  outcomes <- MicroSim(n_i, df_X =df_X, seed = iteration*100, k)
  
  
  ################################ #
  #assign cancer types for 2040
  v_Cancer_new <- outcomes$m_M[,n_t+1] =="C" & outcomes$df_X$dxyear == 2040
  predrates_sumtemp <- predrates_sum
  predrates_sumtemp$age_current <- predrates_sumtemp$age_current +1
  rates_iccc <- inner_join(outcomes$df_X[v_Cancer_new,], predrates_sumtemp, by = c("sexMale", "year_current", "age_current"))%>% dplyr::select(all_of(ICCCgroups))
  outcomes$df_X[v_Cancer_new,"ICCCgroup"] <- samplev(rates_iccc)
  
  ################################ #
  
  outcomes$df_X$year_death[is.na(outcomes$df_X$year_death)] <- 9999
  outcomes$df_X$year_cre[is.na(outcomes$df_X$year_cre)] <- 9999
  outcomes$df_X$age_death[is.na(outcomes$df_X$age_death)] <- 9999
  outcomes$df_X$age_cre[is.na(outcomes$df_X$age_cre)] <- 9999
  outcomes$df_X$year_imm[is.na(outcomes$df_X$year_imm)] <- 1900
  
  # if (k == 1){
  #   saveRDS(outcomes, paste0("outcomes", "_", k, ".rds"))
  # }
  
  print(paste("b =", b, "k =", k, "iteration =", iteration, "finished")) 
  
  # fwrite(adj_fac_m, "adj_fac_m.csv", row.names=FALSE)
  
  
  ## Generate results from each iteration -------------------------------------------
  
  adj_fac <- adj_fac_m[1]
  
  #1970-2040 prevalence counts  (all cancers combined)
  prev_sim <- colSums(outcomes$m_M == "C" | outcomes$m_M == "CRE")*adj_fac
  
  #1970-2040 prevalence rates  (all cancers combined)
  prev_pop <- pop %>% group_by(year) %>% dplyr::summarise(pop = sum(pop)) %>% select(pop)
  prev_rate <- matrix((prev_sim / prev_pop$pop)*10^5)
  prevpop_PSA <- prev_pop$pop
  
  #1970-2040 prevalence - proportion of prevalent people that were diagnosed out of country  (all cancers combined)
  previmm <- which(outcomes$df_X$year_imm != 1900 & outcomes$df_X$dxyear < outcomes$df_X$year_imm) #not born in ON, and were diagnosed before entering ON
  m_M_prev_imm <- outcomes$m_M[previmm,]
  prev_imm <- colSums(m_M_prev_imm == "C" | m_M_prev_imm == "CRE", na.rm= T)*adj_fac
  
  #1970-2040 prevalence - proportion of prevalent people that were born and diagnosed in Ontario  (all cancers combined)
  prevbornON <- which(outcomes$df_X$year_imm == 1900) #born in Ontario
  m_M_prev_bornON <- outcomes$m_M[prevbornON,]
  prev_bornON <- colSums(m_M_prev_bornON == "C" | m_M_prev_bornON == "CRE", na.rm= T)*adj_fac
  
  #1970-2040 prevalence - proportion of prevalent people that were born outside Ontario but diagnosed in Ontario  (all cancers combined)
  previmm_dxON <- which(outcomes$df_X$year_imm != 1900 & outcomes$df_X$dxyear >= outcomes$df_X$year_imm) #not born in ON, were diagnosed in ON
  m_M_prev_imm_dxON <- outcomes$m_M[previmm_dxON,]
  prev_imm_dxON <- colSums(m_M_prev_imm_dxON == "C" | m_M_prev_imm_dxON == "CRE", na.rm= T)*adj_fac
  
  #1990-2019 ICES prevalence comparison (all cancers combined)
  true_prev_id <- which(outcomes$df_X$dxyear >= 1990 & outcomes$df_X$dxyear <= 2020 & outcomes$df_X$dxyear >= outcomes$df_X$year_imm)
  m_M_true_prev <- outcomes$m_M[true_prev_id,]
  prev_sim_true_prev <- colSums(m_M_true_prev == "C" | m_M_true_prev == "CRE", na.rm= T)*adj_fac
  
  #1970-2019 ICES prevalence comparison (all cancers combined)
  true_prev_id_70 <- which(outcomes$df_X$dxyear >= 1970 & outcomes$df_X$dxyear <= 2020 & outcomes$df_X$dxyear >= outcomes$df_X$year_imm)
  m_M_true_prev_70 <- outcomes$m_M[true_prev_id_70,]
  prev_sim_true_prev_70 <- colSums(m_M_true_prev_70 == "C" | m_M_true_prev_70 == "CRE", na.rm= T)*adj_fac
  
  #################################################### #
  #1990-2019 ICES prevalence comparison - by ICCC group
  
  tomodel1  <- outcomes$df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  
  emptydf <- expand.grid(ICCCgroup = ICCCgroups, year = ocr_pogo_years)
  
  microsim_trueprev_ICCC <- vector(mode = "list")
  for(i in 1:dim(tomodel1)[1]){
    
    true_prev_ICCC <- which(outcomes$df_X$dxyear >= 1990 & outcomes$df_X$dxyear <= 2020 & outcomes$df_X$ICCCgroup == tomodel1$ICCCgroup[i] & outcomes$df_X$dxyear >= outcomes$df_X$year_imm)
    m_M_true_prev_ICCC <- outcomes$m_M[true_prev_ICCC,]
    prev_sim_true_prev_byICCC <- colSums(m_M_true_prev_ICCC == "C" | m_M_true_prev_ICCC == "CRE", na.rm= T)*adj_fac
    
    prev_ICCC_df <- data.frame(prev_sim = prev_sim_true_prev_byICCC)
    prev_ICCC_df$ICCCgroup <- paste0(tomodel1$ICCCgroup[i])
    prev_ICCC_df$year <- rep(ocr_pogo_years)
    
    microsim_trueprev_ICCC[[i]] <- prev_ICCC_df
    #names(microsim_trueprev_ICCC)[i] <- paste0(tomodel1$ICCCgroup[i])
  }
  
  microsim_trueprev_ICCC_all <- ldply(microsim_trueprev_ICCC, data.frame)
  microsim_trueprev_ICCC_all_df <- left_join(microsim_trueprev_ICCC_all, emptydf, by = c("ICCCgroup", "year"))
  prev_sim_true_prev_ICCC <- microsim_trueprev_ICCC_all_df$prev_sim
  
  #################################################### #
  #1970-2019 ICES prevalence comparison - by ICCC group
  
  tomodel1  <- outcomes$df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  
  microsim_trueprev_70_ICCC <- vector(mode = "list")
  for(i in 1:dim(tomodel1)[1]){
    
    true_prev_70_ICCC <- which(outcomes$df_X$dxyear >= 1970 & outcomes$df_X$dxyear <= 2020 & outcomes$df_X$ICCCgroup == tomodel1$ICCCgroup[i] & outcomes$df_X$dxyear >= outcomes$df_X$year_imm)
    m_M_true_prev_70_ICCC <- outcomes$m_M[true_prev_70_ICCC,]
    prev_sim_true_prev_70_byICCC <- colSums(m_M_true_prev_70_ICCC == "C" | m_M_true_prev_70_ICCC == "CRE", na.rm= T)*adj_fac
    
    prev_70_ICCC_df <- data.frame(prev_sim = prev_sim_true_prev_70_byICCC)
    prev_70_ICCC_df$ICCCgroup <- paste0(tomodel1$ICCCgroup[i])
    prev_70_ICCC_df$year <- rep(ocr_pogo_years)
    
    microsim_trueprev_70_ICCC[[i]] <- prev_70_ICCC_df
    #names(microsim_trueprev_70_ICCC)[i] <- paste0(tomodel1$ICCCgroup[i])
  }
  
  microsim_trueprev_70_ICCC_all <- ldply(microsim_trueprev_70_ICCC, data.frame)
  microsim_trueprev_70_ICCC_all_df <- left_join(microsim_trueprev_70_ICCC_all, emptydf, by = c("ICCCgroup", "year"))
  prev_sim_true_prev_70_ICCC <- microsim_trueprev_70_ICCC_all_df$prev_sim
  
  
  #################################################### #
  
  #1970-2040 prevalence counts  (by ICCC group)
  
  tomodel1  <- outcomes$df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  
  microsimprev_ICCC <- vector(mode = "list")
  for(i in 1:dim(tomodel1)[1]){
    
    dfX_prevICCC <- which(outcomes$df_X$ICCCgroup == tomodel1$ICCCgroup[i])
    m_M_prevICCC <- outcomes$m_M[dfX_prevICCC,]
    prevICCC <- colSums(m_M_prevICCC == "C" | m_M_prevICCC == "CRE")*adj_fac
    
    prevICCC_df <- data.frame(prev_sim = prevICCC)
    prevICCC_df$ICCCgroup <- paste0(tomodel1$ICCCgroup[i])
    prevICCC_df$year <- rep(ocr_pogo_years)
    
    microsimprev_ICCC[[i]] <- prevICCC_df
    #names(microsimprev_ICCC)[i] <- paste0(tomodel1$ICCCgroup[i])
  }
  
  microsimprev_ICCC_all <- ldply(microsimprev_ICCC, data.frame)
  microsimprev_ICCC_all_df <- left_join(microsimprev_ICCC_all, emptydf, by = c("ICCCgroup", "year"))
  prev_sim_ICCC <- microsimprev_ICCC_all_df$prev_sim
  
  #1970-2040 prevalence rates (by ICCC group)
  
  prev_rate_ICCC <- matrix((prev_sim_ICCC / prev_pop$pop)*10^5)
  

  #################################################### #
  
  #1970-2040 prevalence counts  (by sex, decade of diagnosis - all cancers combined)
  
  cancer_id <- which(is.na(outcomes$df_X$dxyear )==F)
  cancer_df_X <- outcomes$df_X[cancer_id,]
  cancer_df_X$decade <- case_when(cancer_df_X$dxyear>=1970 & cancer_df_X$dxyear<1980 ~"70s",
                                  cancer_df_X$dxyear>=1980 & cancer_df_X$dxyear<1990 ~"80s",
                                  cancer_df_X$dxyear>=1990 & cancer_df_X$dxyear<2000 ~"90s",
                                  cancer_df_X$dxyear>=2000 & cancer_df_X$dxyear<2010 ~"00s",
                                  cancer_df_X$dxyear>=2010 & cancer_df_X$dxyear<2020 ~"10s",
                                  cancer_df_X$dxyear>=2020 & cancer_df_X$dxyear<2030 ~"20s",
                                  cancer_df_X$dxyear>=2030 & cancer_df_X$dxyear<=2040 ~"30s")
  cancer_df_X$decade <- factor(cancer_df_X$decade, levels = c("70s", "80s", "90s", "00s", "10s", "20s", "30s"),
                               labels = c("1970-79", "1980-89", "1990-99", "2000-09", "2010-19", "2020-29", "2030-39"))
  cancer_df_X$sex <- factor(cancer_df_X$sexMale, levels = c(0, 1), labels = c("Female", "Male"))
  cancer_m   <- outcomes$m_M[cancer_id ,]
  
  tomodel2 <- cancer_df_X %>% dplyr::count(decade, sex)
  emptydf_prev_char <- expand.grid(sex = sex_var_f,
                                   decade = decade_dx,
                                   year = ocr_pogo_years) #,
                                   #ICCCgroup = ICCCgroups)
  
  microsimprev_char <- vector(mode = "list")
  
  for(i in 1:dim(tomodel2)[1]){
    
    dfX_prev <- which(cancer_df_X$decade == tomodel2$decade[i] & cancer_df_X$sex == tomodel2$sex[i])
    m_M_prev <- cancer_m[dfX_prev,]
    prev <- colSums(m_M_prev == "C" | m_M_prev == "CRE", na.rm = T)*adj_fac
    
    prev_df <- data.frame(prev_sim = prev)
    
    prev_df$sex <- rep(tomodel2$sex[i])
    prev_df$year <- rep(ocr_pogo_years)
    prev_df$decade <- rep(tomodel2$decade[i])
    
    microsimprev_char[[i]] <- prev_df
  }
  
  microsimprev_char_all <- ldply(microsimprev_char, data.frame)
  microsimprev_char_all_df <- left_join(emptydf_prev_char, microsimprev_char_all, by = c("sex", "decade", "year"))
  microsimprev_char_all_df[is.na(microsimprev_char_all_df)] <- 0
  prev_sim_char <- microsimprev_char_all_df$prev_sim


#################################################### #
  #1970-2040 prevalence counts  (by ICCC group, sex, decade of diagnosis)
  
  tomodel2 <- cancer_df_X %>% dplyr::count(ICCCgroup, decade, sex) %>% filter(is.na(ICCCgroup) !=1)
  
  emptydf_prev_char <- expand.grid(sex = sex_var_f, 
                                   decade = decade_dx, 
                                   year = ocr_pogo_years, 
                                   ICCCgroup = ICCCgroups)
  
  microsimprev_ICCC_char <- vector(mode = "list")
  
  for(i in 1:dim(tomodel2)[1]){
   print(i)
    
    dfX_prevICCC <- which(cancer_df_X$ICCCgroup == tomodel2$ICCCgroup[i] & cancer_df_X$decade == tomodel2$decade[i] & cancer_df_X$sex == tomodel2$sex[i])
    m_M_prevICCC <- cancer_m[dfX_prevICCC,]
    
    if ((is.vector(m_M_prevICCC)==TRUE) | is.null(dim(m_M_prevICCC))) {

      m_M_prevICCC <- as.data.frame(t(m_M_prevICCC))
    
    } else{
      
      prevICCC <- colSums(m_M_prevICCC == "C" | m_M_prevICCC == "CRE", na.rm = T)*adj_fac  
    }
      
    prevICCC_df <- data.frame(prev_sim = prevICCC)
    
    prevICCC_df$ICCCgroup <- rep(tomodel2$ICCCgroup[i])
    prevICCC_df$sex <- rep(tomodel2$sex[i])
    prevICCC_df$year <- rep(ocr_pogo_years)
    prevICCC_df$decade <- rep(tomodel2$decade[i])

    microsimprev_ICCC_char[[i]] <- prevICCC_df
    #names(microsimprev_ICCC_char)[i] <- paste0(tomodel2$ICCCgroup[i],"_",tomodel2$decade[i],"_",tomodel2$sex[i])
  }

  microsimprev_ICCC_char_all <- ldply(microsimprev_ICCC_char, data.frame)
  microsimprev_ICCC_char_df <- left_join(emptydf_prev_char, microsimprev_ICCC_char_all, by = c("ICCCgroup", "sex", "decade", "year"))
  microsimprev_ICCC_char_df[is.na(microsimprev_ICCC_char_df)] <- 0
  prev_sim_ICCC_char <- microsimprev_ICCC_char_df$prev_sim
  
  
  #################################################### #
  
  #Overall incidence counts per year
  childpop <- pop %>% filter(age < 15) %>% group_by(year) %>% dplyr::summarise(pop = sum(pop)) %>% dplyr::rename(dxyear = year)
  childpop_PSA <- childpop$pop
  
  childpop_agesex <- pop %>% filter(age < 15) %>% group_by(year, age, sex) %>% dplyr::summarise(pop = sum(pop)) %>% dplyr::rename(dxyear = year)
  childpop_agesex_PSA <- childpop_agesex$pop
  
  dxyear_sim_imm_df <- outcomes$df_X %>%
    filter(dxyear >= year_imm) %>% #to look at incidence diagnosed in Ontario (among residents and immigrants)
    group_by(dxyear) %>%
    dplyr::count() %>%
    filter(is.na(dxyear) !=1) %>%
    mutate(dxyear = dxyear-1, n_adj = round(n * adj_fac))
  dxyear_sim_imm_df <- left_join(dxyear_sim_imm_df, childpop) %>% mutate(rate = ((n_adj/pop)*10^6))
  dxyear_sim_imm <- dxyear_sim_imm_df$n_adj
  dxyear_sim_imm_rate <- dxyear_sim_imm_df$rate
  
  dxyear_sim_imm5 <- dxyear_sim_imm_df %>% filter(dxyear >= start_year_pred) %>%
    dplyr::mutate(fiveyrperiod = case_when(dxyear >= 2020 & dxyear <= 2024 ~ '2020-2024',
                                           dxyear >= 2025 & dxyear <= 2029 ~ '2025-2029',
                                           dxyear >= 2030 & dxyear <= 2034 ~ '2030-2034',
                                           dxyear >= 2035 & dxyear <= 2039 ~ '2035-2039'))
  
  dxyear_sim_imm5 <- aggregate(n_adj ~ fiveyrperiod, data = dxyear_sim_imm5, FUN = "sum")
  dxyear_sim_imm5_df <- dxyear_sim_imm5$n_adj
  
  childpop_sum <- pop %>% filter(age < 15) %>% group_by(year) %>% dplyr::summarise(pop = sum(pop)) %>% dplyr::rename(dxyear = year) %>% filter(dxyear >= start_year_pred) %>%
    dplyr::mutate(fiveyrperiod = case_when(dxyear >= 2020 & dxyear <= 2024 ~ '2020-2024',
                                           dxyear >= 2025 & dxyear <= 2029 ~ '2025-2029',
                                           dxyear >= 2030 & dxyear <= 2034 ~ '2030-2034',
                                           dxyear >= 2035 & dxyear <= 2039 ~ '2035-2039'))
  childpop_sum <- aggregate(pop ~ fiveyrperiod, data = childpop_sum, FUN = "sum")
  
  dxyear_sim_imm5_rate <- left_join(dxyear_sim_imm5, childpop_sum, by = c("fiveyrperiod")) %>% mutate(rate = (n_adj/pop)*10^6)
  dxyear_sim_imm5_rate_df <- dxyear_sim_imm5_rate$rate
  
  #Incidence counts per year by cancer type
  dxyeargroup_sim_df <- outcomes$df_X %>%
    filter(is.na(ICCCgroup) !=1) %>%
    filter(dxyear >= year_imm) %>%
    group_by(dxyear, ICCCgroup) %>%
    dplyr::count() %>%
    mutate(dxyear = dxyear-1, n_adj = n * adj_fac)
  #merge with df that includes all years/iccc groups
  emptydf <- data.frame(dxyear = rep(ocr_pogo_years, each = length(ICCCgroups)),
                        ICCCgroup = rep(ICCCgroups, times = length(ocr_pogo_years)),
                        n = as.numeric(rep(NA)))
  emptydf <- filter(emptydf, dxyear >= min(dxyeargroup_sim_df$dxyear, na.rm = TRUE) & dxyear <= max(dxyeargroup_sim_df$dxyear, na.rm = TRUE))
  dxyeargroup_sim_df <- left_join(emptydf, dxyeargroup_sim_df, by = c("ICCCgroup", "dxyear")) %>% select(!n.x) %>% dplyr::rename(n = n.y)
  dxyeargroup_sim_df[is.na(dxyeargroup_sim_df)] <- 0
  dxyeargroup_sim <- dxyeargroup_sim_df$n_adj
  
  #Incidence counts per 5year time periods by cancer type
  dxyeargroup_iccc <- dxyeargroup_sim_df %>% filter(dxyear >= start_year_pred) %>%
    mutate(fiveyrperiod = case_when(dxyear >= 2020 & dxyear <= 2024 ~ '2020-2024',
                                    dxyear >= 2025 & dxyear <= 2029 ~ '2025-2029',
                                    dxyear >= 2030 & dxyear <= 2034 ~ '2030-2034',
                                    dxyear >= 2035 & dxyear <= 2039 ~ '2035-2039'))
  dxyeargroup_iccc <- aggregate(n_adj ~ fiveyrperiod + ICCCgroup, data = dxyeargroup_iccc, FUN = "sum")
  dxyeargroup_iccc <- left_join(dxyeargroup_iccc, childpop_sum, by = c("fiveyrperiod")) %>% mutate(n_avg = n_adj/5, rate = round((n_adj/pop)*10^6, 1)) %>% select(!pop)
  dxyeargroup_iccc_df <- dxyeargroup_iccc$n_adj
  dxyeargroup_iccc_rate_df <- dxyeargroup_iccc$rate
  
  #Incidence counts per year by cancer type, sex, age at diagnosis
  dxyeargroup_sim_age_df <- outcomes$df_X %>%
    filter(is.na(ICCCgroup) !=1) %>%
    filter(dxyear >= year_imm) %>%
    group_by(dxyear, ICCCgroup, sexMale, age) %>%
    dplyr::count() %>%
    mutate(dxyear = dxyear-1, age = age-1, n_adj = n * adj_fac)

  # Create the dataframe with all combinations
  emptydf_inc_agesex <- expand.grid(sexMale = sex_var, age = dxages, dxyear = ocr_pogo_years, ICCCgroup = ICCCgroups)
  emptydf_inc_agesex$n <- as.numeric(rep(NA))
  emptydf_inc_agesex <- filter(emptydf_inc_agesex, dxyear <= 2039) #remove yr2040 for incidence reporting
  
  dxyeargroup_sim_age_df <- left_join(emptydf_inc_agesex, dxyeargroup_sim_age_df, by = c("ICCCgroup", "dxyear", "sexMale", "age")) %>% select(!n.x) %>% dplyr::rename(n = n.y)
  dxyeargroup_sim_age_df[is.na(dxyeargroup_sim_age_df)] <- 0
  dxyeargroup_sim_age <- dxyeargroup_sim_age_df$n_adj

  ############ #
  #overall prevalence by attained age for 2020 and 2040
  prev_ICCC_age <- which(outcomes$df_X$dxyear > init_year & outcomes$df_X$dxyear <= max_year)
  m_M_prev_ICCC_age <- outcomes$m_M[prev_ICCC_age,]
  m_M_prev_ICCC_age[m_M_prev_ICCC_age == "tobeborn"] = 0
  m_M_prev_ICCC_age[m_M_prev_ICCC_age == "noC"] = 0
  m_M_prev_ICCC_age[m_M_prev_ICCC_age == "outON"] = 0
  m_M_prev_ICCC_age[m_M_prev_ICCC_age == "leftON"] = NA
  m_M_prev_ICCC_age[m_M_prev_ICCC_age == "candead"] = NA
  m_M_prev_ICCC_age[m_M_prev_ICCC_age == "dead"] = NA
  attainedage <- apply(m_M_prev_ICCC_age, 1, function(x){ y = x != "0"; cumsum(y)})
  attainedage <- t(attainedage)
  attainedage[attainedage == 0] <- NA
  attainedage2 <- outcomes$df_X[prev_ICCC_age,]
  colnames(attainedage) <- init_year:max_year
  
  attainedage2$age2 <- attainedage2$age #age at diagnosis for all
  attainedage2$age2[attainedage2$dxyear < attainedage2$year_imm & attainedage2$year_imm != 1900] = (attainedage2$year_imm - attainedage2$year_birth)[attainedage2$dxyear < attainedage2$year_imm & attainedage2$year_imm != 1900] #for those who had a diagnosis before immigrating, age2 variable displays their age at immigration
  
  attainedage <- attainedage + attainedage2$age2 - 1 #need the minus 1 to account for individuals starting with age 0 not being reflected in the attainedage matrix
  
  ############ #
  
  attainedage_group <- c("0-19", "20-39", "40-59", "60+")
  emptydf_prev <- data.frame(attainedage = rep(attainedage_group))
  
  age <- attainedage[,(length(v_years) - 20)] #2020 is col 51
  table <-  table(age)
  table <- data.table(table)
  table$age <- as.numeric(table$age)
  table <- table %>% mutate(attainedage = case_when(age <= 19 ~ "0-19",
                                                    age >= 20 & age <= 39 ~ "20-39",
                                                    age >= 40 & age <= 59 ~ "40-59",
                                                    age >= 60 ~ "60+"))
  table <- table %>% group_by(attainedage) %>% dplyr::summarise(n = sum(N))
  table$n <- round(table$n*adj_fac)
  table <- left_join(emptydf_prev, table, by = c("attainedage"))
  table$year <- rep(2020)
  table[is.na(table)] <- 0
  table <- data.table(table)
  table2020 <- table
  
  age <- attainedage[,(length(v_years) - 0)] #2040 is col 71
  table <-  table(age)
  table <- data.table(table)
  table$age <- as.numeric(table$age)
  table <- table %>% mutate(attainedage = case_when(age <= 19 ~ "0-19",
                                                    age >= 20 & age <= 39 ~ "20-39",
                                                    age >= 40 & age <= 59 ~ "40-59",
                                                    age >= 60 ~ "60+"))
  table <- table %>% group_by(attainedage) %>% dplyr::summarise(n = sum(N))
  table$n <- round(table$n*adj_fac)
  table <- left_join(emptydf_prev, table, by = c("attainedage"))
  table$year <- rep(2040)
  table[is.na(table)] <- 0
  table <- data.table(table)
  table2040 <- table
  
  attainedage_year_df <- rbind(table2020, table2040)
  attainedage_year_df_m <- attainedage_year_df$n
  
  #prevalence by attained age and cancer type for 2020 and 2040
  tomodel1  <- outcomes$df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  attainedage_ICCC2020 <- data.frame(attainedage_group = NA, n = NA, percent = NA, year = NA, ICCCgroup = NA)
  attainedage_ICCC2040 <- data.frame(attainedage_group = NA, n = NA, percent = NA, year = NA, ICCCgroup = NA)
  emptydf_prev <- data.frame(attainedage_group = rep(attainedage_group))
  
  for(i in 1:dim(tomodel1)[1]){
    
    #2020
    attainedage_ICCC_2020 <- attainedage[attainedage2$ICCCgroup == tomodel1$ICCCgroup[i] & attainedage2$dxyear <= 2020,]
    attainedage_ICCC_2020 <- data.table(attainedage_ICCC_2020)
    attainedage_ICCC_2020 <- attainedage_ICCC_2020 %>% select("2020")
    colnames(attainedage_ICCC_2020)[1] <- "attainedage"
    attainedage_ICCC_2020 <- na.omit(attainedage_ICCC_2020)
    attainedage_ICCC_2020 <- attainedage_ICCC_2020 %>% mutate(attainedage_group = case_when(attainedage <= 19 ~ "0-19",
                                                                                            attainedage >= 20 & attainedage <= 39 ~ "20-39",
                                                                                            attainedage >= 40 & attainedage <= 59 ~ "40-59",
                                                                                            attainedage >= 60 ~ "60+"))
    attainedage_ICCC_2020 <- attainedage_ICCC_2020 %>% group_by(attainedage_group) %>% dplyr::count()
    attainedage_ICCC_2020$n <- round(attainedage_ICCC_2020$n*adj_fac)
    attainedage_ICCC_2020$percent <- round((attainedage_ICCC_2020$n / sum(attainedage_ICCC_2020$n))*100,2)
    attainedage_ICCC_2020 <- left_join(emptydf_prev, attainedage_ICCC_2020, by = c("attainedage_group"))
    attainedage_ICCC_2020$year <- rep(2020)
    attainedage_ICCC_2020$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    attainedage_ICCC_2020 <- data.table(attainedage_ICCC_2020)
    attainedage_ICCC2020 <- rbind(attainedage_ICCC2020, attainedage_ICCC_2020)
    
    #2040
    attainedage_ICCC_2040 <- attainedage[attainedage2$ICCCgroup == tomodel1$ICCCgroup[i] & attainedage2$dxyear <= 2040,]
    attainedage_ICCC_2040 <- data.table(attainedage_ICCC_2040)
    attainedage_ICCC_2040 <- attainedage_ICCC_2040 %>% select("2040")
    colnames(attainedage_ICCC_2040)[1] <- "attainedage"
    attainedage_ICCC_2040 <- na.omit(attainedage_ICCC_2040)
    attainedage_ICCC_2040 <- attainedage_ICCC_2040 %>% mutate(attainedage_group = case_when(attainedage <= 19 ~ "0-19",
                                                                                            attainedage >= 20 & attainedage <= 39 ~ "20-39",
                                                                                            attainedage >= 40 & attainedage <= 59 ~ "40-59",
                                                                                            attainedage >= 60 ~ "60+"))
    attainedage_ICCC_2040 <- attainedage_ICCC_2040 %>% group_by(attainedage_group) %>% dplyr::count()
    attainedage_ICCC_2040$n <- round(attainedage_ICCC_2040$n*adj_fac)
    attainedage_ICCC_2040$percent <- round((attainedage_ICCC_2040$n / sum(attainedage_ICCC_2040$n))*100,2)
    attainedage_ICCC_2040 <- left_join(emptydf_prev, attainedage_ICCC_2040, by = c("attainedage_group"))
    attainedage_ICCC_2040$year <- rep(2040)
    attainedage_ICCC_2040$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    attainedage_ICCC_2040 <- data.table(attainedage_ICCC_2040)
    attainedage_ICCC2040 <- rbind(attainedage_ICCC2040, attainedage_ICCC_2040)
  }
  attainedage_ICCC <- rbind(attainedage_ICCC2020, attainedage_ICCC2040)
  attainedage_ICCC <- filter(attainedage_ICCC, rowSums(is.na(attainedage_ICCC)) != ncol(attainedage_ICCC))
  attainedage_ICCC[is.na(attainedage_ICCC)] <- 0
  #attainedage_ICCC_m[,k] <- attainedage_ICCC$n
  attainedage_ICCC_m <- attainedage_ICCC$n
  
  ########################################################################################################## #
  
  #overall prevalence by attained age group for ASPR calculation
 
  attainedage_group2 <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
                          "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
  
  emptydf_prev_aspr <- expand.grid(attainedage = attainedage_group2)
  
  age <- attainedage[,(length(v_years) - 20)] #2020 is col 51
  table <-  table(age)
  table <- data.table(table)
  table$age <- as.numeric(table$age)
  table <- table %>% mutate(attainedage = case_when(age <= 4 ~ "0-4",
                                                    
                                                    age >= 5 & age <= 9 ~ "5-9",
                                                    age >= 10 & age <= 14 ~ "10-14",
                                                    
                                                    age >= 15 & age <= 19 ~ "15-19",
                                                    age >= 20 & age <= 24 ~ "20-24",
                                                    
                                                    age >= 25 & age <= 29 ~ "25-29",
                                                    age >= 30 & age <= 34 ~ "30-34",
                                                    
                                                    age >= 35 & age <= 39 ~ "35-39",
                                                    age >= 40 & age <= 44 ~ "40-44",
                                                    
                                                    age >= 45 & age <= 49 ~ "45-49",
                                                    age >= 50 & age <= 54 ~ "50-54",
                                                    
                                                    age >= 55 & age <= 59 ~ "55-59",
                                                    age >= 60 & age <= 64 ~ "60-64",
                                                    
                                                    age >= 65 & age <= 69 ~ "65-69",
                                                    age >= 70 & age <= 74 ~ "70-74",
                                                    
                                                    age >= 75 & age <= 79 ~ "75-79",
                                                    age >= 80 & age <= 84 ~ "80-84",
  
                                                    age >= 85 ~ "85+"))
  table <- table %>% group_by(attainedage) %>% dplyr::summarise(n = sum(N))
  table$n <- round(table$n*adj_fac)
  table <- left_join(emptydf_prev_aspr, table, by = c("attainedage"))
  table$year <- rep(2020)
  table[is.na(table)] <- 0
  table <- data.table(table)
  table2020 <- table
  
  age <- attainedage[,(length(v_years) - 0)] #2040 is col 71
  table <-  table(age)
  table <- data.table(table)
  table$age <- as.numeric(table$age)
  table <- table %>% mutate(attainedage = case_when(age <= 4 ~ "0-4",
                                                    
                                                    age >= 5 & age <= 9 ~ "5-9",
                                                    age >= 10 & age <= 14 ~ "10-14",
                                                    
                                                    age >= 15 & age <= 19 ~ "15-19",
                                                    age >= 20 & age <= 24 ~ "20-24",
                                                    
                                                    age >= 25 & age <= 29 ~ "25-29",
                                                    age >= 30 & age <= 34 ~ "30-34",
                                                    
                                                    age >= 35 & age <= 39 ~ "35-39",
                                                    age >= 40 & age <= 44 ~ "40-44",
                                                    
                                                    age >= 45 & age <= 49 ~ "45-49",
                                                    age >= 50 & age <= 54 ~ "50-54",
                                                    
                                                    age >= 55 & age <= 59 ~ "55-59",
                                                    age >= 60 & age <= 64 ~ "60-64",
                                                    
                                                    age >= 65 & age <= 69 ~ "65-69",
                                                    age >= 70 & age <= 74 ~ "70-74",
                                                    
                                                    age >= 75 & age <= 79 ~ "75-79",
                                                    age >= 80 & age <= 84 ~ "80-84",
                                                    
                                                    age >= 85 ~ "85+"))
  table <- table %>% group_by(attainedage) %>% dplyr::summarise(n = sum(N))
  table$n <- round(table$n*adj_fac)
  table <- left_join(emptydf_prev_aspr, table, by = c("attainedage"))
  table$year <- rep(2040)
  table[is.na(table)] <- 0
  table <- data.table(table)
  table2040 <- table
  
  attainedage_aspr_df <- rbind(table2020, table2040)
  attainedage_aspr_df_m <- attainedage_aspr_df$n 
  
  
  
  ########################################################################################################## #
  
  #prevalence by attained age and cancer type for ASPR calculation - 2020 and 2040
  tomodel1  <- outcomes$df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  attainedage_ICCC2020 <- data.frame(attainedage_group = NA, n = NA,  year = NA, ICCCgroup = NA)
  attainedage_ICCC2040 <- data.frame(attainedage_group = NA, n = NA,  year = NA, ICCCgroup = NA)
  emptydf_prev_aspr <- data.frame(attainedage_group = rep(attainedage_group2))
  
  for(i in 1:dim(tomodel1)[1]){
    
    #2020
    attainedage_ICCC_2020 <- attainedage[attainedage2$ICCCgroup == tomodel1$ICCCgroup[i] & attainedage2$dxyear <= 2020,]
    attainedage_ICCC_2020 <- data.table(attainedage_ICCC_2020)
    attainedage_ICCC_2020 <- attainedage_ICCC_2020 %>% select("2020")
    colnames(attainedage_ICCC_2020)[1] <- "attainedage"
    attainedage_ICCC_2020 <- na.omit(attainedage_ICCC_2020)
    attainedage_ICCC_2020 <- attainedage_ICCC_2020 %>% mutate(attainedage_group = case_when(attainedage <= 4 ~ "0-4",
                                                      
                                                      attainedage >= 5 & attainedage <= 9 ~ "5-9",
                                                      attainedage >= 10 & attainedage <= 14 ~ "10-14",
                                                      
                                                      attainedage >= 15 & attainedage <= 19 ~ "15-19",
                                                      attainedage >= 20 & attainedage <= 24 ~ "20-24",
                                                      
                                                      attainedage >= 25 & attainedage <= 29 ~ "25-29",
                                                      attainedage >= 30 & attainedage <= 34 ~ "30-34",
                                                      
                                                      attainedage >= 35 & attainedage <= 39 ~ "35-39",
                                                      attainedage >= 40 & attainedage <= 44 ~ "40-44",
                                                      
                                                      attainedage >= 45 & attainedage <= 49 ~ "45-49",
                                                      attainedage >= 50 & attainedage <= 54 ~ "50-54",
                                                      
                                                      attainedage >= 55 & attainedage <= 59 ~ "55-59",
                                                      attainedage >= 60 & attainedage <= 64 ~ "60-64",
                                                      
                                                      attainedage >= 65 & attainedage <= 69 ~ "65-69",
                                                      attainedage >= 70 & attainedage <= 74 ~ "70-74",
                                                      
                                                      attainedage >= 75 & attainedage <= 79 ~ "75-79",
                                                      attainedage >= 80 & attainedage <= 84 ~ "80-84",
                                                      
                                                      attainedage >= 85 ~ "85+"))
    
    attainedage_ICCC_2020 <- attainedage_ICCC_2020 %>% group_by(attainedage_group) %>% dplyr::count()
    attainedage_ICCC_2020$n <- round(attainedage_ICCC_2020$n*adj_fac)
    #attainedage_ICCC_2020$percent <- round((attainedage_ICCC_2020$n / sum(attainedage_ICCC_2020$n))*100,2)
    attainedage_ICCC_2020 <- left_join(emptydf_prev_aspr, attainedage_ICCC_2020, by = c("attainedage_group"))
    attainedage_ICCC_2020$year <- rep(2020)
    attainedage_ICCC_2020$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    attainedage_ICCC_2020 <- data.table(attainedage_ICCC_2020)
    attainedage_ICCC2020 <- rbind(attainedage_ICCC2020, attainedage_ICCC_2020)
    
    #2040
    attainedage_ICCC_2040 <- attainedage[attainedage2$ICCCgroup == tomodel1$ICCCgroup[i] & attainedage2$dxyear <= 2040,]
    attainedage_ICCC_2040 <- data.table(attainedage_ICCC_2040)
    attainedage_ICCC_2040 <- attainedage_ICCC_2040 %>% select("2040")
    colnames(attainedage_ICCC_2040)[1] <- "attainedage"
    attainedage_ICCC_2040 <- na.omit(attainedage_ICCC_2040)
    attainedage_ICCC_2040 <- attainedage_ICCC_2040 %>% mutate(attainedage_group = case_when(attainedage <= 4 ~ "0-4",
                                                                                            
                                                                                            attainedage >= 5 & attainedage <= 9 ~ "5-9",
                                                                                            attainedage >= 10 & attainedage <= 14 ~ "10-14",
                                                                                            
                                                                                            attainedage >= 15 & attainedage <= 19 ~ "15-19",
                                                                                            attainedage >= 20 & attainedage <= 24 ~ "20-24",
                                                                                            
                                                                                            attainedage >= 25 & attainedage <= 29 ~ "25-29",
                                                                                            attainedage >= 30 & attainedage <= 34 ~ "30-34",
                                                                                            
                                                                                            attainedage >= 35 & attainedage <= 39 ~ "35-39",
                                                                                            attainedage >= 40 & attainedage <= 44 ~ "40-44",
                                                                                            
                                                                                            attainedage >= 45 & attainedage <= 49 ~ "45-49",
                                                                                            attainedage >= 50 & attainedage <= 54 ~ "50-54",
                                                                                            
                                                                                            attainedage >= 55 & attainedage <= 59 ~ "55-59",
                                                                                            attainedage >= 60 & attainedage <= 64 ~ "60-64",
                                                                                            
                                                                                            attainedage >= 65 & attainedage <= 69 ~ "65-69",
                                                                                            attainedage >= 70 & attainedage <= 74 ~ "70-74",
                                                                                            
                                                                                            attainedage >= 75 & attainedage <= 79 ~ "75-79",
                                                                                            attainedage >= 80 & attainedage <= 84 ~ "80-84",
                                                                                            
                                                                                            attainedage >= 85 ~ "85+"))
    attainedage_ICCC_2040 <- attainedage_ICCC_2040 %>% group_by(attainedage_group) %>% dplyr::count()
    attainedage_ICCC_2040$n <- round(attainedage_ICCC_2040$n*adj_fac)
    #attainedage_ICCC_2040$percent <- round((attainedage_ICCC_2040$n / sum(attainedage_ICCC_2040$n))*100,2)
    attainedage_ICCC_2040 <- left_join(emptydf_prev_aspr, attainedage_ICCC_2040, by = c("attainedage_group"))
    attainedage_ICCC_2040$year <- rep(2040)
    attainedage_ICCC_2040$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    attainedage_ICCC_2040 <- data.table(attainedage_ICCC_2040)
    attainedage_ICCC2040 <- rbind(attainedage_ICCC2040, attainedage_ICCC_2040)
  }
  attainedage_ICCC <- rbind(attainedage_ICCC2020, attainedage_ICCC2040)
  attainedage_ICCC <- filter(attainedage_ICCC, rowSums(is.na(attainedage_ICCC)) != ncol(attainedage_ICCC))
  attainedage_ICCC[is.na(attainedage_ICCC)] <- 0
  attainedage_ICCC_aspr_m <- attainedage_ICCC$n
  
################################################################################ #  
  
  #counting prevalence by X-year survivors (e.g. 5, 10, 20...) - time since diagnosis in years
  surv_status <- which(outcomes$df_X$dxyear > init_year & outcomes$df_X$dxyear <= max_year)
  m_M_surv_status <- outcomes$m_M[surv_status,]
  m_M_surv_status[m_M_surv_status == "tobeborn"] = 0
  m_M_surv_status[m_M_surv_status == "noC"] = 0
  m_M_surv_status[m_M_surv_status == "outON"] = 0
  m_M_surv_status[m_M_surv_status == "leftON"] = NA
  m_M_surv_status[m_M_surv_status == "candead"] = NA
  m_M_surv_status[m_M_surv_status == "dead"] = NA
  survstatus <- apply(m_M_surv_status, 1, function(x){ y = x != "0"; cumsum(y)})
  survstatus <- t(survstatus)
  survstatus[survstatus == 0] <- NA
  survstatus2 <- outcomes$df_X[surv_status,]
  colnames(survstatus) <- init_year:max_year
  
  survstatus2$survtime <- 0 #time from diagnosis placeholder var.
  survstatus2$survtime[survstatus2$dxyear < survstatus2$year_imm & survstatus2$year_imm != 1900] = (survstatus2$year_imm - survstatus2$dxyear)[survstatus2$dxyear < survstatus2$year_imm & survstatus2$year_imm != 1900] #for those who had a diagnosis before immigrating, survtime variable displays their time since diagnosis in yrs prior to entering ON (which then gets added on in the next line)
  
  survstatus <- survstatus + survstatus2$survtime - 1 #people with survstatus of 0 have been diagnosed in that year
  
  survstatus_group <- c("<5", "5-9", "10-14", "15-19", "20+")
  emptydf_surv <- data.frame(survstatus = rep(survstatus_group))
  
  survstatus_group2 <- c("0-1 years", "2 years", "3 years", "4 years", "5+ years")
  emptydf_surv2 <- data.frame(survstatus2 = rep(survstatus_group2))
  
  survtable <- survstatus[,(length(v_years) - 20)] #2020 is col 51
  table <-  table(survtable)
  table <- data.table(table)
  table$surv <- as.numeric(table$surv)
  table$surv <- table$surv - 1 #changing from factor to numeric var above mistakenly changed range of values from 0-50 to 1-51
  table_all <- table %>% mutate(survstatus = case_when(surv < 5 ~ "<5",
                                                       surv >= 5 & surv <= 9 ~ "5-9",
                                                       surv >= 10 & surv <= 14 ~ "10-14",
                                                       surv >= 15 & surv <= 19 ~ "15-19",
                                                       surv >= 20 ~ "20+"),
                                survstatus2 = case_when (surv <= 1 ~ "0-1 years",
                                                         surv <= 2 ~ "2 years",
                                                         surv <= 3 ~ "3 years",
                                                         surv <= 4 ~ "4 years",
                                                         surv >= 5 ~ "5+ years"))
  #FU in 1 yr periods until 5 yrs
  table <- table_all %>% group_by(survstatus) %>% dplyr::summarise(n = sum(N))
  table$n <- (table$n*adj_fac)
  table <- left_join(emptydf_surv, table, by = c("survstatus"))
  table$year <- rep(2020)
  table[is.na(table)] <- 0
  table <- data.table(table)
  table2020 <- table
  
  #FU in 1 yr periods until 5 yrs
  table2 <- table_all %>% group_by(survstatus2) %>% dplyr::summarise(n = sum(N))
  table2$n <- (table2$n*adj_fac)
  table2 <- left_join(emptydf_surv2, table2, by = c("survstatus2"))
  table2$year <- rep(2020)
  table2[is.na(table2)] <- 0
  table2 <- data.table(table2)
  table22020 <- table2
  
  emptydf_surv <- data.frame(survstatus = rep(survstatus_group))
  emptydf_surv2 <- data.frame(survstatus2 = rep(survstatus_group2))
  
  survtable <- survstatus[,(length(v_years) - 0)] #2040 is col 71
  table <-  table(survtable)
  table <- data.table(table)
  table$surv <- as.numeric(table$surv)
  table$surv <- table$surv - 1 #changing from factor to numeric var above mistakenly changed range of values from 0-50 to 1-51
  table_all <- table %>% mutate(survstatus = case_when(surv < 5 ~ "<5",
                                                       surv >= 5 & surv <= 9 ~ "5-9",
                                                       surv >= 10 & surv <= 14 ~ "10-14",
                                                       surv >= 15 & surv <= 19 ~ "15-19",
                                                       surv >= 20 ~ "20+"),
                                survstatus2 = case_when (surv <= 1 ~ "0-1 years",
                                                         surv <= 2 ~ "2 years",
                                                         surv <= 3 ~ "3 years",
                                                         surv <= 4 ~ "4 years",
                                                         surv >= 5 ~ "5+ years"))
  #FU in 5 year periods
  table <- table_all %>% group_by(survstatus) %>% dplyr::summarise(n = sum(N))
  table$n <- (table$n*adj_fac)
  table <- left_join(emptydf_surv, table, by = c("survstatus"))
  table$year <- rep(2040)
  table[is.na(table)] <- 0
  table <- data.table(table)
  table2040 <- table
  
  survstatusdf <- rbind(table2020, table2040)
  survstatus_df <- survstatusdf$n
  
  #FU in 1 yr periods until 5 yrs
  table2 <- table_all %>% group_by(survstatus2) %>% dplyr::summarise(n = sum(N))
  table2$n <- (table2$n*adj_fac)
  table2 <- left_join(emptydf_surv2, table2, by = c("survstatus2"))
  table2$year <- rep(2040)
  table2[is.na(table2)] <- 0
  table2 <- data.table(table2)
  table22040 <- table2
  
  survstatusdf2 <- rbind(table22020, table22040)
  survstatus_df2 <- survstatusdf2$n
  
  #prevalence by surv status and cancer type for 2020 and 2040 - 5yr periods
  tomodel1  <- outcomes$df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  survstatus_ICCC2020 <- data.frame(survstatus_group = NA, n = NA, percent = NA, year = NA, ICCCgroup = NA)
  survstatus_ICCC2040 <- data.frame(survstatus_group = NA, n = NA, percent = NA, year = NA, ICCCgroup = NA)
  emptydf_surv <- data.frame(survstatus_group = rep(survstatus_group))
  
  survstatus_ICCC20202 <- data.frame(survstatus_group2 = NA, n = NA, percent = NA, year = NA, ICCCgroup = NA)
  survstatus_ICCC20402 <- data.frame(survstatus_group2 = NA, n = NA, percent = NA, year = NA, ICCCgroup = NA)
  emptydf_surv2 <- data.frame(survstatus_group2 = rep(survstatus_group2))
  
  for(i in 1:dim(tomodel1)[1]){
    
    ###2020
    survstatus_ICCC_2020 <- survstatus[survstatus2$ICCCgroup == tomodel1$ICCCgroup[i] & survstatus2$dxyear <= 2020,]
    survstatus_ICCC_2020 <- data.table(survstatus_ICCC_2020)
    survstatus_ICCC_2020 <- survstatus_ICCC_2020 %>% select("2020")
    colnames(survstatus_ICCC_2020)[1] <- "survstatus"
    survstatus_ICCC_2020 <- na.omit(survstatus_ICCC_2020)
    survstatus_ICCC_2020_all <- survstatus_ICCC_2020 %>% mutate(survstatus_group = case_when(survstatus < 5 ~ "<5",
                                                                                             survstatus >= 5 & survstatus <= 9 ~ "5-9",
                                                                                             survstatus >= 10 & survstatus <= 14 ~ "10-14",
                                                                                             survstatus >= 15 & survstatus <= 19 ~ "15-19",
                                                                                             survstatus >= 20 ~ "20+"),
                                                                
                                                                survstatus_group2 = case_when (survstatus <= 1 ~ "0-1 years",
                                                                                               survstatus <= 2 ~ "2 years",
                                                                                               survstatus <= 3 ~ "3 years",
                                                                                               survstatus <= 4 ~ "4 years",
                                                                                               survstatus >= 5 ~ "5+ years"))
    
    #2020 - #FU in 5 year periods
    survstatus_ICCC_2020 <- survstatus_ICCC_2020_all %>% group_by(survstatus_group) %>% dplyr::count()
    survstatus_ICCC_2020$n <- (survstatus_ICCC_2020$n)*adj_fac
    survstatus_ICCC_2020$percent <- round((survstatus_ICCC_2020$n / sum(survstatus_ICCC_2020$n))*100,2)
    survstatus_ICCC_2020 <- left_join(emptydf_surv, survstatus_ICCC_2020, by = c("survstatus_group"))
    survstatus_ICCC_2020$year <- rep(2020)
    survstatus_ICCC_2020$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    survstatus_ICCC_2020 <- data.table(survstatus_ICCC_2020)
    survstatus_ICCC2020 <- rbind(survstatus_ICCC2020, survstatus_ICCC_2020)
    
    #2020 - FU in 1yr periods until 5yrs
    survstatus_ICCC_20202 <- survstatus_ICCC_2020_all %>% group_by(survstatus_group2) %>% dplyr::count()
    survstatus_ICCC_20202$n <- (survstatus_ICCC_20202$n)*adj_fac
    survstatus_ICCC_20202$percent <- round((survstatus_ICCC_20202$n / sum(survstatus_ICCC_20202$n))*100,2)
    survstatus_ICCC_20202 <- left_join(emptydf_surv2, survstatus_ICCC_20202, by = c("survstatus_group2"))
    survstatus_ICCC_20202$year <- rep(2020)
    survstatus_ICCC_20202$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    survstatus_ICCC_20202 <- data.table(survstatus_ICCC_20202)
    survstatus_ICCC20202 <- rbind(survstatus_ICCC20202, survstatus_ICCC_20202)
    
    
    ###2040
    survstatus_ICCC_2040 <- survstatus[survstatus2$ICCCgroup == tomodel1$ICCCgroup[i] & survstatus2$dxyear <= 2040,]
    survstatus_ICCC_2040 <- data.table(survstatus_ICCC_2040)
    survstatus_ICCC_2040 <- survstatus_ICCC_2040 %>% select("2040")
    colnames(survstatus_ICCC_2040)[1] <- "survstatus"
    survstatus_ICCC_2040 <- na.omit(survstatus_ICCC_2040)
    survstatus_ICCC_2040_all <- survstatus_ICCC_2040 %>% mutate(survstatus_group = case_when(survstatus < 5 ~ "<5",
                                                                                             survstatus >= 5 & survstatus <= 9 ~ "5-9",
                                                                                             survstatus >= 10 & survstatus <= 14 ~ "10-14",
                                                                                             survstatus >= 15 & survstatus <= 19 ~ "15-19",
                                                                                             survstatus >= 20 ~ "20+"),
                                                                survstatus_group2 = case_when (survstatus <= 1 ~ "0-1 years",
                                                                                               survstatus <= 2 ~ "2 years",
                                                                                               survstatus <= 3 ~ "3 years",
                                                                                               survstatus <= 4 ~ "4 years",
                                                                                               survstatus >= 5 ~ "5+ years"))
    #2040 - #FU in 5 year periods
    survstatus_ICCC_2040 <- survstatus_ICCC_2040_all %>% group_by(survstatus_group) %>% dplyr::count()
    survstatus_ICCC_2040$n <- (survstatus_ICCC_2040$n)*adj_fac
    survstatus_ICCC_2040$percent <- round((survstatus_ICCC_2040$n / sum(survstatus_ICCC_2040$n))*100,2)
    survstatus_ICCC_2040 <- left_join(emptydf_surv, survstatus_ICCC_2040, by = c("survstatus_group"))
    survstatus_ICCC_2040$year <- rep(2040)
    survstatus_ICCC_2040$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    survstatus_ICCC_2040 <- data.table(survstatus_ICCC_2040)
    survstatus_ICCC2040 <- rbind(survstatus_ICCC2040, survstatus_ICCC_2040)
    
    #2020 - FU in 1yr periods until 5yrs
    survstatus_ICCC_20402 <- survstatus_ICCC_2040_all %>% group_by(survstatus_group2) %>% dplyr::count()
    survstatus_ICCC_20402$n <- (survstatus_ICCC_20402$n)*adj_fac
    survstatus_ICCC_20402$percent <- round((survstatus_ICCC_20402$n / sum(survstatus_ICCC_20402$n))*100,2)
    survstatus_ICCC_20402 <- left_join(emptydf_surv2, survstatus_ICCC_20402, by = c("survstatus_group2"))
    survstatus_ICCC_20402$year <- rep(2040)
    survstatus_ICCC_20402$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    survstatus_ICCC_20402 <- data.table(survstatus_ICCC_20402)
    survstatus_ICCC20402 <- rbind(survstatus_ICCC20402, survstatus_ICCC_20402)
  }
  
  #FU in 5 year periods
  survstatus_ICCC <- rbind(survstatus_ICCC2020, survstatus_ICCC2040)
  survstatus_ICCC <- filter(survstatus_ICCC, rowSums(is.na(survstatus_ICCC)) != ncol(survstatus_ICCC))
  survstatus_ICCC[is.na(survstatus_ICCC)] <- 0
  survstatus_ICCC_m <- survstatus_ICCC$n
  
  #FU in 1yr periods until 5yrs
  survstatus2_ICCC <- rbind(survstatus_ICCC20202, survstatus_ICCC20402)
  survstatus2_ICCC <- filter(survstatus2_ICCC, rowSums(is.na(survstatus2_ICCC)) != ncol(survstatus2_ICCC))
  survstatus2_ICCC[is.na(survstatus2_ICCC)] <- 0
  survstatus2_ICCC_m <- survstatus2_ICCC$n
  
  #Survival curve by decade of diagnosis
  
  cancer_id <- which(is.na(outcomes$df_X$dxyear )==F)
  cancer_df_X <- outcomes$df_X[cancer_id,]

  cancer_df_X$dead = ifelse(cancer_df_X$year_death < 9999, 1, 0) #use year of death to determine status instead of m_M 
  cancer_df_X$cre = ifelse(cancer_df_X$year_cre < 9999, 1, 0) #use year of cre to determine event status instead of m_M
  
  cancer_df_X$t_cancer = cancer_df_X$v_dx0 + cancer_df_X$v_cre0 #use time in cancer/CRE states to determine FU time instead of m_M 

  
  cancer_df_X$decade <- case_when(cancer_df_X$dxyear>=1970 & cancer_df_X$dxyear<1980 ~"70s",
                                  cancer_df_X$dxyear>=1980 & cancer_df_X$dxyear<1990 ~"80s",
                                  cancer_df_X$dxyear>=1990 & cancer_df_X$dxyear<2000 ~"90s",
                                  cancer_df_X$dxyear>=2000 & cancer_df_X$dxyear<2010 ~"00s",
                                  cancer_df_X$dxyear>=2010 & cancer_df_X$dxyear<2020 ~"10s",
                                  cancer_df_X$dxyear>=2020 & cancer_df_X$dxyear<2030 ~"20s",
                                  cancer_df_X$dxyear>=2030 & cancer_df_X$dxyear<=2040 ~"30s")
  
  cancer_df_X$period <- case_when(cancer_df_X$dxyear>=1970 & cancer_df_X$dxyear<1975 ~"1970-74",
                                        cancer_df_X$dxyear>=1975 & cancer_df_X$dxyear<1980 ~"1975-79",
                                        cancer_df_X$dxyear>=1980 & cancer_df_X$dxyear<1985 ~"1980-84",
                                        cancer_df_X$dxyear>=1985 & cancer_df_X$dxyear<1990 ~"1985-89",
                                        cancer_df_X$dxyear>=1990 & cancer_df_X$dxyear<1995 ~"1990-94",
                                        cancer_df_X$dxyear>=1995 & cancer_df_X$dxyear<2000 ~"1995-99",
                                        cancer_df_X$dxyear>=2000 & cancer_df_X$dxyear<2005 ~"2000-04",
                                        cancer_df_X$dxyear>=2005 & cancer_df_X$dxyear<2010 ~"2005-09",
                                        cancer_df_X$dxyear>=2010 & cancer_df_X$dxyear<2015 ~"2010-14",
                                        cancer_df_X$dxyear>=2015 & cancer_df_X$dxyear<2020 ~"2015-19",
                                        cancer_df_X$dxyear>=2020 & cancer_df_X$dxyear<2025 ~"2020-24",
                                        cancer_df_X$dxyear>=2025 & cancer_df_X$dxyear<2030 ~"2025-29",
                                        cancer_df_X$dxyear>=2030 & cancer_df_X$dxyear<2035 ~"2030-34",
                                        cancer_df_X$dxyear>=2035 & cancer_df_X$dxyear<2040 ~"2035-39")
  
  cancer_df_X$decade <- factor(cancer_df_X$decade, levels = c("70s", "80s", "90s", "00s", "10s", "20s", "30s"),
                               labels = c("1970-79", "1980-89", "1990-99", "2000-09", "2010-19", "2020-29", "2030-39"))
  
  cancer_df_X$period <- factor(cancer_df_X$period, levels = c("1970-74", "1975-79", "1980-84", "1985-89", "1990-94", "1995-99", "2000-04", "2005-09", "2010-14", "2015-19", "2020-24", "2025-29", "2030-34", "2035-39"),
                               labels = c("1970-74", "1975-79", "1980-84", "1985-89", "1990-94", "1995-99", "2000-04", "2005-09", "2010-14", "2015-19", "2020-24", "2025-29", "2030-34", "2035-39"))
  
  cancer_df_X$t_cancer[cancer_df_X$t_cancer==0] = 0.5
  
  #Overall survival - all cancers combined
  m1_decade <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X)
  sum_decade <- summary(m1_decade, times = 0:n_t)
  survtime <- data.frame(time = sum_decade$time, strata = sum_decade$strata, surv = sum_decade$surv)
  survtime <- left_join(template_surv, survtime, by = c("time", "strata"))
  surv_kk <- survtime$surv
  
  #Overall survival - ALL
  m1_decade_ALL <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "ALL"))
  sum_decade_ALL <-summary(m1_decade_ALL, times = 0:n_t)
  survtime_ALL <- data.frame(time = sum_decade_ALL$time, strata = sum_decade_ALL$strata, surv = sum_decade_ALL$surv)
  survtime_ALL <- left_join(template_surv, survtime_ALL, by = c("time", "strata"))
  surv_kk_ALL <- survtime_ALL$surv
  
  #Overall survival - AML and other leukemias
  m1_decade_AML <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "AML and other leukemias"))
  sum_decade_AML <- summary(m1_decade_AML, times = 0:n_t)
  survtime_AML <- data.frame(time = sum_decade_AML$time, strata = sum_decade_AML$strata, surv = sum_decade_AML$surv)
  survtime_AML <- left_join(template_surv, survtime_AML, by = c("time", "strata"))
  surv_kk_AML <- survtime_AML$surv
  
  #Overall survival - Non-Hodgkin lymphomas and other lymphomas
  m1_decade_NHL <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Non-Hodgkin lymphomas and other lymphomas"))
  sum_decade_NHL <-summary(m1_decade_NHL, times = 0:n_t)
  survtime_NHL <- data.frame(time = sum_decade_NHL$time, strata = sum_decade_NHL$strata, surv = sum_decade_NHL$surv)
  survtime_NHL <- left_join(template_surv, survtime_NHL, by = c("time", "strata"))
  surv_kk_NHL <- survtime_NHL$surv
  
  #Overall survival - Hodgkin lymphomas
  m1_decade_HL <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Hodgkin lymphomas"))
  sum_decade_HL <-summary(m1_decade_HL, times = 0:n_t)
  survtime_HL <- data.frame(time = sum_decade_HL$time, strata = sum_decade_HL$strata, surv = sum_decade_HL$surv)
  survtime_HL <- left_join(template_surv, survtime_HL, by = c("time", "strata"))
  surv_kk_HL <- survtime_HL$surv
  
  #Overall survival - Astrocytoma
  m1_decade_AST <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Astrocytoma"))
  sum_decade_AST <-summary(m1_decade_AST, times = 0:n_t)
  survtime_AST <- data.frame(time = sum_decade_AST$time, strata = sum_decade_AST$strata, surv = sum_decade_AST$surv)
  survtime_AST <- left_join(template_surv, survtime_AST, by = c("time", "strata"))
  surv_kk_AST <- survtime_AST$surv
  
  #Overall survival - Other CNS neoplasms
  m1_decade_OCNS <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Other CNS neoplasms"))
  sum_decade_OCNS <-summary(m1_decade_OCNS, times = 0:n_t)
  survtime_OCNS <- data.frame(time = sum_decade_OCNS$time, strata = sum_decade_OCNS$strata, surv = sum_decade_OCNS$surv)
  survtime_OCNS <- left_join(template_surv, survtime_OCNS, by = c("time", "strata"))
  surv_kk_OCNS <- survtime_OCNS$surv
  
  #Overall survival - Neuroblastoma
  m1_decade_NEU <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Neuroblastoma"))
  sum_decade_NEU <-summary(m1_decade_NEU, times = 0:n_t)
  survtime_NEU <- data.frame(time = sum_decade_NEU$time, strata = sum_decade_NEU$strata, surv = sum_decade_NEU$surv)
  survtime_NEU <- left_join(template_surv, survtime_NEU, by = c("time", "strata"))
  surv_kk_NEU <- survtime_NEU$surv
  
  #Overall survival - Retinoblastoma
  m1_decade_RET <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Retinoblastoma"))
  sum_decade_RET <-summary(m1_decade_RET, times = 0:n_t)
  survtime_RET <- data.frame(time = sum_decade_RET$time, strata = sum_decade_RET$strata, surv = sum_decade_RET$surv)
  survtime_RET <- left_join(template_surv, survtime_RET, by = c("time", "strata"))
  surv_kk_RET <- survtime_RET$surv
  
  #Overall survival - Renal tumours
  m1_decade_REN <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Renal tumours"))
  sum_decade_REN <-summary(m1_decade_REN, times = 0:n_t)
  survtime_REN <- data.frame(time = sum_decade_REN$time, strata = sum_decade_REN$strata, surv = sum_decade_REN$surv)
  survtime_REN <- left_join(template_surv, survtime_REN, by = c("time", "strata"))
  surv_kk_REN <- survtime_REN$surv
  
  #Overall survival - Hepatic tumours
  m1_decade_HEP <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Hepatic tumours"))
  sum_decade_HEP <-summary(m1_decade_HEP, times = 0:n_t)
  survtime_HEP <- data.frame(time = sum_decade_HEP$time, strata = sum_decade_HEP$strata, surv = sum_decade_HEP$surv)
  survtime_HEP <- left_join(template_surv, survtime_HEP, by = c("time", "strata"))
  surv_kk_HEP <- survtime_HEP$surv
  
  #Overall survival - Bone tumours
  m1_decade_BONE <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Bone tumours"))
  sum_decade_BONE <-summary(m1_decade_BONE, times = 0:n_t)
  survtime_BONE <- data.frame(time = sum_decade_BONE$time, strata = sum_decade_BONE$strata, surv = sum_decade_BONE$surv)
  survtime_BONE <- left_join(template_surv, survtime_BONE, by = c("time", "strata"))
  surv_kk_BONE <- survtime_BONE$surv
  
  #Overall survival - Soft tissue sarcomas
  m1_decade_SAR <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Soft tissue sarcomas"))
  sum_decade_SAR <-summary(m1_decade_SAR, times = 0:n_t)
  survtime_SAR <- data.frame(time = sum_decade_SAR$time, strata = sum_decade_SAR$strata, surv = sum_decade_SAR$surv)
  survtime_SAR <- left_join(template_surv, survtime_SAR, by = c("time", "strata"))
  surv_kk_SAR <- survtime_SAR$surv
  
  #Overall survival - Germ cell tumours
  m1_decade_GERM <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Germ cell tumours"))
  sum_decade_GERM <-summary(m1_decade_GERM, times = 0:n_t)
  survtime_GERM <- data.frame(time = sum_decade_GERM$time, strata = sum_decade_GERM$strata, surv = sum_decade_GERM$surv)
  survtime_GERM <- left_join(template_surv, survtime_GERM, by = c("time", "strata"))
  surv_kk_GERM <- survtime_GERM$surv
  
  #Overall survival - Other epithelial and unspecified neoplasms
  m1_decade_OTHER <- survfit(Surv(t_cancer, dead)~ decade, data = cancer_df_X%>%filter(ICCCgroup == "Other epithelial and unspecified neoplasms"))
  sum_decade_OTHER <-summary(m1_decade_OTHER, times = 0:n_t)
  survtime_OTHER <- data.frame(time = sum_decade_OTHER$time, strata = sum_decade_OTHER$strata, surv = sum_decade_OTHER$surv)
  survtime_OTHER <- left_join(template_surv, survtime_OTHER, by = c("time", "strata"))
  surv_kk_OTHER <- survtime_OTHER$surv
  
  #Overall survival - all cancers combined and 5-year period of diagnosis
  m1_decade <- survfit(Surv(t_cancer, dead)~ period, data = cancer_df_X)
  sum_decade <- summary(m1_decade, times = 0:n_t)
  survtime <- data.frame(time = sum_decade$time, strata = sum_decade$strata, surv = sum_decade$surv)
  survtime <- left_join(template_surv5, survtime, by = c("time", "strata"))
  surv_kk_5 <- survtime$surv
  
  ###### Loop approach - #Overall survival - by cancer type and 5-year period of diagnosis
  tomodel1  <- cancer_df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  survprobs_cancertype <- vector(mode  = "list") #create the list by cancer type
  
  for(i in 1:dim(tomodel1)[1]){ 
    
    cancer_df_X_subset  <- filter(cancer_df_X, ICCCgroup == tomodel1$ICCCgroup[i])
    
    m1_decade_type <- survfit(Surv(t_cancer, dead)~ period, data = cancer_df_X_subset)
    sum_decade_type <-summary(m1_decade_type, times = 0:n_t)
    survtime_type <- data.frame(time = sum_decade_type$time, strata = sum_decade_type$strata, surv = sum_decade_type$surv)
    survtime_type <- left_join(template_surv5, survtime_type, by = c("time", "strata"))
    survtime_type$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    
    survprobs_cancertype[[i]] <- survtime_type
    names(survprobs_cancertype)[i] <- paste0(tomodel1$ICCCgroup[i])
  }
  survprobs_cancertype_all <- ldply(survprobs_cancertype, data.frame)
  surv_kk_ICCC_5 <- survprobs_cancertype_all$surv
  
  ##################Event-free survival
 
  cancer_df_X$pfs_status <- ifelse(cancer_df_X$dead == 1 | cancer_df_X$cre == 1, 1, 0)
  
  #for those who only spend 1 year in CRE state (v_cre0 == 1), year_cre is always 2040 (incorrect) - line below fixes
  cancer_df_X$year_cre2 <- ifelse(cancer_df_X$v_cre0 == 1 & cancer_df_X$year_cre == 2040, (cancer_df_X$dxyear + cancer_df_X$v_dx0), cancer_df_X$year_cre)
  
  cancer_df_X <- cancer_df_X %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(year_pfs_event = min(year_cre2, year_death, na.rm = TRUE)) %>% #find the first occurring event of cre or death
  ungroup()
  
  cancer_df_X$time_pfs_event <- ifelse(cancer_df_X$pfs_status == 1, (cancer_df_X$year_pfs_event - cancer_df_X$dxyear), cancer_df_X$t_cancer) #for people with an PFS event, calculate years between event/dx. for those without event, use time in cancer state (Event free)

  # cancer_df_X <- cancer_df_X %>%
  # rowwise() %>%
  # mutate(year_event = min(year_cre, year_death, year_emm, na.rm = TRUE),
  #        event = ifelse(dead == 1 | cre == 1 | !is.na(year_emm), 1, 0)) #find the first occurring event
  # 
  # cancer_df_X <- cancer_df_X %>%
  #   rowwise() %>%
  #   mutate(time_event = ifelse(event == 1, (year_event - dxyear), t_cancer))
  
  #Event-free survival - all cancers combined
  m1_decade <- survfit(Surv(time_pfs_event, pfs_status)~ decade, data = cancer_df_X)
  sum_decade <- summary(m1_decade, times = 0:n_t)
  survtime <- data.frame(time = sum_decade$time, strata = sum_decade$strata, surv = sum_decade$surv)
  survtime <- left_join(template_surv, survtime, by = c("time", "strata"))
  surv_pfs_kk <- survtime$surv
  
  ###### Loop approach - #Event-free survival - by cancer type and decade of diagnosis
  tomodel1  <- cancer_df_X %>% dplyr::count(ICCCgroup) %>% filter(is.na(ICCCgroup) !=1)
  survprobs_cancertype <- vector(mode  = "list") #create the list by cancer type
  
  for(i in 1:dim(tomodel1)[1]){ 
    
    cancer_df_X_subset  <- filter(cancer_df_X, ICCCgroup == tomodel1$ICCCgroup[i])
    
    m1_decade <- survfit(Surv(time_pfs_event, pfs_status)~ decade, data = cancer_df_X_subset)
    sum_decade <- summary(m1_decade, times = 0:n_t)
    survtime <- data.frame(time = sum_decade$time, strata = sum_decade$strata, surv = sum_decade$surv)
    survtime <- left_join(template_surv, survtime, by = c("time", "strata"))
    survtime$ICCCgroup <- rep(tomodel1$ICCCgroup[i])
    
    survprobs_cancertype[[i]] <- survtime
  
    names(survprobs_cancertype)[i] <- paste0(tomodel1$ICCCgroup[i])
  }
  survprobs_cancertype_all <- ldply(survprobs_cancertype, data.frame)
  surv_pfs_kk_ICCC <- survprobs_cancertype_all$surv
  
  ##### #
  
  
  return(list(
    prev_sim                 = prev_sim,
    prev_rate                = prev_rate,
    prevpop_PSA              = prevpop_PSA,
    prev_imm                 = prev_imm,
    prev_bornON              = prev_bornON,
    prev_imm_dxON            = prev_imm_dxON,
    prev_sim_true_prev       = prev_sim_true_prev,
    prev_sim_true_prev_70    = prev_sim_true_prev_70,
    prev_sim_true_prev_ICCC  = prev_sim_true_prev_ICCC,
    prev_sim_true_prev_70_ICCC = prev_sim_true_prev_70_ICCC,
    prev_sim_ICCC            = prev_sim_ICCC,
    prev_rate_ICCC           = prev_rate_ICCC,
    prev_sim_char            = prev_sim_char,
    prev_sim_ICCC_char       = prev_sim_ICCC_char,
    childpop_PSA             = childpop_PSA,
    dxyear_sim_imm           = dxyear_sim_imm,
    dxyeargroup_sim_age      = dxyeargroup_sim_age,
    dxyear_sim_imm_rate      = dxyear_sim_imm_rate,
    dxyear_sim_imm5_df       = dxyear_sim_imm5_df,
    dxyear_sim_imm5_rate_df  = dxyear_sim_imm5_rate_df,
    dxyeargroup_sim          = dxyeargroup_sim,
    dxyeargroup_iccc_df      = dxyeargroup_iccc_df,
    dxyeargroup_iccc_rate_df = dxyeargroup_iccc_rate_df,
    attainedage_year_df_m    = attainedage_year_df_m,
    attainedage_aspr_df_m    = attainedage_aspr_df_m,
    attainedage_ICCC_m       = attainedage_ICCC_m,
    attainedage_ICCC_aspr_m  = attainedage_ICCC_aspr_m,
    survstatus_df            = survstatus_df,
    survstatus_df2           = survstatus_df2,
    survstatus_ICCC_m        = survstatus_ICCC_m,
    survstatus2_ICCC_m       = survstatus2_ICCC_m,
    surv_kk                  = surv_kk,
    surv_kk_ALL              = surv_kk_ALL,
    surv_kk_AML              = surv_kk_AML,
    surv_kk_NHL              = surv_kk_NHL,
    surv_kk_HL               = surv_kk_HL,
    surv_kk_AST              = surv_kk_AST,
    surv_kk_OCNS             = surv_kk_OCNS,
    surv_kk_NEU              = surv_kk_NEU,
    surv_kk_RET              = surv_kk_RET,
    surv_kk_REN              = surv_kk_REN,
    surv_kk_HEP              = surv_kk_HEP,
    surv_kk_BONE             = surv_kk_BONE,
    surv_kk_SAR              = surv_kk_SAR,
    surv_kk_GERM             = surv_kk_GERM,
    surv_kk_OTHER            = surv_kk_OTHER,
    surv_kk_ICCC_5           = surv_kk_ICCC_5,
    surv_kk_5                = surv_kk_5,
    surv_pfs_kk              = surv_pfs_kk,
    surv_pfs_kk_ICCC         = surv_pfs_kk_ICCC)
  )

} #end of the simulation loop

## Store results from each iteration ---------------------------------------------------------------------------

for (k in 1:n_sim) {
  sim_results_iter             <- sim_results[[k]]
  prev_sim[,k]                 <- sim_results_iter$prev_sim
  prev_rate[,k]                <- sim_results_iter$prev_rate
  prevpop_PSA[,k]              <- sim_results_iter$prevpop_PSA
  prev_imm[,k]                 <- sim_results_iter$prev_imm
  prev_bornON[,k]              <- sim_results_iter$prev_bornON
  prev_imm_dxON[,k]            <- sim_results_iter$prev_imm_dxON
  prev_sim_true_prev[,k]       <- sim_results_iter$prev_sim_true_prev
  prev_sim_true_prev_70[,k]    <- sim_results_iter$prev_sim_true_prev_70
  prev_sim_true_prev_ICCC[,k]  <- sim_results_iter$prev_sim_true_prev_ICCC
  prev_sim_true_prev_70_ICCC[,k] <- sim_results_iter$prev_sim_true_prev_70_ICCC
  prev_sim_ICCC[,k]          <- sim_results_iter$prev_sim_ICCC
  prev_rate_ICCC[,k]         <- sim_results_iter$prev_rate_ICCC
  prev_sim_char[,k]           <- sim_results_iter$prev_sim_char
  prev_sim_ICCC_char[,k]      <- sim_results_iter$prev_sim_ICCC_char
  childpop_PSA[,k]             <- sim_results_iter$childpop_PSA
  dxyear_sim_imm[,k]           <- sim_results_iter$dxyear_sim_imm
  dxyeargroup_sim_age[,k]      <- sim_results_iter$dxyeargroup_sim_age
  dxyear_sim_imm_rate[,k]      <- sim_results_iter$dxyear_sim_imm_rate
  dxyear_sim_imm5_df[,k]       <- sim_results_iter$dxyear_sim_imm5_df
  dxyear_sim_imm5_rate_df[,k]  <- sim_results_iter$dxyear_sim_imm5_rate_df
  dxyeargroup_sim[,k]          <- sim_results_iter$dxyeargroup_sim
  dxyeargroup_iccc_df[,k]      <- sim_results_iter$dxyeargroup_iccc_df
  dxyeargroup_iccc_rate_df[,k] <- sim_results_iter$dxyeargroup_iccc_rate_df
  attainedage_year_df_m[,k]    <- sim_results_iter$attainedage_year_df_m
  attainedage_aspr_df_m[,k]    <- sim_results_iter$attainedage_aspr_df_m
  attainedage_ICCC_m[,k]       <- sim_results_iter$attainedage_ICCC_m
  attainedage_ICCC_aspr_m[,k]  <- sim_results_iter$attainedage_ICCC_aspr_m
  survstatus_df[,k]            <- sim_results_iter$survstatus_df
  survstatus_df2[,k]           <- sim_results_iter$survstatus_df2
  survstatus_ICCC_m[,k]        <- sim_results_iter$survstatus_ICCC_m
  survstatus2_ICCC_m[,k]       <- sim_results_iter$survstatus2_ICCC_m
  surv_kk[,k]                  <- sim_results_iter$surv_kk
  surv_kk_ALL[,k]              <- sim_results_iter$surv_kk_ALL
  surv_kk_AML[,k]              <- sim_results_iter$surv_kk_AML
  surv_kk_NHL[,k]              <- sim_results_iter$surv_kk_NHL
  surv_kk_HL[,k]               <- sim_results_iter$surv_kk_HL
  surv_kk_AST[,k]              <- sim_results_iter$surv_kk_AST
  surv_kk_OCNS[,k]             <- sim_results_iter$surv_kk_OCNS
  surv_kk_NEU[,k]              <- sim_results_iter$surv_kk_NEU
  surv_kk_RET[,k]              <- sim_results_iter$surv_kk_RET
  surv_kk_REN[,k]              <- sim_results_iter$surv_kk_REN
  surv_kk_HEP[,k]              <- sim_results_iter$surv_kk_HEP
  surv_kk_BONE[,k]             <- sim_results_iter$surv_kk_BONE
  surv_kk_SAR[,k]              <- sim_results_iter$surv_kk_SAR
  surv_kk_GERM[,k]             <- sim_results_iter$surv_kk_GERM
  surv_kk_OTHER[,k]            <- sim_results_iter$surv_kk_OTHER
  surv_kk_ICCC_5[,k]           <- sim_results_iter$surv_kk_ICCC_5
  surv_kk_5[,k]                <- sim_results_iter$surv_kk_5
  surv_pfs_kk[,k]              <- sim_results_iter$surv_pfs_kk
  surv_pfs_kk_ICCC[,k]         <- sim_results_iter$surv_pfs_kk_ICCC
}

## Save model output from each iteration -----------------------------------------------------------------------------

#Save locally
print(paste("Ready to save batch =", b)) #print which batch has finished running

write.csv(prev_sim, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/prev_sim_kk.csv"), row.names = F)
write.csv(prev_sim_char, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/prev_sim_char_kk.csv"), row.names = F)
write.csv(prev_sim_ICCC_char, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/prev_sim_ICCC_char_kk.csv"), row.names = F)

write.csv(attainedage_aspr_df_m, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/attainedage_aspr_df_m_kk.csv"), row.names = F)
write.csv(attainedage_ICCC_aspr_m, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/attainedage_ICCC_aspr_m_kk.csv"), row.names = F)
write.csv(survstatus_df, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/survstatus_df_kk.csv"), row.names = F)
write.csv(survstatus_ICCC_m, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/survstatus_ICCC_m_kk.csv"), row.names = F)

write.csv(surv_kk, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/surv_kk.csv"), row.names = F)
write.csv(surv_kk_ICCC_5, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/surv_kk_ICCC_5.csv"), row.names = F)
write.csv(surv_kk_5, paste0("~/GitHub/POSIM-prev/analysis/SA_surv_output/batch",b,"rawdataoutput/surv_kk_5.csv"), row.names = F)

} #end batch loop

