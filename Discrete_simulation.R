
library(tidyverse)
library(here)
library(parallel)
library(rUCLs)
library(EnvStats)

h <- here::here

#########################################################
no_cores <- detectCores()
clust <- makeCluster(no_cores) 
clusterCall(clust, function() library(rUCLs))
options(stringsAsFactors = FALSE)
#########################################################
# User defined assumptions: 
## Static ##

# Concentration (parent)
mean_high <- 600
sd_high <- 1800
mean_low <- 100
sd_low <- 300
sd_high/mean_high == sd_low/mean_low #check CVs match
cv <- sd_high/mean_high
cv <- sprintf("%.1f", cv)

# Coverage
total_area <- 100 #acres
du_size <- 1 #acres
sample_size <- 15
total_du_count <- round(total_area/du_size, 0)   

## Ranges ##

# Probability that a DU is in "high" category
prob_high_min <- 0.05 
prob_high_max <- 0.25

# Proportion of site sampled in pilot study
prop_sampled_min <-0.1
prop_sampled_max <-0.9

# Action level
action_level_min <- 0.1 * mean_low
action_level_max <- 10 * mean_high

## Monte Carlo ##
iterations <- 10000 

#########################################################
# Assign each DU to either the “Low” or “High” category based on random sample of proportion "high" 
prob_high <- rep(runif(n = iterations,
                       min = prob_high_min,
                       max = prob_high_max),
                 each = total_du_count) #random sample from uniform distribution for "high" probability

high_designation <- rbinom(n = total_du_count*iterations, 
                           size = 1, 
                           prob = prob_high) #based on probability of "high", sample high/low

site_dat <- tibble(du = rep(1:total_du_count, iterations),
                   iteration = rep(1:iterations,each=total_du_count),
                   prob_high = prob_high,
                   concentration_class = ifelse(high_designation==1, "High", "Low"),
                   parent_mean = ifelse(high_designation==1, mean_high, mean_low),
                   parent_sd = ifelse(high_designation==1, sd_high, sd_low)) #compile site data

#########################################################
# Assign Action level
site_dat$action_level <- rep(runif(n = iterations,
                                   min = action_level_min,
                                   max = action_level_max),
                             each = total_du_count) #random sample from uniform distribution for "high" probability

#########################################################
# Sample for the "true" mean/sd (outer loop)

# Assume underlying site data follow log-normal distribution
meanlog <- log(site_dat$parent_mean^2 / sqrt(site_dat$parent_sd^2 + site_dat$parent_mean^2)) #arithmetic mean to mean log
sdlog <- sqrt(log(1 + (site_dat$parent_sd^2 / site_dat$parent_mean^2))) #sd to sd log

site_true <- t(mapply(rlnorm,
                      n = sample_size,
                      meanlog = meanlog,
                      sdlog = sdlog))  #generate n samples per DU for each iteration

dimnames(site_true) <- list(
  d1=rep(paste0("DU",1:total_du_count, ".i", rep(1:iterations, each=total_du_count))),
  d2=paste0("n",1:sample_size)) #name array

site_true <- as.data.frame.table(site_true) #convert to data frame and clean up names
names(site_true) <- c("du", "sample", "concentration")

site_dat <- site_true %>%
  mutate(iteration =  as.numeric(gsub("i", "", gsub(".*\\.","", du))),
         du = as.numeric(gsub("DU", "", gsub("\\..*","", du)))) %>%
  group_by(iteration, du) %>%
  summarize(true_mean = mean(concentration),
            true_sd = sd(concentration)) %>%
  left_join(site_dat,.)

#########################################################
# Randomly choose subset of DUs to be the Pilot Study subjects.
prop_sampled <- runif(iterations,
                      min = prop_sampled_min,
                      max = prop_sampled_max) #random sample from uniform  for prop. of site sampled

pilot_du_count <- round((total_area*prop_sampled)/du_size,0) #calculate number of DUs sampled

sampled_du <- do.call("rbind", 
                      (lapply(1:length(pilot_du_count),
                              function(i)
                                tibble(du = 1:total_du_count, 
                                       iteration = i, 
                                       pilot_du_count = pilot_du_count[i],
                                       pilot = du %in% sample(x=1:total_du_count, 
                                                              size=pilot_du_count, 
                                                              replace=FALSE))))) #assign DUs to sampled/unsampled

site_dat <- left_join(site_dat, sampled_du)

#########################################################
# Sample for the "pilot" mean/sd (inner loop)

#Assume true site data follow log-normal distribution
meanlog2 <- log(site_dat$true_mean^2 / sqrt(site_dat$true_sd^2 + site_dat$true_mean^2)) #arithmetic mean to mean log
sdlog2 <- sqrt(log(1 + (site_dat$true_sd^2 / site_dat$true_mean^2))) #sd to sd log

site_pilot<- t(mapply(rlnorm,
                      n=sample_size,
                      meanlog = meanlog2,
                      sdlog = sdlog2))  #generate n samples per DU for each iteration

dimnames(site_pilot) <- list(
  d1=rep(paste0("DU",1:total_du_count, ".i", rep(1:iterations, each=total_du_count))),
  d2=paste0("n",1:sample_size)) #name array

site_pilot <- as.data.frame.table(site_pilot) #convert to data frame and clean up names
names(site_pilot) <- c("du", "sample", "concentration")
site_pilot <- site_pilot %>%
  mutate(iteration =  as.numeric(gsub("i", "", gsub(".*\\.","", du))),
         du = as.numeric(gsub("DU", "", gsub("\\..*","", du)))) %>%
  left_join(., sampled_du)

#########################################################
# Summary statistics for each DU & iteration (uses Integral rUCL package)
g <-unique(site_pilot[site_pilot$pilot==TRUE,c("du", "iteration")]) #dataframe of unique DU and iteration combinations

clusterExport(clust, "g") #for parallel processing
clusterExport(clust, "site_pilot") #for parallel processing

# NEEDS PRIVATE rUCL package to run
gof <-parLapply(clust,
                1:nrow(g),
                function(i) 
                  proucl_gof(
                    site_pilot$concentration[site_pilot$du==g$du[i] &
                                               site_pilot$iteration==g$iteration[i]], 
                    det = TRUE,
                    quals= "", 
                    include_nd = TRUE)) #fits distribution to n samples in each DU for each iteration

clusterExport(clust, "gof") #for parallel processing

ucl_output <- do.call("rbind",
                      parLapply(clust,
                                1:nrow(g),
                                function(i)                  
                                  ucl(site_pilot$concentration[site_pilot$du==g$du[i] &
                                                                 site_pilot$iteration==g$iteration[i]],
                                      det=TRUE,
                                      gof[[i]],
                                      boot=FALSE, 
                                      max.sub=FALSE,
                                      method="Select"))) #sum stats based on GOF

stopCluster(clust)

ucl_output <- bind_cols(g,ucl_output) %>%
  select(du, iteration, epc)

site_dat <- left_join(site_dat, ucl_output) %>%
  mutate(true_mean_exceeds = ifelse(true_mean > action_level, 1, 0),
         pilot_ucl_exceeds = ifelse(epc > action_level, 1, 0)) #assign exceedances


#########################################################
# Decision logic
site_dat_sum <- site_dat %>%
  group_by(iteration, prob_high) %>%
  summarize(du_count=n(),
            count_high = length(concentration_class[concentration_class=="High"]),
            parent_mean_high =  mean_high,
            parent_sd_high = sd_high,
            parent_mean_low =  mean_low,
            parent_sd_low = sd_low,
            pilot_du_count = unique(pilot_du_count),
            pilot_ucl =  mean(epc, na.rm = T),
            action_level = unique(action_level),
            pilot_ratio = mean(epc, na.rm = T)/action_level,
            count_pilot_ucl_exceeds = sum(pilot_ucl_exceeds, na.rm = T),
            count_true_mean_exceeds = sum(true_mean_exceeds),
            prop_pilot_ucl_exceeds = sum(pilot_ucl_exceeds, na.rm = T)/pilot_du_count,
            prop_true_mean_exceeds = sum(true_mean_exceeds)/n()) %>%
  mutate(a= ifelse(count_pilot_ucl_exceeds>0, 0, 1), #0 = non-compliance (do not walk away), 1 = compliance (walk away)
         b= ifelse(count_true_mean_exceeds>0, 0, 1)) %>% #0 = non-compliance (true), 1 = compliance (true)
  mutate(match = ifelse(a==0 & b==0, 1, 
                        ifelse(a<b, 2,
                               ifelse(a>b, 3, 
                                      4))),
         match_desc = ifelse(a==0 & b==0, "Not in Compliance - Correct", 
                             ifelse(a<b, "Not in Compliance - Error",
                                    ifelse(a>b, "Compliance - Error", 
                                           "Compliance - Correct"))))

write_csv(site_dat_sum, h("output", paste0("SimulationDecision_CV", cv, ".csv")))

#########################################################
#Output summary for excel figures - by Ratio
dat <- read_csv(file = h("output", paste0("SimulationDecision_CV", cv, ".csv")))

dat <- dat %>%
  mutate(ratio_bin = .bincode(pilot_ratio, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, Inf), right = TRUE, include.lowest = TRUE))

bin_cat <- tibble(ratio_bin = 1:11,
                  ratio_bin_label = factor(c("R <= 0.1",
                                             paste(seq(0.1, 0.9, 0.1), "< R <=", seq(0.2, 1.0, 0.1)),
                                             "R > 1.0")))

dat <- left_join(dat, bin_cat)

dat_sum1 <- dat %>%
  select(ratio_bin_label, match)

dat_sum1 <-table(dat_sum1)
dat_sum1 <-dat_sum1[c(10, 1:9, 11),]

write.csv(dat_sum1, h("output", paste0("MatchSummary_byRatio_CV", cv, ".csv")), row.names = T)

#########################################################
#Output summary for excel figures - by % site
dat <- read_csv(file = h("output", paste0("SimulationDecision_CV", cv, ".csv")))

dat <- dat %>%
  mutate(pilot_du_perc = pilot_du_count/100) %>%
  mutate(du_bin = .bincode(pilot_du_perc, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), right = FALSE, include.lowest = TRUE))

bin_cat <- tibble(du_bin = 2:9,
                  du_bin_label = factor(c(paste(seq(0.1, 0.7, 0.1), "<= Pilot <", seq(0.2, 0.8, 0.1)),
                                          "0.8 <= Pilot <= 0.9"
                  )))

dat <- left_join(dat, bin_cat)

dat_sum1 <- dat %>%
  select(du_bin_label, match)

dat_sum1 <- table(dat_sum1)

write.csv(dat_sum1, h("output", paste0("MatchSummary_byDUCount_CV", cv, ".csv")), row.names = T)
