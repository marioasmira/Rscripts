# Set up ------------------------------------------------------------------
library(rethinking)
library(brms)
library(ggplot2)
library(cowplot)
setwd("~/Code/Rscripts")
set.seed(1234)


# Intro -------------------------------------------------------------------

# test two genotypes for egg volume
# we know:
# - females lay Poisson(18.5) eggs per day
# - the volume amounts are roughly 8~10 mm^3 x 10^-3
# - there may be a high correlation within mothers

# what we want to know:
# - do we need to measure individual mothers?
# - if so, how many in order to distinguish small differences?

# variables
egg_lambda <- 18.5  # measured
std_vol <- 8        # measured
std_sd <- 1      # measured

# questions:
# 1) assuming the same volume, how many mothers do we need to measure to know if there's a high within mother correlation?
# 2) assuming several levels of correlation within mothers, how accurate is a bunched vs individual measurement when trying to differentiate two volumes?

# what we will vary
sd_diffs <- c(0.5, 0, -0.5)
vol_diffs <- c(0.5, 1, 2)
n_mothers <- c(10, 30, 60, 100)

# here we will ignore the vol_diffs variable and focus on trying to find the maternal variance

# Question 1 --------------------------------------------------------------

# now we want to simulate our experiment
# we will have n mothers each in it's own petri dish laying eggs
# we will collect 5, 10 or 20 egg from each mother
n_eggs <- c(5, 10, 20)
# we will measure the volume of those eggs analyze the data in the following way vol ~ 1 + (1 | mother_id)
# so we will have 3 measurements to test:
# 1) how accurate is the intercept
# 2) how accurate is the maternal sd
# 3) how accurate is the population sd

# to simulate we create a data.frame and then fill it with all the combination of the relevant parameters:
sim_maternal_effect <- data.frame(
  sd_diff = NA,
  n_mothers = NA,
  n_eggs = NA,
  mother_ID = NA,
  measure = NA
)
sim_maternal_effect <- na.omit(sim_maternal_effect)

# simulate
# for each value of sd difference
for(i in 1:length(sd_diffs)){
  
  # for each number of eggs
  for(l in 1:length(n_eggs)){
    
    # for each number of mothers
    for(j in 1:length(n_mothers)){
      
      # for each mother on the respective number
      for(k in 1:n_mothers[j]){
        
        # each mother has it's own mean value of volume which is
        mean_mother_k <- std_vol + rnorm(1, mean = 0, sd = std_sd - sd_diffs[i])
        
        # we make a temporary table to hold each mother's data (this is equal to the one we built before)
        temp_eggs <- data.frame(
          sd_diff = NA,
          n_mothers = NA,
          n_eggs = NA,
          mother_ID = NA,
          measure = NA
        )
        
        # for each mother we collect n_eggs (5) eggs
        for(e in 1:n_eggs[l]){
          
          # each egg is equal to the volume of the mother plus the population level error
          egg_vol <- mean_mother_k + rnorm(1, mean = 0, sd = std_sd + sd_diffs[i])
          
          # now we fill the temporary table
          temp_eggs[e,1] <- sd_diffs[i]
          temp_eggs[e,2] <- n_mothers[j]
          temp_eggs[e,3] <- n_eggs[l]
          temp_eggs[e,4] <- k
          temp_eggs[e,5] <- egg_vol
        }
        
        # we add the temporary table to the end of the large data frame
        sim_maternal_effect <- rbind(sim_maternal_effect, temp_eggs)
      }
    }
  }
}

# next we set up all combinations of parameters
reps <- length(sd_diffs) * length(n_mothers) * length(n_eggs)
diff_values <- sort(rep(sd_diffs, reps / length(sd_diffs)))
n_mom_values <- rep(n_mothers, reps / length(n_mothers))
n_eggs_values <- rep(sort(rep(n_eggs, reps / (length(n_eggs) * length(sd_diffs)))), length(sd_diffs))
  
# and check if combinations are correctly set up
unique(paste(diff_values, n_mom_values, n_eggs_values, sep = ""))

# if so, build a list od models to store everything
models <- vector("list", length(diff_values))

# analyse each amount - takes a looong time
for(i in 1:length(models)){
  models[[i]] <- brm(
    measure ~ 1 + (1 | mother_ID),
    data = subset(sim_maternal_effect,
                  n_mothers == n_mom_values[i] &
                    sd_diff == diff_values[i] &
                    n_eggs == n_eggs_values[i]),
    
    cores = 8,
    chains = 8,
    iter = 5000,
    warmup = 1000
  )
  
  print(paste("Now finished loop:", as.character(i), sep = " "))
}

# Save an object to a file
saveRDS(models, file = "sim_models_example_1.rds")

# Restore the object
models <- readRDS(file = "sim_models_example_1.rds")


# making a data table for the output
sim_result <- data.frame(
  diff = NA,
  n_mom = NA,
  n_eggs = NA,
  sim_ICC = NA,
  sim_diff = NA,
  sim_int = NA
)
sim_result <- na.omit(sim_result)

# function to calculate intraclass correlation coefficient
# ICC = var(between moms)/(var(between moms)+var(within moms))
ICC <- function(alpha, epsilon) {
  sig_u <- alpha^2
  sig_u/(sig_u+(epsilon)^2)
}

for(i in 1:length(models)){
  sd <- posterior_samples(models[[i]], pars = "sd_mother_ID__Intercept")
  sigma <- posterior_samples(models[i], pars = "sigma")
  interc <- std_vol - posterior_samples(models[i], pars = "b_Intercept")
  diff <- diff_values[i]
  n_mom <- n_mom_values[i]
  n_egg <- n_eggs_values[i]

  temp_diff <- data.frame(
    diff,
    n_mom,
    n_egg,
    ICC(sd, sigma),
    (sigma - sd ),
    interc
  )
  
  colnames(temp_diff) <- c("diff", "n_mom", "n_eggs", "sim_ICC", "sim_diff", "sim_int")
  sim_result <- rbind(sim_result, temp_diff)
}


# plot differences between population sigma and mother sd
sim_result$comb <- paste(sim_result$diff, sim_result$n_eggs, sim_result$n_mom, sep = "")
treatments_summ <- unique(sim_result$comb)
summary <- data.frame(comb = treatments_summ,
                      diff = NA,
                      n_mom = NA,
                      n_eggs = NA,
                      mean = NA,
                      SEM = NA,
                      SD = NA,
                      UHPDI = NA,
                      LHPDI = NA,
                      highbound = NA,
                      lowbound = NA
)

for(i in 1:length(treatments_summ)){
  data_temp <- subset(sim_result, comb == summary$comb[i])
  summary$diff[i] <- ifelse(data_temp$diff[1] == 0.5, "Mom > Res", ifelse(data_temp$diff[1] == 0, "Mom = Res", "Mom < Res"))
  summary$n_mom[i] <- data_temp$n_mom[1]
  summary$n_eggs[i] <- data_temp$n_eggs[1]
  summary$mean[i] <- mean(data_temp$sim_diff, na.rm = TRUE)
  summary$SEM[i] <-sd(data_temp$sim_diff, na.rm = TRUE) /
    sqrt(length(data_temp$sim_diff))
  summary$SD[i] <- sd(data_temp$sim_diff, na.rm = TRUE)
  summary$UHPDI[i] <- HPDI(data_temp$sim_diff, prob = 0.95)[2]
  summary$LHPDI[i] <- HPDI(data_temp$sim_diff, prob = 0.95)[1]
  summary$highbound[i] <- summary$UHPDI[i] - summary$mean[i]
  summary$lowbound[i] <- summary$mean[i] - summary$LHPDI[i]
}

g1 <- ggplot(data = summary, aes(x = n_mom,
                                 y = mean,
                                 colour = as.factor(diff),
                                 shape = as.factor(n_eggs),
                                 linetype = as.factor(diff))) +
  geom_point(position = position_dodge(width = 10), size = 2) +
  geom_errorbar(aes(ymin = mean - lowbound,
                    ymax = mean + highbound), width = 0.3,
                position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 10)) +
  ylab("Sigma diff") +
  xlab("Number of mothers") +
  scale_colour_discrete(name = "Sigma diff") +
  scale_linetype_discrete(name = "Sigma diff") +
  scale_shape_discrete(name = "Eggs collected") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = -1) +
  scale_x_continuous(breaks = n_mothers)


# plotting ICC (intraclass correlation coefficient)
summary_ICC <- data.frame(comb = treatments_summ,
                          diff = NA,
                          n_mom = NA,
                          n_eggs = NA,
                          mean = NA,
                          SEM = NA,
                          SD = NA,
                          UHPDI = NA,
                          LHPDI = NA,
                          highbound = NA,
                          lowbound = NA
)

for(i in 1:length(treatments_summ)){
  data_temp <- subset(sim_result, comb == summary_ICC$comb[i])
  summary_ICC$diff[i] <- ifelse(data_temp$diff[1] == 0.5, "Mom > Res", ifelse(data_temp$diff[1] == 0, "Mom = Res", "Mom < Res"))
  summary_ICC$n_mom[i] <- data_temp$n_mom[1]
  summary_ICC$n_eggs[i] <- data_temp$n_eggs[1]
  summary_ICC$mean[i] <- mean(data_temp$sim_ICC, na.rm = TRUE)
  summary_ICC$SEM[i] <-sd(data_temp$sim_ICC, na.rm = TRUE) /
    sqrt(length(data_temp$sim_ICC))
  summary_ICC$SD[i] <- sd(data_temp$sim_ICC, na.rm = TRUE)
  summary_ICC$UHPDI[i] <- HPDI(data_temp$sim_ICC, prob = 0.95)[2]
  summary_ICC$LHPDI[i] <- HPDI(data_temp$sim_ICC, prob = 0.95)[1]
  summary_ICC$highbound[i] <- summary_ICC$UHPDI[i] - summary_ICC$mean[i]
  summary_ICC$lowbound[i] <- summary_ICC$mean[i] - summary_ICC$LHPDI[i]
}

g2 <- ggplot(data = summary_ICC, aes(x = n_mom,
                                     y = mean,
                                     colour = as.factor(diff),
                                     shape = as.factor(n_eggs),
                                     linetype = as.factor(diff))) +
  geom_point(position = position_dodge(width = 10)) +
  geom_errorbar(aes(ymin = mean - lowbound,
                    ymax = mean + highbound),
                width = 0.3,
                position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 10)) + 
  ylab("ICC") +
  xlab("Number of mothers") +
  scale_colour_discrete(name = "Sigma diff") +
  scale_linetype_discrete(name = "Sigma diff") +
  scale_shape_discrete(name = "Eggs collected") +
  geom_hline(yintercept = 0.5) +
  scale_x_continuous(breaks = n_mothers)



# plot differences between seed intercept and simulated output
summary_int <- data.frame(comb = treatments_summ,
                          diff = NA,
                          n_mom = NA,
                          n_eggs = NA,
                          mean = NA,
                          SEM = NA,
                          SD = NA,
                          UHPDI = NA,
                          LHPDI = NA,
                          highbound = NA,
                          lowbound = NA
)

for(i in 1:length(treatments_summ)){
  data_temp <- subset(sim_result, comb == summary_int$comb[i])
  summary_int$diff[i] <- ifelse(data_temp$diff[1] == 0.5, "Mom > Res", ifelse(data_temp$diff[1] == 0, "Mom = Res", "Mom < Res"))
  summary_int$n_mom[i] <- data_temp$n_mom[1]
  summary_int$n_eggs[i] <- data_temp$n_eggs[1]
  summary_int$mean[i] <- mean(data_temp$sim_int, na.rm = TRUE)
  summary_int$SEM[i] <-sd(data_temp$sim_int, na.rm = TRUE) /
    sqrt(length(data_temp$sim_int))
  summary_int$SD[i] <- sd(data_temp$sim_int, na.rm = TRUE)
  summary_int$UHPDI[i] <- HPDI(data_temp$sim_int, prob = 0.95)[2]
  summary_int$LHPDI[i] <- HPDI(data_temp$sim_int, prob = 0.95)[1]
  summary_int$highbound[i] <- summary_int$UHPDI[i] - summary_int$mean[i]
  summary_int$lowbound[i] <- summary_int$mean[i] - summary_int$LHPDI[i]
}

g3 <- ggplot(data = summary_int, aes(x = n_mom,
                                     y = mean,
                                     colour = as.factor(diff),
                                     shape = as.factor(n_eggs),
                                     linetype = as.factor(diff))) +
  geom_point(position = position_dodge(width = 10)) +
  geom_errorbar(aes(ymin = mean - lowbound,
                    ymax = mean + highbound),
                width = 0.3,
                position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 10)) +
  ylab("Intercept diff") +
  xlab("Number of mothers") +
  scale_colour_discrete(name = "Intercept diff") +
  scale_linetype_discrete(name = "Intercept diff") +
  scale_shape_discrete(name = "Eggs collected") +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = n_mothers)

gt <- plot_grid(g1, g2, g3, labels = c("A", "B", "C"), ncol = 1)  

save_plot("sim_results_1.png", gt,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 6
)

# so to answer the first question we can go through the variables we tested

# 1) how accurate is the intercept
# the intercept is very hard to specify when the difference is small,
# if the difference is large (if mothers have highly correlated egg volumes)
# a number of 60 to 100 mothers would be more than enough to safely estimate
# the intercept

# 2) how accurate is the maternal sd
# 3) how accurate is the population sd
# these two can be looked at at the same time - both graphs show that the estimation
# of standard deviations is close even with 30 mother IDs, but only with 60 do we start
# having some confidence that the 0.1 difference would be actually different


# Question 2 --------------------------------------------------------------

# 2) assuming several levels of correlation within mothers, how accurate is a bunched vs individual measurement when trying to differentiate two volumes?
# here we simulate both conditions (bunched and individual mothers) to see which provides a more accurate estimate
# for the individual mothers, we'll use the 10 egg collection from beforeand for bunched we'll simulate again

# bulding data frame to store simulated data
sim_bunched <- data.frame(
  sd_diff = NA,
  vol_diffs = NA,
  n_eggs = NA,
  n_mothers = NA,
  maternal_effects = NA,
  mother_ID = NA,
  measure = NA
)
sim_bunched <- na.omit(sim_bunched)

# we'll put different numbers of mothers in a single cage and randomly collect eggs from pool
# we'll collect either 50 or 100 per cage
n_eggs_cage <- c(50, 100)
# for each standard deviation difference
for(i in 1:length(sd_diffs)){
  
  # for each difference in volume
  for(j in 1:length(vol_diffs)){
    
    # for each number of mothers (equivalent to different cages)
    for(k in 1:length(n_mothers)){
      
      # for each collection number
      for(l in 1:length(n_eggs_cage)){
        
        # building a vector to store all the eggs
        egg_cage <- vector()
        
        # for each mother in the group
        for(m in 1:n_mothers[k]){
          # number of eggs mother k lays
          laid_eggs <- rpois(1, lambda = egg_lambda)
          
          # calculating the volume for each egg of this mother
          eggs_mother_k <-
            std_vol +
            vol_diffs[j] +
            rnorm(laid_eggs, mean = 0, sd = std_sd - sd_diffs[i]) +
            rnorm(laid_eggs, mean = 0, sd = std_sd + sd_diffs[i])
          
          # adding these eggs to the cage
          egg_cage <- c(egg_cage, eggs_mother_k)
        }
        
        sim_bunched_temp <- data.frame(
          sd_diff = sd_diffs[i],
          vol_diffs = vol_diffs[j],
          n_eggs_cage = n_eggs_cage[l],
          n_mothers = n_mothers[k],
          maternal_effects = F,
          mother_ID = NA,
          measure = sample(egg_cage, n_eggs_cage[l])
        )
        
        sim_bunched <- rbind(sim_bunched, sim_bunched_temp)
      }
    }
  }
}

# for each standard deviation difference
for(i in 1:length(sd_diffs)){
  
  # for each difference in volume
  for(j in 1:length(vol_diffs)){
    
    # for each number of mothers (equivalent to different cages)
    for(k in 1:length(n_mothers)){
      
      # for each mother in the group
      for(m in 1:n_mothers[k]){


        sim_bunched_temp <- data.frame(
          sd_diff = sd_diffs[i],
          vol_diffs = vol_diffs[j],
          n_eggs_cage = 10,
          n_mothers = n_mothers[k],
          maternal_effects = T,
          mother_ID = m,
          measure =
            std_vol +
            vol_diffs[j] +
            rnorm(10, mean = 0, sd = std_sd - sd_diffs[i]) +
            rnorm(10, mean = 0, sd = std_sd + sd_diffs[i])
        )
        
        # adding these eggs to the pile
        sim_bunched <- rbind(sim_bunched, sim_bunched_temp)
      }
    }
  }
}

# next we set up all combinations of parameters
reps <- length(sd_diffs) * length(n_mothers) * length(n_eggs_cage)
diff_values <- sort(rep(sd_diffs, reps / length(sd_diffs)))
n_mom_values <- rep(n_mothers, reps / length(n_mothers))
n_eggs_values <- rep(sort(rep(n_eggs_cage, reps / (length(n_eggs_cage) * length(sd_diffs)))), length(sd_diffs))

# and check if combinations are correctly set up
length(unique(paste(diff_values, n_mom_values, n_eggs_values, sep = "")))

# if so, build a list od models to store everything
models2 <- vector("list", length(diff_values))

# analyse each amount - takes a looong time
for(i in 1:length(models2)){
  models2[[i]] <- brm(
    measure ~ as.factor(vol_diffs),
    data = subset(sim_bunched,
                  n_mothers == n_mom_values[i] &
                    sd_diff == diff_values[i] &
                    n_eggs_cage == n_eggs_values[i] &
                    maternal_effects == FALSE),
    
    cores = 8,
    chains = 8,
    iter = 5000,
    warmup = 1000
  )
  
  print(paste("Now finished loop:", as.character(i), sep = " "))
}

# Save an object to a file
saveRDS(models2, file = "sim_models_example_2.rds")

# Restore the object
models2 <- readRDS(file = "sim_models_example_2.rds")

# next we set up all combinations of parameters
reps_ME <- length(sd_diffs) * length(n_mothers)
diff_values_ME <- sort(rep(sd_diffs, reps_ME / length(sd_diffs)))
n_mom_values_ME <- rep(n_mothers, reps_ME / length(n_mothers))

# and check if combinations are correctly set up
length(unique(paste(diff_values_ME, n_mom_values_ME, sep = "")))

# if so, build a list od models to store everything
models3 <- vector("list", length(diff_values_ME))

# analyse each amount - takes a looong time
for(i in 1:length(models3)){
  models3[[i]] <- brm(
    measure ~ as.factor(vol_diffs) + (1 | mother_ID),
    data = subset(sim_bunched,
                  n_mothers == n_mom_values[i] &
                    sd_diff == diff_values[i] &
                    maternal_effects == TRUE),
    
    cores = 8,
    chains = 8,
    iter = 5000,
    warmup = 1000,
    control = list(adapt_delta = 0.99)
  )
  
  print(paste("Now finished loop:", as.character(i), sep = " "))
}

# Save an object to a file
saveRDS(models3, file = "sim_models_example_3.rds")

# Restore the object
models3 <- readRDS(file = "sim_models_example_3.rds")

# making a data table for the output
sim_result_2 <- data.frame(
  diff = NA,
  n_mom = NA,
  n_eggs = NA,
  vol_diffs = NA,
  ind_mothers = NA,
  sim_int = NA
)
sim_result_2 <- na.omit(sim_result_2)


for(i in 1:length(models2)){
  interc_0 <- posterior_samples(models2[i], pars = "b_Intercept")
  interc_1 <- interc_0 + posterior_samples(models2[i], pars = "b_as.factorvol_diffs1")
  interc_2 <- interc_0 + posterior_samples(models2[i], pars = "b_as.factorvol_diffs2")
  diff <- diff_values[i]
  n_mom <- n_mom_values[i]
  n_egg <- n_eggs_values[i]
  vol_differences <- sort(rep(vol_diffs, length(interc_0[,1])))
  
  temp_diff <- data.frame(
    diff,
    n_mom,
    n_egg,
    vol_differences,
    ind_mothers = FALSE,
    rbind(interc_0, interc_1, interc_2)
  )
  
  colnames(temp_diff) <- c("diff", "n_mom", "n_eggs", "vol_diffs", "ind_mothers", "sim_int")
  sim_result_2 <- rbind(sim_result_2, temp_diff)
}

for(i in 1:length(models3)){
  interc_0 <- posterior_samples(models3[i], pars = "b_Intercept")
  interc_1 <- interc_0 + posterior_samples(models3[i], pars = "b_as.factorvol_diffs1")
  interc_2 <- interc_0 + posterior_samples(models3[i], pars = "b_as.factorvol_diffs2")
  diff <- diff_values_ME[i]
  n_mom <- n_mom_values_ME[i]
  n_egg <- 10
  vol_differences <- sort(rep(vol_diffs, length(interc_0[,1])))
  
  temp_diff <- data.frame(
    diff,
    n_mom,
    n_egg,
    vol_differences,
    ind_mothers = TRUE,
    rbind(interc_0, interc_1, interc_2)
  )
  
  colnames(temp_diff) <- c("diff", "n_mom", "n_eggs", "vol_diffs", "ind_mothers", "sim_int")
  sim_result_2 <- rbind(sim_result_2, temp_diff)
}

# plot differences between population sigma and mother sd
sim_result_2$comb <- paste(sim_result_2$diff,
                           sim_result_2$n_eggs,
                           sim_result_2$n_mom,
                           sim_result_2$vol_diffs,
                           sim_result_2$ind_mothers, sep = "")
treatments_summ <- unique(sim_result_2$comb)
summary <- data.frame(comb = treatments_summ,
                      diff = NA,
                      n_mom = NA,
                      n_eggs = NA,
                      vol_diffs = NA,
                      ind_mother = NA,
                      mean = NA,
                      SEM = NA,
                      SD = NA,
                      UHPDI = NA,
                      LHPDI = NA,
                      highbound = NA,
                      lowbound = NA
)

for(i in 1:length(treatments_summ)){
  data_temp <- subset(sim_result_2, comb == summary$comb[i])
  summary$diff[i] <- ifelse(data_temp$diff[1] == 0.5, "Mom > Res", ifelse(data_temp$diff[1] == 0, "Mom = Res", "Mom < Res"))
  summary$n_mom[i] <- data_temp$n_mom[1]
  summary$n_eggs[i] <- data_temp$n_eggs[1]
  summary$ind_mother[i] <-ifelse(data_temp$ind_mother[1] == T, "ME", "No ME")
  summary$vol_diffs[i] <- ifelse(data_temp$vol_diffs[1] == 0.5, "A", ifelse(data_temp$vol_diffs[1] == 1, "B", "C"))
  summary$mean[i] <- mean(data_temp$sim_int, na.rm = TRUE)
  summary$SEM[i] <-sd(data_temp$sim_int, na.rm = TRUE) /
    sqrt(length(data_temp$sim_int))
  summary$SD[i] <- sd(data_temp$sim_int, na.rm = TRUE)
  summary$UHPDI[i] <- HPDI(data_temp$sim_int, prob = 0.95)[2]
  summary$LHPDI[i] <- HPDI(data_temp$sim_int, prob = 0.95)[1]
  summary$highbound[i] <- summary$UHPDI[i] - summary$mean[i]
  summary$lowbound[i] <- summary$mean[i] - summary$LHPDI[i]
}

hline_05 <- data.frame(
  diff = rep(c("Mom > Res", "Mom = Res", "Mom < Res"), 2),
  ind_mother = rep(c("ME", "No ME"), 3),
  five = 8.5,
  nine = 9,
  ten = 10
)

g4 <- ggplot(data = summary, aes(x = n_mom,
                                 y = mean,
                                 colour = as.factor(vol_diffs),
                                 shape = as.factor(n_eggs),
                                 linetype = as.factor(vol_diffs))) +
  geom_point(position = position_dodge(width = 10), size = 2) +
  geom_errorbar(aes(ymin = mean - lowbound,
                    ymax = mean + highbound), width = 0.3,
                position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 10)) +
  ylab("Volume") +
  xlab("Number of mothers") +
  scale_colour_discrete(name = "Treatment") +
  scale_linetype_discrete(name = "Treatment") +
  scale_shape_discrete(name = "Eggs collected") +
  geom_hline(data = hline_05, aes(yintercept = five)) +
  geom_hline(data = hline_05, aes(yintercept = nine)) +
  geom_hline(data = hline_05, aes(yintercept = ten)) +
  scale_x_continuous(breaks = n_mothers) +
  facet_grid(vars(as.factor(diff)), vars(as.factor(ind_mother)))

save_plot("sim_results_2.png", g4,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 10
)



# Example plots -----------------------------------------------------------

# same sd for both errors
sd_mother <- 0.3
sd_egg <- 0.3

# 10 mothers with 10 eggs each
n_eggs_example <- 10
n_moms_example <- 10

# start a data frame to fill with random data
example_data <- data.frame(
  mom_ID = NA,
  egg_ID = NA,
  measure = NA,
  egg_number = NA
)

# simulate data
for(i in 1:n_moms_example){
    mother_mean <- rnorm(1, mean = 8, sd = sd_mother)
  
    temp_example <- data.frame(
      mom_ID = as.factor(i),
      egg_ID = seq(1,10,1),
      measure = mother_mean + rnorm(n_eggs_example, mean = 0, sd = sd_egg),
      egg_number = NA
    )
    example_data <- rbind(example_data, temp_example)
}

# add a number to each egg and a different sequence to plot in a random order
example_data$egg_number <- seq(1,length(example_data[,1])) - 1
example_data$random <- sample(example_data$egg_number)

# remove NAs and calculate distance to the mean of the population
example_data <- na.omit(example_data)
example_data$error <-example_data$measure - mean(example_data$measure)

# plot distances to population
g5 <- ggplot(data = example_data, aes(x = random, y = measure)) +
  geom_point() +
  geom_hline(yintercept = mean(example_data$measure)) +
  geom_errorbar(data = example_data, aes(ymin = measure - error, ymax = measure), width = 0.1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Volume")

# save plot
save_plot("sim_example_pop.png", g5,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 10
)

# a vector with all unique mother IDs
moms <- unique(example_data$mom_ID)

# data frame to fill with means
example_mom_means <- data.frame(
  mom_ID = rep(1,length(moms)),
  mea = rep(1,length(moms))
)

# means of each mother
for(i in 1:length(moms)){
  temp_example <- subset(example_data, mom_ID == moms[i])
  example_mom_means$mom_ID[i] <- temp_example$mom_ID[1]
  example_mom_means$mea[i] <- mean(temp_example$measure)
}

# calculate the distance of each egg from the mean of the respective mother
for(i in 1:nrow(example_data)){
  example_data$error_mom[i] <- example_data$measure[i] - subset(example_mom_means, mom_ID == example_data$mom_ID[i])$mea
}

# plot distances to mother mean
g6 <- ggplot(data = example_data, aes(x = egg_ID, y = measure)) +
  geom_point() +
  geom_hline(data = example_mom_means, aes(yintercept = mea)) +
  geom_hline(yintercept = mean(example_data$measure), linetype = "dashed") +
  geom_errorbar(data = example_data, aes(ymin = measure - error_mom, ymax = measure), width = 0.1) +
  facet_grid(~mom_ID, space = 'free_x', scales = 'free_x', switch = 'x') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Volume")

# save plot
save_plot("sim_example_mom.png", g6,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 10
)
