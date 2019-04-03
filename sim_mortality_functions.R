library(ggplot2)
library(cowplot)

#### juvenile ####
juv_optimum <- 22

juvenileMortality <- function(envValue){
  baseMort <- 0.01
  misMatch <- juv_optimum - envValue
  deathRate <- ifelse(misMatch > 0, 10, 6)
  
  mortality <- exp(-((envValue - juv_optimum) ^ 2)/(2 * deathRate ^ 2)) #baseMort + deathRate * misMatch ^ 2 
  mortality <- ifelse(mortality > 1, 1, mortality)
  mortality <- ifelse(mortality < baseMort, baseMort, mortality)
  return(mortality)
}

juv <- data.frame(env = seq(5,35,0.01), mort = NA)

juv$mort <- juvenileMortality(juv$env)

ggplot(juv, aes(x = env, y = mort)) + geom_point(size = 0.1) + geom_vline(xintercept = juv_optimum) +
  annotate("text", y = 0.8, x = juv_optimum - 13, label = "1 - (0.01 + dRate x (Temp - 22)^2)") +
  annotate("text", y = 0.5, x = juv_optimum - 3, label = "dRate = 0.003") +
  annotate("text", y = 0.5, x = juv_optimum + 3, label = "dRate = 0.01") +
  ylab("Survival probability") + xlab("Temperature") +
  scale_x_continuous("Temperature", breaks = seq(5,35,5))

#### adult ####
adu_optimum <- 22
adultMortality <- function(envValue, age){
  baseMort <- 0.01
  misMatch <- adu_optimum - envValue
  deathRate <- ifelse(misMatch > 0, 25, 8)
    
  mortality <- exp(-((envValue - adu_optimum) ^ 2)/(2 * deathRate ^ 2)) * (1 - 0.05)^age #(baseMort + deathRate * misMatch ^ 2)
  mortality <- ifelse(mortality > 1, 1, mortality)
  mortality <- ifelse(mortality < baseMort, baseMort, mortality)
  return(mortality)
}

env <- seq(0,40,0.01)

age <- NA
for (i in 1:6){
  tempage <- as.data.frame(rep(i, length(env)))
  age <- rbind(age,tempage)
}

colnames(age) <- c("age")
age <- na.omit(age)

adu <- data.frame(env = rep(env,6), age = age, mort = NA)

adu$mort <- adultMortality(adu$env,adu$age)

ggplot(adu, aes(x = env, y = mort, colour = as.factor(age))) + geom_point(size = 0.1) + geom_vline(xintercept = adu_optimum) +
  annotate("text", y = 0.1, x = 10, label = "1 - (0.01 + (dRate x (Temp - 22)^2) * 1.05^age)") +
  annotate("text", y = 0.5, x = adu_optimum-4, label = "dRate = 0.001") +
  annotate("text", y = 0.5, x = adu_optimum+4, label = "dRate = 0.004") +
  ylab("Survival probability") + xlab("Temperature") +
  scale_x_continuous("Temperature", breaks = seq(0,40,5)) +
  labs(colour = "Age")
  
#### Fertility ####
fer_optimum <- 28
gaus <- function(envValue){
  stddev <- ifelse(fer_optimum - envValue > 0, 8, 4)
  fert <- exp(-((envValue - fer_optimum) ^ 2.0)/(2.0 * stddev ^ 2)) * 150.0;
  return(fert)
}
fertility <- function(envValue, age){
  fert <- gaus(envValue) * (1 - 0.05)^age
  return(fert)
}

env <- seq(0,40,0.01)

age <- NA
for (i in 1:6){
  tempage <- as.data.frame(rep(i, length(env)))
  age <- rbind(age,tempage)
}

colnames(age) <- c("age")
age <- na.omit(age)

fer <- data.frame(env = rep(env,6), age = age, fertility = NA)

fer$fertility <- fertility(fer$env,fer$age)

ggplot(fer, aes(x = env, y = fertility, colour = as.factor(age))) + geom_point(size = 0.1) + geom_vline(xintercept = fer_optimum) +
  annotate("text", y = 100, x = 10, label = "exp(-((Temp - 28)^2)/(2 * sd^2))\n* 150.0 * 1.05^age)") +
  annotate("text", y = 20, x = fer_optimum-4, label = "sd = 8") +
  annotate("text", y = 20, x = fer_optimum+4, label = "sd = 4") +
  ylab("Number of eggs") + xlab("Temperature") +
  scale_x_continuous("Temperature", breaks = seq(0,40,5)) +
  labs(colour = "Age")
