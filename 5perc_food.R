setwd("~/Collect_data/Exp_5perc_food")

library(rethinking)
library(brms)

d <- read.csv("Egg collection.csv")

d$Adults <- d$Male + d$Female
d$check <- d$Adults - d$Eggs
d$comb <- paste(d$Mother.temperature,d$Offspring.temperature,d$Offspring.food, sep = "")

d <- subset(d, Offspring.food == "R")

x <- unique(levels(d$ID))
y <- unlist(lapply(x, function(x) substr(x, 1, 1)))
cor <- data.frame(id = x, treat = y, m_egg = NA, m_adul = NA,  mm_egg = NA, mm_adul = NA)
for(i in 1:length(cor$id)){
  temp_m <- subset(d, ID == as.character(cor$id[i]) & Offspring.temperature == as.character(cor$treat[i]))
  temp_mm <- subset(d, ID == as.character(cor$id[i]) & Offspring.temperature != as.character(cor$treat[i]))
  cor$m_egg[i] = ifelse(length(temp_m$Adults) == 0, NA, temp_m$Eggs)
  cor$m_adul[i] = ifelse(length(temp_m$Adults) == 0, NA, temp_m$Adults)
  cor$mm_egg[i] = ifelse(length(temp_mm$Adults) == 0, NA, temp_mm$Eggs)
  cor$mm_adul[i] = ifelse(length(temp_mm$Adults) == 0, NA, temp_mm$Adults)
}
cor <- na.omit(cor)

match <- bf(m_adul | trials(m_egg) ~ treat + (1 | p | id)) + binomial()
mismatch <- bf(mm_adul | trials(mm_egg) ~ treat + (1 | p | id)) + binomial()

mm1 <- brm(match + mismatch, data = cor, cores = 4, chains = 4, iter = 5000, warmup = 1000)

h_mm <- c("madul_Intercept - (madul_Intercept + madul_treatH) < 0",
          "madul_Intercept - mmadul_Intercept < 0",
          "madul_Intercept - (mmadul_Intercept + mmadul_treatH) < 0",
          "(madul_Intercept + madul_treatH) - mmadul_Intercept < 0",
          "(madul_Intercept + madul_treatH) - (mmadul_Intercept + mmadul_treatH) < 0",
          "mmadul_Intercept - (mmadul_Intercept + mmadul_treatH) < 0")
# testing them
hypothesis(mm1, h_mm)

c_cor <- c("id__madul_Intercept__mmadul_Intercept > 0")

hypothesis(mm1, c_cor, class = "cor")
