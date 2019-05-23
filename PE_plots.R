source("~/Code/Rscripts/plotPE.R")

#### NoPEvsPE_separate ####
setwd("~/Collect_data/NoPEvsPE_separate")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:1){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}


n <- 0
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_0 <- read.csv(dist_num, header = FALSE)

colnames(d_0) <- c("replicate", "gen", "dist", "choice")

d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
d_0$choice <- ifelse(d_0$choice > 6, d_0$choice - 6, d_0$choice)
d_0$choice <- as.factor(d_0$choice)
d_0$replicate <- as.factor(d_0$replicate)


geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_color_discrete(name = c("Genotype"),
                       breaks = c("0", "1","2","3","4","5",
                                  "6","7","8","9","10","11"),
                       labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                  "NoPE 3","NoPE 4", "NoPE 5",
                                  "PE 0", "PE 1","PE 2",
                                  "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_0 <- read.csv(val_num, header = FALSE)

g_0$diff <- g_0[,6]-g_0[,4]
colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_0[c(1,2,3)] #g
y <- g_0[c(1,2,4)] #p
w <- g_0[c(1,2,5)] #pe
z <- g_0[c(1,2,6)] #e
a_0 <- g_0[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_0$type <- "diff"

g_0 <- rbind(x,y,z)
g_0 <- subset(g_0,  gen > 1)

g_0$type <- as.factor(g_0$type)

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())


n <- 1
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_1 <- read.csv(dist_num, header = FALSE)

colnames(d_1) <- c("replicate", "gen", "dist", "choice")

d_1$PE <- ifelse(d_1$choice < 6, "No PE", "PE")
d_1$choice <- ifelse(d_1$choice > 5, d_1$choice - 6, d_1$choice)
d_1$choice <- as.factor(d_1$choice)
d_1$replicate <- as.factor(d_1$replicate)


geno_1 <- ggplot(data = d_1, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_color_discrete(name = c("Genotype"),
                       breaks = c("0", "1","2","3","4","5",
                                  "6","7","8","9","10","11"),
                       labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                  "NoPE 3","NoPE 4", "NoPE 5",
                                  "PE 0", "PE 1","PE 2",
                                  "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_1 <- read.csv(val_num, header = FALSE)

g_1$diff <- g_1[,6]-g_1[,4]
colnames(g_1) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_1[c(1,2,3)] #g
y <- g_1[c(1,2,4)] #p
w <- g_1[c(1,2,5)] #pe
z <- g_1[c(1,2,6)] #e
a_1 <- g_1[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_1$type <- "diff"

g_1 <- rbind(x,y,z)
g_1 <- subset(g_1,  gen > 1)

g_1$type <- as.factor(g_1$type)

val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

NoPEvsPE <- plot_grid(geno_0, geno_1, val_0,val_1, labels = c("A", "B", "C", "D"), ncol = 2)

save_plot("NoPEvsPE.png", NoPEvsPE,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)

a_0 <- subset(a_0,  gen > 1)
a_0$type <- as.factor(a_0$type)

diff_0 <- ggplot(data = a_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

a_1 <- subset(a_1,  gen > 1)
a_1$type <- as.factor(a_1$type)

diff_1 <- ggplot(data = a_1, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

NoPEvsPE_diff <- plot_grid(geno_0, geno_1, val_0, val_1, diff_0, diff_1,
                           labels = c("A", "B", "C", "D", "E", "F"), ncol = 2)

save_plot("NoPEvsPE_diff.png", NoPEvsPE_diff,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 8
)

#### NoPEvsPE_evol_stab ####

setwd("~/Collect_data/NoPEvsPE_evol_stab")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

rsamples <- sampleDist(1, 1, 10)
n <- 0
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_0 <- read.csv(dist_num, header = FALSE)

colnames(d_0) <- c("replicate", "gen", "dist", "choice")

d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
d_0$choice <- as.factor(d_0$choice)
d_0 <- d_0[d_0$replicate %in% rsamples, ]
d_0$replicate <- as.factor(d_0$replicate)


geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 2) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_fill_discrete(name = c("Genotype"),
                       breaks = c("0", "1","2","3","4","5",
                                  "6","7","8","9","10","11"),
                       labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                  "NoPE 3","NoPE 4", "NoPE 5",
                                  "PE 0", "PE 1","PE 2",
                                  "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_0 <- read.csv(val_num, header = FALSE)

g_0$diff <- g_0[,6]-g_0[,4]
colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_0[c(1,2,3)] #g
y <- g_0[c(1,2,4)] #p
w <- g_0[c(1,2,5)] #pe
z <- g_0[c(1,2,6)] #e
a_0 <- g_0[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_0$type <- "diff"

g_0 <- rbind(x,y,z)
g_0 <- subset(g_0,  gen > 1)
g_0 <- g_0[g_0$replicate %in% rsamples, ]
g_0$type <- as.factor(g_0$type)

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  scale_colour_discrete(name = c("Values"),
                      breaks = c("env", "gene", "pheno"),
                      labels = c("Env", "Avg Gene", "Avg Pheno")) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylim(10,30)

NoPEvsPE_evolve <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)

save_plot("NoPEvsPE_evolve.png", NoPEvsPE_evolve,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 6
)

#### NoPEvsPE_evol_var2 ####

setwd("~/Collect_data/NoPEvsPE_evol_var2")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

rsamples <- sampleDist(1, 1, 10)
n <- 0
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_0 <- read.csv(dist_num, header = FALSE)

colnames(d_0) <- c("replicate", "gen", "dist", "choice")

d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
d_0$choice <- as.factor(d_0$choice)
d_0 <- d_0[d_0$replicate %in% rsamples, ]
d_0$replicate <- as.factor(d_0$replicate)


geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 2) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_fill_discrete(name = c("Genotype"),
                      breaks = c("0", "1","2","3","4","5",
                                 "6","7","8","9","10","11"),
                      labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                 "NoPE 3","NoPE 4", "NoPE 5",
                                 "PE 0", "PE 1","PE 2",
                                 "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_0 <- read.csv(val_num, header = FALSE)

g_0$diff <- g_0[,6]-g_0[,4]
colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_0[c(1,2,3)] #g
y <- g_0[c(1,2,4)] #p
w <- g_0[c(1,2,5)] #pe
z <- g_0[c(1,2,6)] #e
a_0 <- g_0[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_0$type <- "diff"

g_0 <- rbind(x,y,z)
g_0 <- subset(g_0,  gen > 1)

g_0 <- g_0[g_0$replicate %in% rsamples, ]
g_0$type <- as.factor(g_0$type)
g_0$type <- relevel(g_0$type, "gene")
g_0$type <- relevel(g_0$type, "pheno")
g_0$type <- relevel(g_0$type, "env")

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  scale_colour_manual(name = c("Values"),
                      values = c("#F8766D", "#619DFF", "#00BA38"),
                        breaks = c("env", "gene", "pheno"),
                        labels = c("Env", "Avg Gene", "Avg Pheno")) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylim(10,30)

NoPEvsPE_var2 <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)

save_plot("NoPEvsPE_evol_var2.png", NoPEvsPE_var2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 6
)

#### NoPEvsPE_evol_var10 ####

setwd("~/Collect_data/NoPEvsPE_evol_var10")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

rsamples <- sampleDist(1, 1, 10)
n <- 0
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_0 <- read.csv(dist_num, header = FALSE)

colnames(d_0) <- c("replicate", "gen", "dist", "choice")

d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
d_0$choice <- as.factor(d_0$choice)
d_0 <- d_0[d_0$replicate %in% rsamples, ]
d_0$replicate <- as.factor(d_0$replicate)


geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 2) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_fill_discrete(name = c("Genotype"),
                      breaks = c("0", "1","2","3","4","5",
                                 "6","7","8","9","10","11"),
                      labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                 "NoPE 3","NoPE 4", "NoPE 5",
                                 "PE 0", "PE 1","PE 2",
                                 "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_0 <- read.csv(val_num, header = FALSE)

g_0$diff <- g_0[,6]-g_0[,4]
colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_0[c(1,2,3)] #g
y <- g_0[c(1,2,4)] #p
w <- g_0[c(1,2,5)] #pe
z <- g_0[c(1,2,6)] #e
a_0 <- g_0[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_0$type <- "diff"

g_0 <- rbind(x,y,z)
g_0 <- subset(g_0,  gen > 1)

g_0 <- g_0[g_0$replicate %in% rsamples, ]
g_0$type <- as.factor(g_0$type)
g_0$type <- relevel(g_0$type, "gene")
g_0$type <- relevel(g_0$type, "pheno")
g_0$type <- relevel(g_0$type, "env")

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  scale_colour_manual(name = c("Values"),
                        values = c("#F8766D", "#619DFF", "#00BA38"),
                        breaks = c("env", "gene", "pheno"),
                        labels = c("Env", "Avg Gene", "Avg Pheno")) +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

NoPEvsPE_var10 <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)

save_plot("NoPEvsPE_evol_var10.png", NoPEvsPE_var10,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 6
)

#### NoPEvsPE_block3 ####

setwd("~/Collect_data/NoPEvsPE_block3")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

n <- 0
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_0 <- read.csv(dist_num, header = FALSE)

colnames(d_0) <- c("replicate", "gen", "dist", "choice")

d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
d_0$choice <- as.factor(d_0$choice)
d_0$replicate <- as.factor(d_0$replicate)


geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_color_discrete(name = c("Genotype"),
                       breaks = c("0", "1","2","3","4","5",
                                  "6","7","8","9","10","11"),
                       labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                  "NoPE 3","NoPE 4", "NoPE 5",
                                  "PE 0", "PE 1","PE 2",
                                  "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_0 <- read.csv(val_num, header = FALSE)

g_0$diff <- g_0[,6]-g_0[,4]
colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_0[c(1,2,3)] #g
y <- g_0[c(1,2,4)] #p
w <- g_0[c(1,2,5)] #pe
z <- g_0[c(1,2,6)] #e
a_0 <- g_0[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_0$type <- "diff"

g_0 <- rbind(x,y,z)
g_0 <- subset(g_0,  gen > 1)

g_0$type <- as.factor(g_0$type)
g_0$type <- relevel(g_0$type, "gene")
g_0$type <- relevel(g_0$type, "pheno")
g_0$type <- relevel(g_0$type, "env")

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

NoPEvsPE_block3 <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)

save_plot("NoPEvsPE_block3.png", NoPEvsPE_block3,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 6
)

#### NoPEvsPE_block2 ####

setwd("~/Collect_data/NoPEvsPE_block2")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

n <- 0
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_0 <- read.csv(dist_num, header = FALSE)

colnames(d_0) <- c("replicate", "gen", "dist", "choice")

d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
d_0$choice <- as.factor(d_0$choice)
d_0$replicate <- as.factor(d_0$replicate)


geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_color_discrete(name = c("Genotype"),
                       breaks = c("0", "1","2","3","4","5",
                                  "6","7","8","9","10","11"),
                       labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                  "NoPE 3","NoPE 4", "NoPE 5",
                                  "PE 0", "PE 1","PE 2",
                                  "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_0 <- read.csv(val_num, header = FALSE)

g_0$diff <- g_0[,6]-g_0[,4]
colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_0[c(1,2,3)] #g
y <- g_0[c(1,2,4)] #p
w <- g_0[c(1,2,5)] #pe
z <- g_0[c(1,2,6)] #e
a_0 <- g_0[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_0$type <- "diff"

g_0 <- rbind(x,y,z)
g_0 <- subset(g_0,  gen > 1)

g_0$type <- as.factor(g_0$type)
g_0$type <- relevel(g_0$type, "gene")
g_0$type <- relevel(g_0$type, "pheno")
g_0$type <- relevel(g_0$type, "env")

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

NoPEvsPE_block2 <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)

save_plot("NoPEvsPE_block2.png", NoPEvsPE_block2,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 6
)

#### NoPEvsPE_block4_5_10 ####

setwd("~/Collect_data/NoPEvsPE_block4_5_10")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:2){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

for(n in 0:2){
  #n <- 0
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")
  
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- ifelse(n == 0, 4, ifelse(n == 1, 5, 10))
  png_name <- paste("NoPEvsPE_block", as.character(block), ".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

#### PE_block3 ####

setwd("~/Collect_data/PE_block3")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

n <- 0
dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
d_0 <- read.csv(dist_num, header = FALSE)

colnames(d_0) <- c("replicate", "gen", "dist", "choice")

d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
d_0$choice <- as.factor(d_0$choice)
d_0$replicate <- as.factor(d_0$replicate)


geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
  geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
  ylab("Distribution") + xlab("Generation") +
  scale_color_discrete(name = c("Genotype"),
                       breaks = c("0", "1","2","3","4","5",
                                  "6","7","8","9","10","11"),
                       labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                  "NoPE 3","NoPE 4", "NoPE 5",
                                  "PE 0", "PE 1","PE 2",
                                  "PE 3","PE 4", "PE 5")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
g_0 <- read.csv(val_num, header = FALSE)

g_0$diff <- g_0[,6]-g_0[,4]
colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")


x <- g_0[c(1,2,3)] #g
y <- g_0[c(1,2,4)] #p
w <- g_0[c(1,2,5)] #pe
z <- g_0[c(1,2,6)] #e
a_0 <- g_0[c(1,2,7)] #diff

x$type <- "gene"
y$type <- "pheno"
w$type <- "peff"
z$type <- "env"
a_0$type <- "diff"

g_0 <- rbind(x,y,z)
g_0 <- subset(g_0,  gen > 1)

g_0$type <- as.factor(g_0$type)
g_0$type <- relevel(g_0$type, "gene")
g_0$type <- relevel(g_0$type, "pheno")
g_0$type <- relevel(g_0$type, "env")

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

PE_block3 <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)

save_plot("PE_block3.png", PE_block3,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3,
          base_height = 6
)

#### PE_block2_4_5_10 ####

setwd("~/Collect_data/PE_block2_4_5_10")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:3){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

for(n in 0:3){
  #n <- 0
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  colnames(g_0) <- c("replicate", "gen", "value", "value", "value", "value","value")
  
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- ifelse(n == 0, 2, ifelse(n == 1, 4, ifelse(n == 2, 5, 10)))
  png_name <- paste("NoPEvsPE_block", as.character(block), ".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}
#### PE_mut_blocks ####

setwd("~/Collect_data/PE_mut_blocks")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:4){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

for(n in 0:4){
  #n <- 0
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  d_0 <- subset(d_0, gen < 31)
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  g_0 <- subset(g_0, V2 < 31)
  
  colnames(g_0) <- c("replicate", "gen",
                     "value", "value", "value", "value","value")
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- ifelse(n == 0, 2,
                  ifelse(n == 1, 3,
                         ifelse(n == 2, 4,
                                ifelse(n == 3, 5,
                                       10))))
  png_name <- paste("PE_block", as.character(block), ".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

for(n in 0:4){
  #n <- 0
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)

  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]

  colnames(g_0) <- c("replicate", "gen",
                     "value", "value", "value", "value","value")
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- ifelse(n == 0, 2,
                  ifelse(n == 1, 3,
                         ifelse(n == 2, 4,
                                ifelse(n == 3, 5,
                                       10))))
  png_name <- paste("PE_block", as.character(block), "_allgen.png", sep = "")
  
  save_plot(png_name, NoPEvsPE_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

#### NoPE_mut_blocks ####

setwd("~/Collect_data/NoPE_mut_blocks")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:4){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

for(n in 0:4){
  #n <- 0
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  d_0 <- subset(d_0, gen < 31)
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  g_0 <- subset(g_0, V2 < 31)
  
  colnames(g_0) <- c("replicate", "gen",
                     "value", "value", "value", "value","value")
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- ifelse(n == 0, 2,
                  ifelse(n == 1, 3,
                         ifelse(n == 2, 4,
                                ifelse(n == 3, 5,
                                       10))))
  png_name <- paste("NoPE_block", as.character(block), ".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

#### NoPEvsPE_blocks1-3 ####

setwd("~/Collect_data/NoPEvsPE_blocks1-3")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:35){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- sampleDist(5, 1, 100)

for(n in 0:35){
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  d_0 <- d_0[d_0$replicate %in% rsamples, ]
  d_0 <- subset(d_0, gen < 31)
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  g_0 <- subset(g_0, V2 < 31)
  
  colnames(g_0) <- c("replicate", "gen",
                     "value", "value", "value", "value","value")
  g_0$replicate <- as.factor(g_0$replicate)
  g_0 <- g_0[g_0$replicate %in% rsamples, ]
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- round(as.numeric(as.character(p$iInterval[n + 1])))
  choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
  init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
  
  png_name <- paste("NoPE_block", as.character(block),"_",init,"_choice",choice,".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

#### NoPEvsPE_evolve_blocks1-3 ####

setwd("~/Collect_data/NoPEvsPE_evolve_blocks1-3")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:17){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- sampleDist(5, 1, 100)

for(n in 0:17){
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  d_0 <- d_0[d_0$replicate %in% rsamples, ]
  #d_0 <- subset(d_0, gen < 31)
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  #g_0 <- subset(g_0, V2 < 31)
  
  colnames(g_0) <- c("replicate", "gen",
                     "value", "value", "value", "value","value")
  g_0$replicate <- as.factor(g_0$replicate)
  g_0 <- g_0[g_0$replicate %in% rsamples, ]
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_evolve_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- round(as.numeric(as.character(p$iInterval[n + 1])))
  choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
  init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
  
  png_name <- paste("NoPE_block", as.character(block),"_",init,"_choice",choice,".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_evolve_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

#### NoPEvsPE_evolve_blocks1-10 ####

setwd("~/Collect_data/NoPEvsPE_evolve_blocks1-10")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:35){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- sampleDist(5, 1, 100)

for(n in 0:35){
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  d_0 <- d_0[d_0$replicate %in% rsamples, ]
  d_0 <- subset(d_0, gen < 31)
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  g_0 <- subset(g_0, V2 < 31)
  
  colnames(g_0) <- c("replicate", "gen",
                     "value", "value", "value", "value","value")
  g_0$replicate <- as.factor(g_0$replicate)
  g_0 <- g_0[g_0$replicate %in% rsamples, ]
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_evolve_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- round(as.numeric(as.character(p$iInterval[n + 1])))
  choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
  init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
  
  png_name <- paste("NoPE_block", as.character(block),"_",init,"_choice",choice,".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_evolve_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

#### Block_evol ####

setwd("~/Collect_data/Block_evol")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:95){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- sampleDist(1, 1, 100)

for(n in 0:95){
  dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
  d_0 <- read.csv(dist_num, header = FALSE)
  
  colnames(d_0) <- c("replicate", "gen", "dist", "choice")
  
  d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
  d_0$choice <- as.factor(d_0$choice)
  d_0$replicate <- as.factor(d_0$replicate)
  d_0 <- d_0[d_0$replicate %in% rsamples, ]
  d_0 <- subset(d_0, gen < 31)
  
  geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
    geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
    ylab("Distribution") + xlab("Generation") +
    scale_color_discrete(name = c("Genotype"),
                         breaks = c("0", "1","2","3","4","5",
                                    "6","7","8","9","10","11"),
                         labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                    "NoPE 3","NoPE 4", "NoPE 5",
                                    "PE 0", "PE 1","PE 2",
                                    "PE 3","PE 4", "PE 5")) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
  g_0 <- read.csv(val_num, header = FALSE)
  
  g_0$diff <- g_0[,6]-g_0[,4]
  g_0 <- subset(g_0, V2 < 31)
  
  colnames(g_0) <- c("replicate", "gen",
                     "value", "value", "value", "value","value")
  g_0$replicate <- as.factor(g_0$replicate)
  g_0 <- g_0[g_0$replicate %in% rsamples, ]
  
  x <- g_0[c(1,2,3)] #g
  y <- g_0[c(1,2,4)] #p
  w <- g_0[c(1,2,5)] #pe
  z <- g_0[c(1,2,6)] #e
  a_0 <- g_0[c(1,2,7)] #diff
  
  x$type <- "gene"
  y$type <- "pheno"
  w$type <- "peff"
  z$type <- "env"
  a_0$type <- "diff"
  
  g_0 <- rbind(x,y,z)
  g_0 <- subset(g_0,  gen > 1)
  
  g_0$type <- as.factor(g_0$type)
  g_0$type <- relevel(g_0$type, "gene")
  g_0$type <- relevel(g_0$type, "pheno")
  g_0$type <- relevel(g_0$type, "env")
  
  val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
    geom_line() + facet_grid(rows = vars(replicate)) +
    ylab("value") + xlab("Generation") +
    theme_minimal() +
    scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  NoPEvsPE_evolve_blocks <- plot_grid(geno_0, val_0, labels = c("A", "B"), ncol = 1)
  
  block <- round(as.numeric(as.character(p$iInterval[n + 1])))
  choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
  init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
  mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n + 1])) == 0, "nomut", "mut")
  
  png_name <- paste("Block", as.character(block),"_",init,"_",mut,"_choice",choice,".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_evolve_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}

#### Overwinter ####

setwd("~/Collect_data/Overwinter")

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:3){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- sampleDist(1, 1, 1)
for(interv in 1:2){
  if(interv == 1){
    intstart <- 1
    intend <- 101}
  else{
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  for(n in 0:3){
    dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
    d_0 <- read.csv(dist_num, header = FALSE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = PE)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
      )
    geno_0
    
    val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
    g_0 <- read.csv(val_num, header = FALSE)
    
    g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    y <- g_0[c(1,2,4)] #p
    w <- g_0[c(1,2,5)] #pe
    z <- g_0[c(1,2,6)] #e
    a_0 <- g_0[c(1,2,10)] #diff
    t <- g_0[c(1,2,7)] #threshold
    o <- g_0[c(1,2,8)] #overwinter
    ag <- g_0[c(1,2,9)] #age
    
    x$type <- "gene"
    y$type <- "pheno"
    w$type <- "peff"
    z$type <- "env"
    a_0$type <- "diff"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    g_0 <- rbind(x,y,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "gene")
    g_0$type <- relevel(g_0$type, "pheno")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank(),
            strip.text = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,ag)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank(),
            strip.text = element_blank())
    
    #val_1
    NoPEvsPE_evolve_blocks <- plot_grid(geno_0, val_0, val_1, ncol = 1)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n + 1]))
    choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
    init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    #mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n + 1])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/init_",init,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    else png_name <- paste("total/init_",init,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 1, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.3,
              base_height = 6
    )
  }
}
#### Fitnessfactor ####

setwd("~/Collect_data/Fitnessfactor")
num_samples <-16

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:(num_samples - 1)){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- sampleDist(3, 1, 10)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 0:(num_samples - 1)){
    dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
    d_0 <- read.csv(dist_num, header = FALSE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    geno_0
    
    val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
    g_0 <- read.csv(val_num, header = FALSE)
    
    g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    y <- g_0[c(1,2,4)] #p
    w <- g_0[c(1,2,5)] #pe
    z <- g_0[c(1,2,6)] #e
    a_0 <- g_0[c(1,2,10)] #diff
    t <- g_0[c(1,2,7)] #threshold
    o <- g_0[c(1,2,8)] #overwinter
    ag <- g_0[c(1,2,9)] #age
    
    x$type <- "gene"
    y$type <- "pheno"
    w$type <- "peff"
    z$type <- "env"
    a_0$type <- "diff"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    g_0 <- rbind(x,y,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "gene")
    g_0$type <- relevel(g_0$type, "pheno")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,ag)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    NoPEvsPE_evolve_blocks <- plot_grid(geno_0, val_0, val_1, ncol = 1)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n + 1]))
    choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n + 1])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 1, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.3,
              base_height = 6
    )
  }
}

#### Fitness factor function ####

Tm <- 32
Psi <- 8.3
Rho <- 0.0
deltaT <- 5.1
Sclfactor <- 49.505

loganFactor <- function(x, y){
  tau <- (32 - x)/5.1
  return((8.3/49.505 * (exp(0.1 * x) - exp(0.1 * 32 - tau))) * (1-0.995)^y)
}

factor <- data.frame(temp = seq(10,40,0.001))
factor$multiplier <- sapply(factor$temp, function(x, y) loganFactor(x, 0))
factor$multiplier1 <- sapply(factor$temp, function(x, y) loganFactor(x, 1))
factor$multiplier2 <- sapply(factor$temp, function(x, y) loganFactor(x, 2))
factor$multiplier3 <- sapply(factor$temp, function(x, y) loganFactor(x, 3))
factor$multiplier4 <- sapply(factor$temp, function(x, y) loganFactor(x, 4))
factor$multiplier11 <- sapply(factor$temp, function(x, y) loganFactor(x, 11))

factor$multiplier <- ifelse(factor$multiplier < 0, 0, factor$multiplier)
factor$multiplier1 <- ifelse(factor$multiplier1 < 0, 0, factor$multiplier1)
factor$multiplier2 <- ifelse(factor$multiplier2 < 0, 0, factor$multiplier2)
factor$multiplier3 <- ifelse(factor$multiplier3 < 0, 0, factor$multiplier3)
factor$multiplier4 <- ifelse(factor$multiplier4 < 0, 0, factor$multiplier4)
factor$multiplier11 <- ifelse(factor$multiplier11 < 0, 0, factor$multiplier11)
plot(factor$temp, factor$multiplier)
plot(factor$temp, factor$multiplier1)
plot(factor$temp, factor$multiplier2)
plot(factor$temp, factor$multiplier3)
plot(factor$temp, factor$multiplier4)
plot(factor$temp, factor$multiplier11)
max(factor$multiplier)

gau <- function(x){
  if(x<= 28) exp(-((x-28)^2)/(2*7^2))*150
  else exp(-((x-28)^2)/(2*3^2))*150
}

gau(18)
factor <- data.frame(temp = seq(10,40,0.001))
factor$multiplier <- sapply(factor$temp, function(x) gau(x))
plot(factor$temp, factor$multiplier)


#### Variable size ####

setwd("~/Code/Parental-Effects")
num_samples <-1

par_num <- paste("parameter_values_", as.character(0), ".csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 1:(num_samples - 1)){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- #seq(1,10,1)
rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 0:(num_samples - 1)){
    dist_num <- paste("population_distribution_", as.character(n), ".csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    geno_0
    
    val_num <- paste("gene_values_", as.character(n), ".csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    y <- g_0[c(1,2,4)] #p
    w <- g_0[c(1,2,5)] #pe
    z <- g_0[c(1,2,6)] #e
    a_0 <- g_0[c(1,2,11)] #diff
    t <- g_0[c(1,2,7)] #threshold
    o <- g_0[c(1,2,8)] #overwinter
    ag <- g_0[c(1,2,9)] #age
    siz <- g_0[c(1,2,10)] # pop size
    
    x$type <- "gene"
    y$type <- "pheno"
    w$type <- "peff"
    z$type <- "env"
    a_0$type <- "diff"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    siz$type <- "popsize"
    g_0 <- rbind(x,y,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "gene")
    g_0$type <- relevel(g_0$type, "pheno")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,ag)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,1000000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    val_2
    NoPEvsPE_evolve_blocks <- plot_grid(geno_0, val_0, val_1, val_2, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n + 1]))
    choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n + 1])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}

#### variable size, long evo ####

setwd("~/Collect_data/Variable_size_long_evo")
num_samples <-4

par_num <- paste("parameter_values_", as.character(9), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 9:9){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    y <- g_0[c(1,2,4)] #p
    w <- g_0[c(1,2,5)] #pe
    z <- g_0[c(1,2,6)] #e
    a_0 <- g_0[c(1,2,11)] #diff
    t <- g_0[c(1,2,7)] #threshold
    o <- g_0[c(1,2,8)] #overwinter
    ag <- g_0[c(1,2,9)] #age
    siz <- g_0[c(1,2,10)] # pop size
    
    x$type <- "gene"
    y$type <- "pheno"
    w$type <- "peff"
    z$type <- "env"
    a_0$type <- "diff"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    siz$type <- "popsize"
    g_0 <- rbind(x,y,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "gene")
    g_0$type <- relevel(g_0$type, "pheno")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,ag)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,1000000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    val_2
    NoPEvsPE_evolve_blocks <- plot_grid(geno_0, val_0, val_1, val_2, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n + 1]))
    choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n + 1])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}

#### variable size, resistances ####

setwd("~/Collect_data/Variable_size_resistances")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n + 1]))
    choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n + 1])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### variable size, resistances, stable ####

setwd("~/Collect_data/Variable_size_resistances_stable")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### variable size, resistances, unstable ####

setwd("~/Collect_data/Variable_size_resistances_unstable")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}

#### variable size, resistances, v unstable ####

setwd("~/Collect_data/Variable_size_resistances_v_unstable")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### variable size, resistances, max ####

setwd("~/Collect_data/Variable_size_resistances_max")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n + 1]))
    choice <- round(as.numeric(as.character(p$iChoice[n + 1])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n + 1])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### variable size, resistances, stable, max ####

setwd("~/Collect_data/Variable_size_resistances_stable_max")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 1001
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### variable size, resistances, unstable, max ####

setwd("~/Collect_data/Variable_size_resistances_unstable_max")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}

#### variable size, resistances, v unstable, max ####

setwd("~/Collect_data/Variable_size_resistances_v_unstable_max")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### variable size, resistances, max, 50k ####

setwd("~/Collect_data/Variable_size_resistances_max_50k")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    #g_0$diff <- g_0[,6]-g_0[,4]
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #heat switch
    fs <- g_0[c(1,2,10)] #freeze switch
    siz <- g_0[c(1,2,11)] # pop size
    
    x$type <- "fpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "heatswitch"
    fs$type <- "freezeswitch"
    siz$type <- "popsize"
    g_0 <- rbind(x,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    g_0$type <- relevel(g_0$type, "fpoint")
    g_0$type <- relevel(g_0$type, "env")
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,hs,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,50000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### heatpoint, stable ####

setwd("~/Collect_data/heatpoint_stable")
num_samples <-8

par_num <- paste("parameter_values_", as.character(1), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:8){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    g_0$diff <- g_0$dHeatPoint-g_0$dTempDiff
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
  
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #temp diff
    fs <- g_0[c(1,2,10)] #porb slope
    siz <- g_0[c(1,2,11)] # pop size
    diff <- g_0[c(1,2,12)] #low point
    x$type <- "hpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "tempdiff"
    fs$type <- "probslope"
    siz$type <- "popsize"
    diff$type <- "lpoint"
    g_0 <- rbind(x,diff,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,1000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### sex, diploid, stable ####

setwd("~/Collect_data/sex_diploid_stable")
num_samples <-8

par_num <- paste("parameter_values_", as.character(2), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:num_samples){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    g_0$diff <- g_0$dHeatPoint-g_0$dTempDiff
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #temp diff
    fs <- g_0[c(1,2,10)] #porb slope
    fem <- g_0[c(1,2,11)] #% of females
    siz <- g_0[c(1,2,12)] # pop size
    diff <- g_0[c(1,2,13)] #low point
    x$type <- "hpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "tempdiff"
    fs$type <- "probslope"
    siz$type <- "popsize"
    diff$type <- "lpoint"
    g_0 <- rbind(x,diff,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,1000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### sex, diploid, evolve ####

setwd("~/Collect_data/sex_diploid_evolve")
num_samples <-4

par_num <- paste("parameter_values_", as.character(2), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:num_samples){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    g_0$diff <- g_0$dHeatPoint-g_0$dTempDiff
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #temp diff
    fs <- g_0[c(1,2,10)] #porb slope
    fem <- g_0[c(1,2,11)] #% of females
    siz <- g_0[c(1,2,12)] # pop size
    diff <- g_0[c(1,2,13)] #low point
    x$type <- "hpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "tempdiff"
    fs$type <- "probslope"
    siz$type <- "popsize"
    diff$type <- "lpoint"
    g_0 <- rbind(x,diff,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,fs)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,1000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}


#### sex, diploid, unstable ####

setwd("~/Collect_data/sex_diploid_unstable")
num_samples <-8

par_num <- paste("parameter_values_", as.character(2), "_0.csv", sep = "")
p <- read.csv(par_num, header = FALSE)
p <- as.data.frame(t(p))
colnames(p) <- as.character(unlist(p[1,]))
p <- p[-1,]
p$per <- 0

for (n in 2:num_samples){
  par_num <- paste("parameter_values_", as.character(n), "_0.csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}

rsamples <- seq(1,1,1)
#rsamples <- sampleDist(1, 1, 1)
intstart <- 1
intend <- 101
for(interv in 1:3){
  if(interv == 1){
    intstart <- 1
    intend <- 101
  }
  else if(interv == 2){
    intstart <- 1
    intend <- max(as.numeric(as.character(p$iGenMax)))
  }
  else if(interv == 3){
    intstart <- (1 + max(as.numeric(as.character(p$iGenMax)))) / 2
    intend <- intstart + 100
  }
  for(n in 1:num_samples){
    dist_num <- paste("population_distribution_", as.character(n), "_0.csv", sep = "")
    d_0 <- read.csv(dist_num, header = TRUE)
    
    colnames(d_0) <- c("replicate", "gen", "dist", "choice")
    
    d_0$PE <- ifelse(d_0$choice < 6, "No PE", "PE")
    d_0$choice <- as.factor(d_0$choice)
    d_0$replicate <- as.factor(d_0$replicate)
    d_0 <- d_0[d_0$replicate %in% rsamples, ]
    d_0 <- subset(d_0, gen <= intend & gen >=intstart)
    
    geno_0 <- ggplot(data = d_0, aes(x = gen, y = dist, fill = choice)) +
      geom_bar(stat = "identity", width = 1) + facet_grid(rows = vars(replicate)) +
      ylab("Distribution") + xlab("Generation") +
      scale_color_discrete(name = c("Genotype"),
                           breaks = c("0", "1","2","3","4","5",
                                      "6","7","8","9","10","11"),
                           labels = c("NoPE 0", "NoPE 1","NoPE 2",
                                      "NoPE 3","NoPE 4", "NoPE 5",
                                      "PE 0", "PE 1","PE 2",
                                      "PE 3","PE 4", "PE 5")) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()
      )
    #geno_0
    
    val_num <- paste("gene_values_", as.character(n), "_0.csv", sep = "")
    g_0 <- read.csv(val_num, header = TRUE)
    
    g_0$diff <- g_0$dHeatPoint-g_0$dTempDiff
    #g_0 <- subset(g_0, V2 < 13)
    
    colnames(g_0) <- c("replicate", "gen","value", "value", "value", "value","value","value","value","value","value","value","value")
    g_0$replicate <- as.factor(g_0$replicate)
    g_0 <- g_0[g_0$replicate %in% rsamples, ]
    
    
    x <- g_0[c(1,2,3)] #g
    w <- g_0[c(1,2,4)] #peweight
    z <- g_0[c(1,2,5)] #e
    t <- g_0[c(1,2,6)] #threshold
    o <- g_0[c(1,2,7)] #overwinter
    ag <- g_0[c(1,2,8)] #age
    hs <- g_0[c(1,2,9)] #temp diff
    fs <- g_0[c(1,2,10)] #porb slope
    fem <- g_0[c(1,2,11)] #% of females
    siz <- g_0[c(1,2,12)] # pop size
    diff <- g_0[c(1,2,13)] #low point
    x$type <- "hpoint"
    w$type <- "peweight"
    z$type <- "env"
    t$type <- "thres"
    o$type <- "overw"
    ag$type <- "age"
    hs$type <- "tempdiff"
    fem$type <- "fem"
    fem$value <- fem$value / siz$value
    fs$type <- "probslope"
    siz$type <- "popsize"
    diff$type <- "lpoint"
    g_0 <- rbind(x,diff,z,t)
    g_0 <- subset(g_0,  gen <= intend & gen >=intstart)
    
    g_0$type <- as.factor(g_0$type)
    
    
    val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_0
    
    g_1 <- rbind(w,o,fs, fem)
    g_1 <- subset(g_1,  gen <= intend & gen >=intstart)
    
    g_1$type <- as.factor(g_1$type)
    
    
    val_1 <- ggplot(data = g_1, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    
    #val_1
    
    g_2 <- rbind(siz)
    g_2 <- subset(g_2,  gen <= intend & gen >=intstart)
    
    g_2$type <- as.factor(g_2$type)
    
    
    val_2 <- ggplot(data = g_2, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      scale_y_continuous(trans = "log10", limits = c(1,1000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_2
    
    
    g_3 <- rbind(ag)
    g_3 <- subset(g_3,  gen <= intend & gen >=intstart)
    
    g_3$type <- as.factor(g_3$type)
    
    
    val_3 <- ggplot(data = g_3, aes(x = gen, y = value, colour = type)) +
      geom_line() + facet_grid(rows = vars(replicate)) +
      ylab("value") + xlab("Generation") +
      theme_minimal() +
      #scale_y_continuous(trans = "log10", limits = c(1,5000)) + 
      #scale_color_manual(values = c("#F8766D", "#619DFF", "#00BA38")) +
      theme(panel.spacing = unit(1, "lines"),
            strip.background = element_blank())
    #val_3
    
    NoPEvsPE_evolve_blocks <- plot_grid(val_0, val_1, val_2, val_3, ncol = 2)
    
    cost <- as.numeric(as.character(p$dCostOfAging[n]))
    choice <- round(as.numeric(as.character(p$iChoice[n])))
    #init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
    mut <- ifelse(as.numeric(as.character(p$dPEMutationRate[n])) == 0, "nomut", "mut")
    
    if(interv == 1) png_name <- paste("initial/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 2) png_name <- paste("total/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    else if(interv == 3) png_name <- paste("mid/",mut,"_choice",choice,"_agecost",cost,"_gen",intstart,"-",intend,"_",as.character(n),".png", sep = "")
    save_plot(png_name, NoPEvsPE_evolve_blocks,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 2, # and 2 rows
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.7,
              base_height = 6
    )
  }
}
