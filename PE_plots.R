library(ggplot2)
library(cowplot)

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

val_0 <- ggplot(data = g_0, aes(x = gen, y = value, colour = type)) +
  geom_line() + facet_grid(rows = vars(replicate)) +
  ylab("value") + xlab("Generation") +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text = element_blank())

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

for (n in 1:5){
  par_num <- paste("parameter_values_", as.character(n), ".csv", sep = "")
  p_temp <- read.csv(par_num, header = FALSE)
  p_temp <- as.data.frame(t(p_temp))
  colnames(p_temp) <- as.character(unlist(p_temp[1,]))
  p_temp <- p_temp[-1,]
  p_temp$per <- n
  p <- rbind(p, p_temp)
}


for(n in 0:5){
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
  
  block <- round(as.numeric(as.character(p$iInterval[n + 1])))
  init <- ifelse(as.numeric(as.character(p$bPEStart[n + 1])) == 1, "PE", "NoPE")
  
  png_name <- paste("NoPE_block", as.character(block),"_",init, ".png", sep = "")
  
  save_plot(png_name, NoPEvsPE_blocks,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2, # and 2 rows
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3,
            base_height = 6
  )
}
