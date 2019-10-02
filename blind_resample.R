# script to randomize a list of samples to create blind experiments
# input is an xlsx file with two columns and output is an xlsx file where the second column is randomly resampled
# Mário Santos Mira - 02/10/2019

# seed value to be used
seed <- 1234

# name of input file
input <- "samples_input.xlsx"

# necessary in order to read xlsx files (this should also work with csv by replacing the respective functions and file names)
library(xlsx)

# read input list
input_list <- read.xlsx(input, sheetIndex = 1)

# retrive sample size
sample_size <- nrow(input_list)

# create output data frame
output_list <- input_list

# set a seed to resapmling is replicable and randomly resample the second column
set.seed(seed)
output_list$sample_number <- sample(input_list$sample_number, size = sample_size, replace = F)

# write the resampled output to file
write.xlsx(output_list, "samples_output.xlsx", row.names = F)
