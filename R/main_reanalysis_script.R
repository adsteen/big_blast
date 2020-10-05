############
# Script to find correct cutoff for "short reads" prior to reanalyzing
############

# Load packages
library(dplyr)
library(tibble)
library(ggplot2)
library(readr)

# Read the file for archaeal isolates, as a test
arc_is <- read_csv("/srv/data/big_blast/results/JGI_arc_isolates_top_bitscores.csv")

# Example plot of the Archaeal isolates
p_arc_is <- ggplot(arc_is, aes(x=alignment_length)) + 
  geom_density() + 
  geom_vline(xintercept = 170) + 
  scale_x_log10()
print(p_arc_is)

# Read in the rest of the files
arc_mags <- read_csv("/srv/data/big_blast/results/JGI_arc_mags_top_bitscores.csv")

# Add columns for the source of data
arc_is <- arc_is %>% mutate(source = "arc_isolates")
arc_mags <- arc_mags %>% mutate(source = "arc_mags")

# Combine all the data sets into one big data frame
all_data <- arc_is %>%
  rbind(arc_mags) # rbind binds columns together by row
# Equivalent to writing all_data <- rbind(arc_is, bac_is)

# Plot of the combined data
p_dens_combined <- ggplot(all_data, aes(x=alignment_length, color=source)) + 
  geom_density() +
  geom_vline(xintercept = 170) +
  scale_x_log10()
print(p_dens_combined)
