

### Alex Wiebe
## Princeton University
## Description: Code to create data for figures showing the relative benefits
## of small-scale conservation.

# Set the working directory
setwd("")

# Load required libraries
library(ggpattern)
library(tidyverse)

# Source additional functions
source("biodiversity_trade_functions.R")

### Run for a single country

# Identical developmental parameters for each country
params = c(s = 1/20, d = 0.001, a = 1/20, h = 1,
           r = 0.04)

# Time sequence from 0 to 300 in steps of 1
times = seq(from = 0, to = 500, by = 1)

# Initial conditions for Country 1
xstart = c(F = 999, A = 1, U = 0, P = 1)

# Remove any existing variable named 'single'
rm(single)

# Run the ODE for the single country model
ode(
  func = single.model,
  y = xstart,
  times = times,
  parms = params
) %>%
  as.data.frame() -> single

### Run with both countries

# Identical developmental parameters for each country, including migration and trade
params = c(s1 = 1/20, s2 = 1/20, d1 = 0.001, d2 = 0.001, a1 = 1/20, a2 = 1/20, h = 1,
           r = 0.04, m1 = 0.001, m2 = 0.001, e = 0.5)

# Time sequence from 0 to 500 in steps of 1
times = seq(from = 0, to = 500, by = 1)

# Delay before the second country starts development
delay = 0

# Initial conditions for both countries, with second country starting later
xstart = c(F1 = single[delay + 1,]$F, A11 = single[delay + 1,]$A, A12 = 0,
           U1 = single[delay + 1,]$U, P1 = single[delay + 1,]$P,
           F2 = 999, A22 = 1, A21 = 0, U2 = 0, P2 = 1)

### Run the ODEs for two countries

# Remove any existing variable named 'two'
rm(two)

# Run the ODE for the two country model
ode(
  func = faup.model,
  y = xstart,
  times = times,
  parms = params
) %>%
  as.data.frame() -> two

# Remove the first row of the result
two = two %>%
  slice(-1)

### Combine the datasets for a plot
full.ag = data.frame(time = rep(0:c(length(c(single[single$time <= delay,]$A, two$A11)) - 1)),
                     A11 = c(single[single$time <= delay,]$A, two$A11),
                     A12 = c(rep(0, delay + 1), two$A12),
                     A22 = c(rep(0, delay + 1), two$A22),
                     A21 = c(rep(0, delay + 1), two$A21),
                     F1 = c(single[single$time <= delay,]$F, two$F1),
                     F2 = c(rep(xstart["F2"], delay + 1), two$F2),
                     P1 = c(single[single$time <= delay,]$P, two$P1),
                     P2 = c(rep(xstart["P2"], delay + 1), two$P2))

### Separate section for calculating species diversity (with any "full.ag" data.frame)
# Define constants for species extinction rate calculation
c1 = 100
z1 = 0.25
# # Beta = 1:
# z2 = z1
# c2 = c1
# # Beta = 0.5
# z2 = 0.3
# c2 = 2*c1*1000^(z1-z2)
# Beta = 0
z2 = 0.35
# # c2 = 5*c1*1000^(z1-z2)
c2 = 250

# # If testing the impact of a different z value:
# z2 = 0.35
# c2 = c1 * (max(full.ag$F1) + 1) ^ (z1 - z2)

# Calculate the species extinction rates
S1 = c1 * full.ag$F1 ^ z1
S2 = c2 * full.ag$F2 ^ z2
full.ag$S1 = S1
full.ag$S2 = S2

# Add columns for extinction rate changes
full.ag = full.ag %>% add_column(dS1 = NA) %>%
  add_column(dS2 = NA) %>%
  add_column(dS_global = NA)

# Calculate the changes in extinction rates
for (i in 1:length(full.ag$time)){
  if (i == 1){
    
    full.ag$dS1[i] = NA
    full.ag$dS2[i] = NA
    full.ag$dS_global[i] = NA

    
  } else {
    
    full.ag$dS1[i] = full.ag$S1[i-1] - full.ag$S1[i]
    full.ag$dS2[i] = full.ag$S2[i-1] - full.ag$S2[i]
    
    if(full.ag$dS1[i] < 0){
      full.ag$dS1[i] = 0
    }
    full.ag$dS_global[i] = sum(full.ag$dS1[i], full.ag$dS2[i], na.rm = T)

  }
}

# Compute the marginal species loss of each region
full.ag$S1_marginal = c1 * z1 * full.ag$F1^(z1 -1)
full.ag$S2_marginal = c2 * z2 * full.ag$F2^(z2 -1)
full.ag$best_conservation = full.ag$S1_marginal/full.ag$S2_marginal

# b1 = full.ag$best_conservation
# b2 = full.ag$best_conservation
# b3 = full.ag$best_conservation
# betas = data.frame(beta1 = b1, beta05 = b2, beta0 = b3)
# # write.csv(betas, "tworegion_relativeconservation.csv", row.names = F)
# write.csv(betas, "tworegion_simultaneousconservation.csv", row.names = F)

# betas = read.csv("tworegion_relativeconservation.csv", header = T)
betas = read.csv("tworegion_simultaneousconservation.csv", header = T)

full.ag = cbind(full.ag, betas)












