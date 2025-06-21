

### Alex Wiebe
## Princeton University
## Description: Code to create figures showing the extinction-time
## relationships in Figure S1.

# Set the working directory
setwd("")

# Load required libraries
library(ggpattern)
library(tidyverse)

# Source additional functions
source("biodiversity_trade_functions.R")

### Run for a single country

# Identical developmental parameters for each country
params = c(s = 1/30, d = 0.001, a = 1/20, h = 1,
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


### Plots

plt = single %>%
  gather(variable, value, c(A, F, U, P)) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(aes(linetype = variable), size = 3)+
  scale_color_manual(values = c("red", "darkgreen", "black", "blue")) +
  scale_linetype_manual(values = c("A" = "solid", "F" = "solid",
                                   "U" = "solid", "P" = "dashed")) +
  theme_classic()+
  geom_line(data = single %>% gather(variable, value, c(F)) %>% filter(variable == "F"),
            aes(x = time, y = value, color = variable), 
            size = 3) +
  geom_line(data = single %>% gather(variable, value, c(A)) %>% filter(variable == "A"),
            aes(x = time, y = value, color = variable), 
            size = 3) +
  geom_line(data = single %>% gather(variable, value, c(P)) %>% filter(variable == "P"),
            aes(x = time, y = value, color = variable), 
            size = 3, linetype = "dashed") +
  xlim(0, 300) +
  labs(x='Time (year)', y = expression(Area ~ (km^2))) +
  theme(
    axis.title = element_text(size = 36),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    legend.position = "none"
  )
plt

# To save:
filename = paste0("one_region_dynamics.png")
png(filename, units = "px", height=3500,
    width=5000, res=600)
plt
dev.off()


### Calculating species extinctions

c = 200
z = 0.15

# # Set 1
# c2 = c
# z2 = z
# S = c * single$F ^ z
# S2 = c2 * single$F ^ z2
# single$S = S
# single$S2 = S2

# # Set 2
# z2 = z
# c2 = c
# S = c * single$F ^ z
# S2 = c2 * (single$F * 10) ^ z2
# single$S = S
# single$S2 = S2

# # Set 3
# z2 = 0.4
# c2 = (c * 1000 ^ z) / (1000 ^ z2)
# S = c * single$F ^ z
# S2 = c2 * single$F ^ z2
# single$S = S
# single$S2 = S2

# Set 4
z2 = z
c2 = 400
S = c * single$F ^ z
S2 = c2 * single$F ^ z2
single$S = S
single$S2 = S2


# Add columns for extinction rate changes
single = single %>% add_column(dS = NA) %>%
  add_column(dS2 = NA)


# Calculate the changes in extinction rates
for (i in 1:length(single$time)){
  if (i == 1){
    
    single$dS[i] = NA
    single$dS2[i] = NA
    
  } else {
    
    single$dS[i] = single$S[i-1] - single$S[i]
    single$dS2[i] = single$S2[i-1] - single$S2[i]
    
  }
}

plt = 
  ggplot(single, aes(x = time, y = dS)) +
  geom_line(col = "maroon", size = 3, linetype = "dashed") +
  geom_line(data = single, aes(x = time, y = dS2), size = 3, col = "maroon") +
  theme_classic()+
  labs(x='Time (year)', y = 'Number of Extinctions') +
  ylim(0, 5) + 
  xlim(0, 300) + 
  theme(
    axis.title = element_text(size = 36),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    legend.position = "none"
  )
plt

# To save:
filename = paste0("one_region_extinctions_higherc.png")
png(filename, units = "px", height=3500,
    width=5000, res=600)
plt
dev.off()








