
## Alex Wiebe
## Princeton University
## Description: Code to create figures showing the species-area
## relationships (SARs) of regions with various parameter values.

setwd("")

library(ggplot2)
library(tidyverse)

# Define the SAR parameters
z = 0.15
c = 200
area = 1000
max_richness = c * area ^ z

# Create a data frame with values of A between 0 and 1
A_values = seq(1, 0, length.out = 1000)
S_values = (max_richness - c * (A_values * area) ^ z)
A_values = abs(A_values - 1)
data = data.frame(A = A_values, S = S_values)

A2 = seq(1, 0, length.out = 1000)

## Choose how a region will be different from the above baseline:
# # No change:
# S2 = (max_richness - c * (A2 * area) ^ z)

# # Larger area:
# S2 = (max_richness - c * (A2 * area) ^ z) * (10^z)

# # Higher z:
# z2 = 0.4
# c2 = max_richness / (area ^ z2)
# S2 = (max_richness - c2 * (A2 * area) ^ z2)

# Higher c:
c2 = 400
S2 = (max_richness * 2 - c2 * (A2 * area) ^ z)

A2 = abs(A2 - 1)
data2 = data.frame(A2 = A2, S2 = S2)

# Plot the function using ggplot2
plt = ggplot(data, aes(x = A, y = S)) +
  geom_line(color = "black", size = 3, linetype = "dashed") +
  geom_line(data = data2, aes(x = A2, y = S2), size = 3) +
  theme_classic()+
  labs(x='Proportional Habitat Loss', y = "Number of Extinctions") +
  theme(
    axis.title = element_text(size = 36),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30)
  )
plt  

# To save:
filename = paste0("one_region_EAR_higherc.png")
png(filename, units = "px", height=3500,
    width=5000, res=600)
plt
dev.off()
























