
### Alex Wiebe
## Princeton University
## Description: Code for analysis of empirical extinction
## data in birds.

setwd("")
library(ggplot2)
library(tidyverse)
library(ggh4x)

bird = read.csv("Wiebe_extinct_bird_data.csv", header = T,
                stringsAsFactors = T)

bird$Decade = as.integer(cut(bird$Year, breaks = seq(1500, 2030, by = 10), labels = seq(1500, 2020, by = 10)))
bird = bird[bird$Year < 2020,] # Remove data from 2020's for this first part because we do not 
# have a complete dataset from that decade yet

tabulate_extinctions = function(region, cause){
  if(region == "all"){ # All countries
    species = bird
  }
  if(region == "first world"){ # "First world" countries/developed countries
    species = bird[bird$Firstworld_atthetime != "",]
  }
  if(region == "third world"){ # "Third world" countries/developing countries
    species = bird[bird$Firstworld_atthetime == "",]
  }
  if(cause == "LUC"){
    species = species[species$LUC_main_cause == "x" | species$LUC_partial_cause == "x",]
  }
  if(cause == "Primarily LUC"){
    species = species[species$LUC_main_cause == "x",]
  }
  
  tabulated_species = table(species$Decade)
  tabulated_species = data.frame(Decade = as.numeric(names(tabulated_species)),
                                 Count = as.numeric(tabulated_species))
  decades = data.frame(Decade = seq(6,52, by = 1),
                       Count = rep.int(0,47))
  tabulated_species = merge(tabulated_species, decades,
                            by = "Decade", all = T)
  tabulated_species[is.na(tabulated_species)] = 0
  tabulated_species$Count = tabulated_species$Count.x + tabulated_species$Count.y
  tabulated_species$Decade = 1490 + tabulated_species$Decade * 10
  
  tabulated_species
}


luc = bird[bird$LUC_main_cause == "x" | bird$LUC_partial_cause == "x",]
# Number of species in the Extinct category (79)
luc %>%
  filter(is.na(Nonextinct_status) | Nonextinct_status == "") %>%
  count()
# Number of species in the Extinct in the Wild category (3)
luc %>%
  filter(Nonextinct_status == "EW") %>%
  count()
# Number of species in the CR PE category (22)
luc %>%
  filter(Nonextinct_status == "CR - Possibly Extinct") %>%
  count()
# Number of species in the CR category (28: 27 + one species in 2024)
luc %>%
  filter(Nonextinct_status == "CR") %>%
  count()
# Only 10 species that went extinct due primarily or partially to LUC went
# extinct before 1800
luc %>%
  filter(Year < 1800) %>%
  count()

# Figure 3c
species = tabulate_extinctions(region = "first world", cause = "LUC")
new.species = tabulate_extinctions(region = "third world", cause = "LUC")
all.species = tabulate_extinctions(region = "all", cause = "LUC")
species = species[species$Decade > 1790,]
new.species = new.species[new.species$Decade > 1790,]
all.species = tabulate_extinctions("all", "LUC")
all.species = all.species[all.species$Decade > 1790,]
plt = ggplot(species, aes(x = Decade, y = Count)) +
  geom_point(stat = "identity", color = "darkred") +
  geom_smooth(method = "loess", se = T, color = "darkred", linewidth = 1.5, span = 1,
              fill = "lightpink") +  
  labs(
       x = "Decade",
       y = "Extinctions") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24)
  ) +
  coord_cartesian(ylim = c(0,11)) +
  geom_point(data = new.species, color = "darkgreen") +
  geom_smooth(data = new.species, method = "loess", se = T, color = "darkgreen",
              fill = "lightgreen", linewidth = 1.5, span = 1) +
  geom_smooth(data = species, method = "loess", se = T, color = "darkred",
              fill = NA, linewidth = 2, span = 1) +
  geom_smooth(data = new.species, method = "loess", se = T, color = "darkgreen",
              fill = NA, linewidth = 2, span = 1) +
  geom_smooth(data = all.species, method = "loess", se = T, color = "darkblue",
              fill = NA, linewidth = 2, span = 1, linetype = "dashed") +
  geom_point(data = new.species, color = "darkgreen") +
  geom_point(stat = "identity", color = "darkred")
plt
# To save:
png(paste0("extinctions_by_bloc.png"), units = "px", height=3500,
    width=5000, res=600)
plt
dev.off()

# Hypothesis testing

# Developed bloc
devd.lm = lm(Count ~ Decade, data = species)
summary(devd.lm) # no significant coefficients
hist(resid(devd.lm)) # bad residuals
devd.quad = lm(Count ~ poly(Decade, 2), data = species)
summary(devd.quad) # significant coefficient on the squared term
hist(resid(devd.quad)) # nice and normal
AIC(devd.lm, devd.quad) # double check with AIC, smaller AIC for quadratic
# model, delta AIC > 2

# Developing bloc
devg.lm = lm(Count ~ Decade, data = new.species)
summary(devg.lm)
hist(resid(devg.lm)) # residuals look ok


### Figure 3d: To plot bird extinctions per year since the 1980s
recent = read.csv("Wiebe_extinct_bird_data.csv", header = T,
                  stringsAsFactors = T)
recent = recent[recent$Year >= 1984,]

recent = recent[recent$LUC_partial_cause == "x" | recent$LUC_main_cause == "x",]

recent_summary = recent %>%
  count(Year) %>%
  rename(Extinctions = n) %>%
  complete(Year = 1984:2024, fill = list(Extinctions = 0))

plt = ggplot(recent_summary, aes(x = Year, y = Extinctions)) +
  geom_bar(stat = "identity", col = "darkgray") +
  geom_smooth(data = recent_summary, method = "lm", se = T, color = "navy", fill = "#c3c6f7",
              linewidth = 2) +
  labs(
       x = "Year",
       y = "Extinctions") +
  coord_cartesian(ylim = c(0, 6)) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 30),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24)
  )
plt
summary(lm(Extinctions ~ Year, data = recent_summary))
summary(lm(sqrt(Extinctions) ~ Year, data = recent_summary))
hist(resid(lm(Extinctions ~ Year, data = recent_summary)))

# To save:
  png(paste0("recent_bird_extinctions.png"), units = "px", height=3500,
      width=5000, res=600)
  plt
  dev.off()







