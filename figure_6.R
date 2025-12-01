# title: figure_6
# author: Sarah Blecksmith
# purpose: Make Figure 6 -- scatterplots of transformed muc2plant and energy adjusted habitual fiber by secretor status

library(tidyverse)
library(ggplot2)
library(bestNormalize)

FL100_data <- read.csv("data/FL100_merged_variables.csv", header = TRUE) 

muc2plant <- read.csv("data/muc2plant_ratioGH_GHPL.csv", header = TRUE) %>%
  select(c(subject_id, muc2plantGH, muc2plantGHPL))

fut2 <- read.csv("data/forR_fut2.csv", header = TRUE)

merged <- merge(x=fut2, y=muc2plant, by = 'subject_id') %>%
  merge(., y=FL100_data[,c("subject_id", "fibe_perKcal", "sol_fibe_perKcal")], by = "subject_id") 

merged$secr_status <- as.factor(merged$secr_status)

# Transform muc2plant
merged$muc2plantGHPL_log <- log_x(merged$muc2plantGHPL)$x.t

secretor_lm <- lm(data= merged, muc2plantGHPL_log ~ fibe_perKcal + secr_status + fibe_perKcal*secr_status)
summary(secretor_lm)


fig6 <- ggplot(merged, aes(x = fibe_perKcal, y = muc2plantGHPL_log)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "muc2plant (log-transformed)") +
  facet_wrap(vars(secr_status))

fig6
ggsave("figure_6.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", fig6)