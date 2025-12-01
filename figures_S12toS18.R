# title: figures_S12toS18
# author: Sarah Blecksmith
# purpose: Make supplmental figures S12 through S18, SCFAs and energy adjusted habitual fiber by secretor status and deviations 
# from expected SCFA ratios and energy adjusted habitual fiber by secretor status

library(tidyverse)
library(ggplot2)
library(bestNormalize)

FL100_data <- read.csv("data/FL100_merged_variables.csv", header = TRUE) 

muc2plant <- read.csv("data/muc2plant_ratioGH_GHPL.csv", header = TRUE) %>%
  select(c(subject_id, muc2plantGH, muc2plantGHPL))

fut2 <- read.csv("data/forR_fut2.csv", header = TRUE)

scfa <- read.csv("data/fecal_SCFA_FL100.csv", header = TRUE) %>%
  rename(subject_id = SubjectID)

merged <- merge(x=fut2, y=muc2plant, by = 'subject_id') %>%
  merge(., y=FL100_data[,c("subject_id", "fibe_perKcal", "sol_fibe_perKcal")], by = "subject_id") %>%
  merge(., y=scfa, by = "subject_id")

merged$secr_status <- as.factor(merged$secr_status)


#Figures S12-S15 - SCFAs and energy adjusted habitual fiber by secretor status

# Transform the SCFA values
bestNormalize(merged$totalSCFA)
merged$totalSCFA_boxcox <- boxcox(merged$totalSCFA)$x.t
bestNormalize(merged$acetate)
merged$acetate_boxcox <- boxcox(merged$acetate)$x.t
bestNormalize(merged$propionate)
merged$propionate_yj <- yeojohnson(merged$propionate)$x.t
bestNormalize(merged$butyrate)
merged$butyrate_boxcox <- boxcox(merged$butyrate)$x.t

# Linear models
totalSCFA_lm <- lm(data= merged, totalSCFA_boxcox ~ fibe_perKcal + secr_status + fibe_perKcal*secr_status)
summary(totalSCFA_lm)
acetate_lm <- lm(data= merged, acetate_boxcox ~ fibe_perKcal + secr_status + fibe_perKcal*secr_status)
summary(acetate_lm)
propionate_lm <- lm(data= merged, propionate_yj ~ fibe_perKcal + secr_status + fibe_perKcal*secr_status)
summary(propionate_lm)
butyrate_lm <- lm(data= merged, butyrate_boxcox ~ fibe_perKcal + secr_status + fibe_perKcal*secr_status)
summary(butyrate_lm)


figS12 <- ggplot(merged, aes(x = fibe_perKcal, y = totalSCFA_boxcox)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "total SCFA (Box Cox-transformed)") +
  facet_wrap(vars(secr_status))
figS12
ggsave("figure_S12.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", figS12)

figS13 <- ggplot(merged, aes(x = fibe_perKcal, y = acetate_boxcox)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "acetate (Box Cox-transformed)") +
  facet_wrap(vars(secr_status))
ggsave("figure_S13.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", figS13)

figS13

figS14 <- ggplot(merged, aes(x = fibe_perKcal, y = propionate_yj)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "propionate (Yeo-Johnson-transformed)") +
  facet_wrap(vars(secr_status))
ggsave("figure_S14.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", figS14)

figS14

figS15 <- ggplot(merged, aes(x = fibe_perKcal, y = butyrate_boxcox)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "butyrate (Box Cox-transformed)") +
  facet_wrap(vars(secr_status))
figS15
ggsave("figure_S15.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", figS15)


# Fig S16-S19 - deviations from expected SCFA ratios and energy adjusted habitual fiber by secretor status

# relative abundance of each SCFA
merged$propionate_ratio <- (merged$propionate/merged$totalSCFA)
merged$butyrate_ratio <- (merged$butyrate/merged$totalSCFA)
merged$acetate_ratio <- (merged$acetate/merged$totalSCFA)

## dist to norm (60:20:20)
merged$acetate_ratio_dist <- merged$acetate_ratio - 0.6
merged$butyrate_ratio_dist <- merged$butyrate_ratio - 0.2
merged$propionate_ratio_dist <- merged$propionate_ratio - 0.2

# Transform
bestNormalize(merged$acetate_ratio_dist)
merged$acetate_ratio_dist_yj <- yeojohnson(merged$acetate_ratio_dist)$x.t
bestNormalize(merged$propionate_ratio_dist)
merged$propionate_ratio_dist_yj <- yeojohnson(merged$propionate_ratio_dist)$x.t
bestNormalize(merged$butyrate_ratio_dist)
merged$butyrate_ratio_dist_log <- log_x(merged$butyrate_ratio_dist)$x.t

figS16 <- ggplot(merged, aes(x = fibe_perKcal, y = acetate_ratio_dist_yj)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "Deviation from acetate ratio (Yeo-Johnson-transformed)") +
  facet_wrap(vars(secr_status))
figS16
ggsave("figure_S16.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", figS16)

figS17 <- ggplot(merged, aes(x = fibe_perKcal, y = propionate_ratio_dist_yj)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "Deviation from propionate ratio (Yeo-Johnson-transformed)") +
  facet_wrap(vars(secr_status))
figS17
ggsave("figure_S17.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", figS17)

figS18 <- ggplot(merged, aes(x = fibe_perKcal, y = butyrate_ratio_dist_log)) +
  geom_point(aes(colour = "#464646"))+
  scale_color_identity() +
  geom_smooth(method = "lm", formula = y~x, se = T, color = "blue") +
  labs(
    x = "habitual fiber intake (grams/1000 kcal)",
    y = "Deviation from butyrate ratio (log-transformed)") +
  facet_wrap(vars(secr_status))
figS18
ggsave("figure_S18.tiff", device = "tiff", dpi = 300, width = 8, height = 4, units = "in", path = "output", figS18)
