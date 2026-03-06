## Clean code to support data analysis in Christensen et al.
## for datawrangling from raw files performed elsewhere (Mary's script)
## begin here with cleaned data (DF3).

## packages
library(tidyverse)
library(nlme)
library(MuMIn)
library(car)
library(ggplot2)
library(effects)


#Read in data ----
DF3 <- read.csv(file = "./data/BCdataforanalysis.csv")
DF_LOI <- read.csv(file = "./data/BlueCarbonData_1.csv")
core_info <- read.csv(file = "./data/core_info.csv")

#set formulae to use later
V_cseg <- pi*(2.38125^2) *1 # volume of 1 cm of core sediment, cm3
A_cseg <- pi*(2.38125^2) # area of top of core, cm2

# Modeling Carbon Density: Manuscript Methods/Statistical Analyses -------------------------------------------------

# Testing random effects
mod0.1 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ 1 | CoreName_2, method = "ML")

mod0.3 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ 1 | Site/CoreName_2, method = "ML")

mod0.0 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ corr_segment_midpoint | CoreName_2, method = "ML")

mod0.2 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ corr_segment_midpoint | Site/CoreName_2, method = "ML")

anova(mod0.2, mod0.3, mod0.0, mod0.1)
## mod0.2 is best

## comparing fixed effects
#mod0a <- lme(log(c_dens) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod0 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

#mod1a <- lme(log(c_dens) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod1 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

## all site level predictors
mod2 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + log(Watercourse_NEAR_DIST.x) + Coast, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod3 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod4 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + log(Watercourse_NEAR_DIST.x)*Type + Coast, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod5 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + log(Watercourse_NEAR_DIST.x), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod22 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + log(Watercourse_NEAR_DIST.x) + Type, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod6 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod7 <- lme(log(c_dens) ~1, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

## all site + core level predictors

mod8 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast + (corr_segment_midpoint), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod9 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod11 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod12 <- lme(log(c_dens) ~1 + Type + log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod13 <- lme(log(c_dens) ~1 + Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod14 <- lme(log(c_dens) ~1 + Type + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod15 <- lme(log(c_dens) ~1 + Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod16 <- lme(log(c_dens) ~1 + Type*log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod17 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod18 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod18a <- lme(log(c_dens) ~1 + sqrt(Percent.Silt.Fraction) + corr_segment_midpoint, data = DF3, random = ~corr_segment_midpoint| Site/CoreName_2, method = "REML")

mod19 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + Type + log(Watercourse_NEAR_DIST.x)  + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod20 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled) + Type + log(Watercourse_NEAR_DIST.x)  + corr_segment_midpoint, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod21 <- lme(log(c_dens) ~1 + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod23 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod24 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod25 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")


model.sel(mod0, mod1, mod2,mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod18a, mod19, mod20, mod21, mod22, mod23, mod24, mod25)

# Table 2 results ---------------------------------------------------------
summary(mod18)
summary(mod11)

## comparing predicted vs observed values
plot(predict(mod18), log(DF3$c_dens))
plot(predict(mod15), log(DF3$c_dens))
plot(predict(mod12), log(DF3$c_dens))
## best with readily available predictors
plot(predict(mod18a), log(DF3$c_dens)) # pretty similar
plot(exp(predict(mod18a)), DF3$c_dens) # not too bad


# plotting residuals ------------------------------------------------------

plot(Effect("REI_Scaled", mod18, residuals = TRUE))
ggsave("REI_resids_C_dens.pdf", path = "./figures/", width = 5, height = 3)

plot(Effect("Watercourse_NEAR_DIST.x", mod18, residuals = TRUE))
ggsave("Watercourse_C_dens.pdf", path = "./figures/", width = 5, height = 3)

plot(Effect("Percent.Silt.Fraction", mod18, residuals = TRUE))
ggsave("Mud_C_dens.pdf", path = "./figures/", width = 5, height = 3)

plot(Effect("Type", mod18, residuals = TRUE))
ggsave("SG_C_dens.pdf", path = "./figures/", width = 5, height = 3)


# Modeling mud and silt content [could cut from this file] -------------------------------------------

> mod18s <- lme(sqrt(Percent.Silt.Fraction) ~1 + sqrt(REI_Scaled) + Type + log(Watercourse_NEAR_DIST.x), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")
> mod1s <- lme(sqrt(Percent.Silt.Fraction) ~1 + sqrt(REI_Scaled) + log(Watercourse_NEAR_DIST.x), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")
> mod2s <- lme(sqrt(Percent.Silt.Fraction) ~1 + sqrt(REI_Scaled) + Type, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")
> mod3s <- lme(sqrt(Percent.Silt.Fraction) ~1 + Type, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")
> mod4s <- lme(sqrt(Percent.Silt.Fraction) ~1, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")
> model.sel(mod18s, mod1s, mod2s, mod3s, mod4s)



# Predicting C_stocks from C_dens models -------------------------------------------------------

#C_density _ model 18

## create the new data frame with the depths we want to predict C_dens and eventually C_stock for: 1 - 100 cm
new_data4 <- data.frame(CoreName_2 = rep(unique(DF3$CoreName_2), each = 100), corr_segment_midpoint = rep(c(1:100), times = length(unique(DF3$CoreName_2)))) 
                      
#merge new_data with DF3 to get site level variables: Site, REI_Raw, Type, Watercourse_NEAR_DIST.x
#how to get Percent.Silt.Fraction? this is a segment level variable. use mean silt val for core. 

Site_info <- DF3 %>% 
  group_by(Site, CoreName_2, REI_Scaled, Type, Watercourse_NEAR_DIST.x, Coast) %>%
  summarise(Percent.Silt.Fraction = mean(Percent.Silt.Fraction)) 

#use site level variables from above
new_data5 <- new_data4 %>%
  left_join(Site_info) %>%
  mutate(corr_segment_midpoint = as.numeric(corr_segment_midpoint))

predicted_vals <- predict(mod18, new_data5, level = 0:2)

pred_vals <- predicted_vals %>%
  rename(log_c_dens_cm3 = predict.CoreName_2) %>%
  separate(CoreName_2, into = c("Site", "CoreName_2"),
           sep ="/") %>%
  mutate(corr_segment_midpoint = rep(c(1:100), times = length(unique(DF3$CoreName_2)))) %>%
  mutate(c_dens = exp(log_c_dens_cm3)) %>%
  mutate(c_stock_cm = c_dens * V_cseg)  # c_stock in g per cm core length
 # mutate(c_stock_bythickness = c_dens*). -> then summarize for each core

new_data6 <- new_data5 %>%
  left_join(pred_vals) %>% 
  group_by(Site, CoreName_2, REI_Scaled, Type, Watercourse_NEAR_DIST.x) %>%
  mutate(c_60 = ifelse(corr_segment_midpoint < 61, c_stock_cm, 0)) %>%
  mutate(c_30 = ifelse(corr_segment_midpoint < 31, c_stock_cm, 0)) %>%
  summarise_at(., c("c_stock_cm", "c_60", "c_30"), sum) %>% #C g / inch core
  rename(c_100 = c_stock_cm) %>%
  mutate(c_60_m2 = 10000*(c_60/A_cseg)) %>% # c_stock / cm2
  mutate(c_100_m2 = 10000*(c_100/A_cseg)) %>%
  mutate(c_30_m2 = 10000*(c_30/A_cseg)) 

write.csv(new_data6, file = "./data/predicted_mod18.csv")

TableS3 <- new_data5 %>%
  left_join(pred_vals) %>% 
  group_by(Coast, Site, CoreName_2, Type, Watercourse_NEAR_DIST.x, Percent.Silt.Fraction) %>%
  mutate(c_30 = ifelse(corr_segment_midpoint < 31, c_stock_cm, 0)) %>%
  mutate(c_60 = ifelse(corr_segment_midpoint < 61, c_stock_cm, 0)) %>%
  summarise_at(., c("c_stock_cm", "c_60", "c_30"), sum) %>% #C g / inch core
  rename(c_100 = c_stock_cm) %>%
  mutate(c_30_m2 = 10000*(c_30/A_cseg)) %>%
  mutate(c_60_m2 = 10000*(c_60/A_cseg)) %>% # c_stock / cm2
  mutate(c_100_m2 = 10000*(c_100/A_cseg))

View(TableS3)
write.csv(TableS3, file = "./data/TableS3.csv"). # this is some information for Table S3. The site information comes from core_info, and the. max depth comes from the data. 

## checking
data7 <- DF3 %>%
  group_by(CoreName_2) %>%
  mutate(c_stock_est = 10000*(c_dens*V_cseg)/A_cseg) %>%
  summarise(Mean_c = mean(c_stock_est)) %>%
  mutate(ln_mnC = log(Mean_c)) %>%
  left_join(new_data6)

hist(data7$c_100) # g carbon / 100 cm core of diameter 17 cm2, or
hist(data7$c_100_m2) 
hist(data7$c_60_m2) 
hist(data7$c_30_m2)
  
## march 2026: these don't have the same theme as in the paper, not sure why. data looks the same. i will fix the theme here.

checking_100 <- ggplot(data7, aes(x = (Mean_c*100), y = c_100_m2)) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) +
  xlab("C Stock est by Mean C (gC / m^2)") +
  xlim(0, 60000) +
  ylim(0, 60000) +
  ylab("C Stock est by model (gC / m^2)") +
  ggtitle("Carbon Stock to 100 cm") +
  geom_abline(intercept = 0, slope = 1)
checking_100
ggsave("Pred vs mean C_Stock 100cm.tiff", path = "./figures/", width = 4, height = 4)

checking_60 <- ggplot(data7, aes(x = (Mean_c*60), y = c_60_m2)) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) +
  xlab("C Stock est by Mean C (gC / m^2)") +
  xlim(0, 35000) +
  ylim(0, 35000) +
  ylab("C Stock est by model (gC / m^2)") +
  ggtitle("Carbon Stock to 60 cm") +
  geom_abline(intercept = 0, slope = 1)
checking_60
ggsave("Pred vs mean C stock 60cm.tiff", path = "./figures/", width = 4, height = 4)

checking_30 <- ggplot(data7, aes(x = (Mean_c*30), y = c_30_m2)) +
    geom_point() +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) +
    xlab("C Stock est by Mean C (gC / m^2)") +
    xlim(0, 18000) +
    ylim(0, 18000) +
    ylab("C Stock est by model (gC / m^2)") +
      ggtitle("Carbon Stock to 30 cm") +
      geom_abline(intercept = 0, slope = 1)
checking_30
ggsave("Pred vs mean C stock 30cm.tiff", path = "./figures/", width = 4, height = 4)

#data7s <- subset(data7, Type == "SG") #i'm not sure why this is back in cm2
#(range(data7s$c_30_cm2))*10000
#(mean(data7s$c_30_cm2))*10000
#(mean(data7s$c_30_cm2)) # g / cm2
#(range(data7s$c_100_cm2))*10000
#(mean(data7s$c_100_cm2))*10000
#(mean(data7s$c_100_cm2))
#(sd(data7s$c_100_cm2))*10000



### Predicting using model 11 - use newdata5 from above

predicted_vals11 <- predict(mod11, new_data5, level = 0:2)

pred_vals11 <- predicted_vals11 %>%
  rename(log_c_dens_cm3 = predict.CoreName_2) %>%
  separate(CoreName_2, into = c("Site", "CoreName_2"),
           sep ="/") %>%
  mutate(corr_segment_midpoint = rep(c(1:100), times = length(unique(DF3$CoreName_2)))) %>%
  mutate(c_dens = exp(log_c_dens_cm3)) %>%
  mutate(c_stock_cm = c_dens * V_cseg) # c_stock in g per cm core length, these are not normally distributed now so be careful with means!

new_data11 <- new_data5 %>%
  left_join(pred_vals11) %>% 
  group_by(Site, CoreName_2, REI_Scaled, Type, Watercourse_NEAR_DIST.x, Coast) %>%
  mutate(c_60 = ifelse(corr_segment_midpoint < 61, c_stock_cm, 0)) %>%
  mutate(c_30 = ifelse(corr_segment_midpoint < 31, c_stock_cm, 0)) %>%
  summarise_at(., c("c_stock_cm", "c_60", "c_30"), sum) %>% #C g / inch core
  rename(c_100 = c_stock_cm) %>%
  mutate(c_60_cm2 = c_60/A_cseg) %>% # c_stock / cm2 
  mutate(c_100_cm2 = c_100/A_cseg) %>%
  mutate(c_30_cm2 = c_30/A_cseg) %>%
  mutate(c_60_m2 = c_60_cm2*10000) %>% # c_stock / cm2 
  mutate(c_100_m2 = c_100_cm2*10000) %>%
  mutate(c_30_m2 = c_30_cm2*10000)

write.csv(new_data11, file = "predicted11.csv")


## checking
data11 <- DF3 %>%
  group_by(CoreName_2) %>%
  mutate(c_stock_est = (c_dens*V_cseg)/A_cseg) %>%
  summarise(Mean_c = mean(c_stock_est)) %>%
  left_join(new_data11)

checking_100_11 <- ggplot(data11, aes(x = I(Mean_c*100), y = c_100_cm2)) +
  geom_point() +
  xlab("C Stock est by Mean C (gC / cm^2)") +
  ylab("C Stock est by model (gC / cm^2)") +
  ggtitle("Carbon Stock to 100 cm Mod11 (line is 1:1)") +
  geom_abline(intercept = 0, slope = 1)
checking_100_11
ggsave("Pred11 vs mean C_Stock 100cm.pdf", path = "./figures/", width = 4, height = 4)

checking_60_11 <- ggplot(data11, aes(x = I(Mean_c*60), y = c_60_cm2)) +
  geom_point() +
  xlab("C Stock est by Mean C (gC / cm^2)") +
  ylab("C Stock est by model (gC / cm^2)") +
  ggtitle("Carbon Stock to 60 cm mod11 (line is 1:1)") +
  geom_abline(intercept = 0, slope = 1)
checking_60_11
ggsave("Pred11 vs mean C stock 60cm.pdf", path = "./figures/", width = 4, height = 4)

checking_30_11 <- ggplot(data11, aes(x = (Mean_c*30), y = c_30_cm2)) +
  geom_point() +
  xlab("C Stock est by Mean C (gC / cm^2)") +
  ylab("C Stock est by model (gC / cm^2)") +
  ggtitle("Carbon Stock to 25 cm mod11 (line is 1:1)") +
  geom_abline(intercept = 0, slope = 1)
checking_30
ggsave("Pred11 vs mean C stock 25cm.pdf", path = "./figures/", width = 4, height = 4)


data11s <- subset(data11, Type == "SG")
(range(data11s$c_100_cm2))*10000
(range(data11s$c_60_cm2))*10000
(range(data11s$c_30_cm2))*10000

sitetypes <- data11 %>%
  group_by(Site, Coast, Mean_c) %>%
  summarise(coretype = (unique(Type)))

SG_100 <- data11 %>%
  group_by(Type) %>%
  summarise(mean_c100 = (mean(log(c_100_m2)))) %>%
  mutate(bt_mean = exp(mean_c100))
SG_100

coasts_100 <- data11 %>%
  group_by(Coast, Type) %>%
  summarise(mean_c100 = (mean(log(c_100_m2)))) 
coasts_100

coasts_100_2 <- data11 %>%
  group_by(Coast, Type) %>%
  summarise(sd_c100 = (sd(log(c_100_m2)))) %>%
  left_join(coasts_100) %>%
  mutate(BT_m = exp(mean_c100)) %>%
  mutate(BT_sdh = exp(mean_c100 + sd_c100)) %>%
  mutate(BT_sdl = exp(mean_c100 - sd_c100))

coasts_100_2

coasts_60 <- data11 %>%
  group_by(Coast, Type) %>%
  summarise(mean_c60 = (mean(log(c_60_m2)))) 

coasts_60_2 <- data11 %>%
  group_by(Coast, Type) %>%
  summarise(sd_c60 = (sd(log(c_60_m2)))) %>%
  left_join(coasts_60) %>%
  mutate(BT_m = exp(mean_c60)) %>%
  mutate(BT_sdh = exp(mean_c60 + sd_c60)) %>%
  mutate(BT_sdl = exp(mean_c60 - sd_c60))

coasts_60_2

coasts_30 <- data11 %>%
  group_by(Coast, Type) %>%
  summarise(mean_c30 = (mean(log(c_30_m2)))) 

coasts_30_2 <- data11 %>%
  group_by(Coast, Type) %>%
  summarise(sd_c30 = (sd(log(c_30_m2)))) %>%
  left_join(coasts_30) %>%
  mutate(BT_m = exp(mean_c30)) %>%
  mutate(BT_sdh = exp(mean_c30 + sd_c30)) %>%
  mutate(BT_sdl = exp(mean_c30 - sd_c30))

coasts_30_2

atlantic_hist <- ggplot(subset(data11, Coast == "Atlantic"), aes(log(Mean_c))) +
                          geom_histogram(aes(fill=factor(Type)))
atlantic_hist 

pacific_hist <- ggplot(subset(data11, Coast == "Pacific"), aes(log(Mean_c))) +
  geom_histogram(aes(fill=factor(Type)))
pacific_hist

gulf_hist <- ggplot(subset(data11, Coast == "Gulf of St.Lawrence"), aes(log(Mean_c))) +
  geom_histogram(aes(fill=factor(Type)), binwidth = 0.0001)
gulf_hist

#table 4 - try with raw values and not predicted? trying to reconcile model with higher carbon in bare, wiht results of more in SG in the coefficient...
mean(data11[data11$Coast == "Atlantic" & data11$Type == "BA", ]$c_100_cm2)*10000
mean(data11[data11$Coast == "Atlantic" & data11$Type == "SG", ]$c_100_cm2)*10000
mean(data11[data11$Coast == "Pacific" & data11$Type == "BA", ]$c_100_cm2)*10000
mean(data11[data11$Coast == "Pacific" & data11$Type == "SG", ]$c_100_cm2)*10000
mean(data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "BA", ]$c_100_cm2)*10000
mean(data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "SG", ]$c_100_cm2)*10000

sd((data11[data11$Coast == "Atlantic" & data11$Type == "BA", ]$c_100_cm2)*10000)
sd((data11[data11$Coast == "Atlantic" & data11$Type == "SG", ]$c_100_cm2)*10000)
sd((data11[data11$Coast == "Pacific" & data11$Type == "BA", ]$c_100_cm2)*10000)
sd((data11[data11$Coast == "Pacific" & data11$Type == "SG", ]$c_100_cm2)*10000)
sd((data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "BA", ]$c_100_cm2)*10000)
sd((data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "SG", ]$c_100_cm2)*10000)

mean(data11[data11$Coast == "Atlantic" & data11$Type == "BA", ]$c_60_cm2)*10000
mean(data11[data11$Coast == "Atlantic" & data11$Type == "SG", ]$c_60_cm2)*10000
mean(data11[data11$Coast == "Pacific" & data11$Type == "BA", ]$c_60_cm2)*10000
mean(data11[data11$Coast == "Pacific" & data11$Type == "SG", ]$c_60_cm2)*10000
mean(data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "BA", ]$c_60_cm2)*10000
mean(data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "SG", ]$c_60_cm2)*10000

sd((data11[data11$Coast == "Atlantic" & data11$Type == "BA", ]$c_60_cm2)*10000)
sd((data11[data11$Coast == "Atlantic" & data11$Type == "SG", ]$c_60_cm2)*10000)
sd((data11[data11$Coast == "Pacific" & data11$Type == "BA", ]$c_60_cm2)*10000)
sd((data11[data11$Coast == "Pacific" & data11$Type == "SG", ]$c_60_cm2)*10000)
sd((data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "BA", ]$c_60_cm2)*10000)
sd((data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "SG", ]$c_60_cm2)*10000)

mean(data11[data11$Coast == "Atlantic" & data11$Type == "BA", ]$c_30_cm2)*10000
mean(data11[data11$Coast == "Atlantic" & data11$Type == "SG", ]$c_30_cm2)*10000
mean(data11[data11$Coast == "Pacific" & data11$Type == "BA", ]$c_30_cm2)*10000
mean(data11[data11$Coast == "Pacific" & data11$Type == "SG", ]$c_30_cm2)*10000
mean(data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "BA", ]$c_30_cm2)*10000
mean(data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "SG", ]$c_30_cm2)*10000

sd((data11[data11$Coast == "Atlantic" & data11$Type == "BA", ]$c_30_cm2)*10000)
sd((data11[data11$Coast == "Atlantic" & data11$Type == "SG", ]$c_30_cm2)*10000)
sd((data11[data11$Coast == "Pacific" & data11$Type == "BA", ]$c_30_cm2)*10000)
sd((data11[data11$Coast == "Pacific" & data11$Type == "SG", ]$c_30_cm2)*10000)
sd((data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "BA", ]$c_30_cm2)*10000)
sd((data11[data11$Coast == "Gulf of St.Lawrence" & data11$Type == "SG", ]$c_30_cm2)*10000)


# LOI ---------------------------------------------------------------------

## merge files to see if the data match
# remove duplicated row
DF_LOI <- DF_LOI %>%
  filter(X != "72")

DFF <- DF_LOI %>%
  left_join(DF3, by = c("CoreName_2", "Extracted_IntervalStart_Depth_cm")) %>%
  mutate(LOI_Percent = LOI_Percent/100) %>%
  mutate(OC_Per = OC_Per.x/100) %>%
  filter(SiteCode.x != "SIp")
  
plot(log(DFF$LOI_Percent) ~ log(DFF$c_dens))
plot(log(DFF$LOI_Percent) ~ log(DFF$OC_Per.x))
plot(log(DFF$LOI_Percent) ~ log(DFF$c_stock))
plot(log(DFF$OC_Per.y) ~ log(DFF$OC_Per.x))   #whew!
   
lm(log(DFF$LOI_Percent) ~ log(DFF$c_dens))
lm(log(DFF$LOI_Percent) ~ log(DFF$OC_Per.x))

## facet plot of LOI as possible predictor of C_org (Percent dry weight). c_dens, c_stock and OC_Per

# c_percent
mod1 <- lm((DFF$OC_Per) ~ (DFF$LOI_Percent))
summary(mod1)

LOI_fun <- function(x) coef(mod1)[2] * (x) + coef(mod1)[1]

LOI_c_per <- ggplot(DFF, aes(x = LOI_Percent, y = OC_Per)) +
  geom_point() + 
  theme_bw() +
  geom_function(fun = LOI_fun) +
  geom_abline(intercept = 0, slope = 1)

LOI_c_per

# c_percent_transformed
mod2 <- lm(sqrt(DFF$OC_Per) ~ sqrt(DFF$LOI_Percent))
summary(mod2)

LOI_fun2 <- function(x) coef(mod2)[2] * (x) + coef(mod2)[1]

LOI_OCPer2 <- ggplot(DFF, aes(x = sqrt(LOI_Percent), y = sqrt(OC_Per))) +
  geom_point() + 
  theme_bw() +
  geom_function(fun = LOI_fun2) +
  geom_abline(intercept = 0, slope = 1)
  
LOI_OCPer2

# c_percent_transformed
mod3 <- lm(log(DFF$OC_Per) ~ log(DFF$LOI_Percent))
summary(mod3)

LOI_fun3 <- function(x) coef(mod3)[2] * (x) + coef(mod3)[1]

LOI_OCPer3 <- ggplot(DFF, aes(x = log(LOI_Percent), y = log(OC_Per))) +
  geom_point() + 
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) +
  geom_smooth(method=lm , color="black", se=TRUE) +
  #geom_function(fun = LOI_fun) +
  annotate("text", x = -4.5, y=-3.5, label = paste("R^2 == 0.88"), parse = TRUE) +
  annotate("text", x = -4.5, y=-3, label = ("y = 1.33x - 0.08")) +
  xlab("Loss on Ignition (LOI) ln(%)") +
  ylab("Organic Carbon from Elemental Analysis ln(%)") #+
 # ggtitle("Carbon estimated by Elemental Analysis and LOI") +
 # geom_function(fun = LOI_fun3, color = "red") #+
  #geom_abline(intercept = 0, slope = 1)

LOI_OCPer3
ggsave("OC_LOI.tiff", path = "./figures/", width = 4, height = 4)

# quick model comparison

#remove SiteCode.x = SIp

mod30 <- lme(log(OC_Per) ~ 1 + log(LOI_Percent), data = DFF, random = ~ 1 | SiteName2/CoreName_2)

mod3a <- lme(log(OC_Per) ~ 1 + log(LOI_Percent) + sqrt(REI_Scaled), data = DFF, random = ~ 1 | SiteName2/CoreName_2)  

mod3b <- lm(log(DFF$OC_Per) ~ 1 + log(LOI_Percent) * sqrt(REI_Scaled), data = DFF, random = ~ 1 | SiteName2/CoreName_2)  

mod3c <- lm(log(DFF$OC_Per) ~ log(DFF$LOI_Percent)*sqrt(DFF$REI_Scaled))

mod3d <- lm(log(DFF$OC_Per) ~ log(DFF$LOI_Percent)*log(DFF$Watercourse_NEAR_DIST.x))

model.sel(mod30, mod3a)
View(DFF)

hist(DFF$LOI_Percent)
     
length(unique(DF3$CoreName_2))
length(unique(DF3$SiteCode))

0.263^2
#estimate marginal SG effect / cm2 surface area

exp(0.088) * V_cseg * 100 / A_cseg
exp(0.088) * V_cseg * 25 / A_cseg

## mean of analyzed c_dens
mean(log(DF3$c_dens)) + 0.088
mean(log(DF3$c_dens))

exp(mean(log(DF3$c_dens)) + 0.088) - exp(mean(log(DF3$c_dens)))

exp(min(log(DF3$c_dens)) + 0.088) - exp(min(log(DF3$c_dens)))

exp(max(log(DF3$c_dens)) + 0.088) - exp(max(log(DF3$c_dens)))

# coef from model 18 for SG is 0.088. so the backtransformed value would be exp(0.088), or 1.09 g C / cm3 with SG vs without. The problem is that mean raw values for gC/cm3 range from 6.9 x 10 ^-6 - 5.66 x 10 ^-2, so this 1 unit addition is way out of proportion. 
## ok so mean (log(c_dens)) = -5.69, and the 0.088 marginal effect on that is reasonable. right so we can't backtransform the effect size so easily. 