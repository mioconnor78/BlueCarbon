## Mary's analysis of Matt's data Feb 2025

## packages
library(tidyverse)
library(nlme)
library(MuMIn)
library(car)
library(ggplot2)
library(effects)


#Read in data ----
DF <- read.csv(file = "./data/WongChristensenDataforModellingOct19.csv")
DF_LOI <- read.csv(file = "./data/BlueCarbonData_1.csv")
CF <- read.csv(file = "./data/compaction.csv")
#DF$intertidal <- ifelse(DF$Elevation > "0", "intertidal", "subtidal")
#DF$Coast_no <- DF[, DF$Coast == "Atlantic"] 
DF$c_dens <- (DF$OC_Per/100)*DF$Corrected_DBD_g_cm3 # C (g / cm3)
DF$Site <- factor(DF$Site)
DF$Type <- factor(DF$Type)
#DF <- DF %>%
 # filter(CoreName_2 != "PEI-CAR-1-BA-") ## DBD values are weird, i can probably figure this out if i try harder; update, i fixed them


## notes: 
# Corrected_DBD_g_cm3 is already corrected for compaction. i checked Matt's PEI file and it looks good. working on Cardigan bare core 1 though that's weird. 

## DBD is a weight, normalized / cm3 by estimating the cylinder size the dirt came from. so, DBD g/cm3. c_dens is then g/cm3 b/c we just multiplied by a percent. So then we have c_dens for each cm of depth. to get the core, we add those c_dens values. so this is still g/cm3. i think i have to multiple the c_dens at each cm by the volume of that 1 cm core segment. so, x g/cm3 * 1 cm x (2*pi*2.38125) cm2 gives the s

V_cseg <- pi*(2.38125^2) # volume of 1 cm of core sediment, cm3
A_cseg <- 2*pi*2.38125 # area of top of core, cm2

## Now, maybe i don't need compaction factors anymore. 


# Data cleaning and reconstruction  ---------------------------------------

#remove row 430, appears to be a duplicate

## reconstructing core compaction correction factors

## we can't find Matt's original core length data, but based on his thesis, he used a single correction factor for each core (as opposed to a different correction factor for each segment). The data appear to use a slightly different correction factor for each segment. For the ones that produce negative values, i propose we adjust those correction factors to match those for the rest of the core, assuming that that a correction factor >1 is an error.

DF$estCfactorS <- ifelse(DF$Extracted_IntervalStart_Depth_c != "0",  DF$Extracted_IntervalStart_Depth_cm/DF$Corrected_IntervalStart_Depth_cm, "0")

#DF$estCfactorE <- ifelse(DF$Extracted_IntervalEnd_Depth_c != "0",  DF$Extracted_IntervalEnd_Depth_cm/DF$Corrected_IntervalEnd_Depth_cm, "0")

#plot(DF$Extracted_IntervalStart_Depth_cm, DF$estCfactorS)

## discovered some weird numbers, looking into it. 
# create a datafile mean and sd est conversion factors to identify possible errors
DF_cores <- DF %>%
  group_by(SiteCode, Core.Number) %>%
  filter(estCfactorS != "0") %>%
  summarise(meanCF = mean(as.numeric(as.character(estCfactorS))), sdCF = sd(estCfactorS))

## replacing values for correction factors: [Mary double check this is still needed - March 4]
DF_test <- DF %>%
  mutate(estCfactorS = as.numeric(estCfactorS)) %>%
  mutate(estCfactorS = replace(estCfactorS, estCfactorS == 0.18018018018018, 0.869565217391304))

DF_cores <- DF_test %>%
  group_by(SiteCode, Core.Number) %>%
  filter(estCfactorS != "0") %>%
  mutate(as.numeric(as.character(estCfactorS))) %>%
  summarise(meanCF = mean(estCfactorS), sdCF = sd(estCfactorS))

## spot check a couple of sites
DF_SAs <- DF_test %>%
  filter(SiteCode == "SAs") %>%
  select(CoreName_2, estCfactorS)
View(DF_SAs)

DF_BOB <- DF_test %>%
  filter(SiteCode == "BOB") %>%
  select(CoreName_2, estCfactorS)
View(DF_BOB)

# merge the mean CF factors back into the main data frame. the CF should be the same for the whole core; i can't find methods or code for how Matt would have estimated them for each segment. So we will proceed with a single CF factor for each core, using the mean CF value estimated above.

DF2 <- DF %>%
  left_join(DF_cores)

plot(DF2$meanCF, DF2$estCfactorS)
hist(DF2$meanCF). # still have some 2s there. 

## having identified problems, we found the original datasheets. bring in Melisa's compaction data and merge
CF2 <- CF %>%
  select(-c(CoreName_2, X, sample, compression))

DFc <- DF2 %>%
  left_join(CF2)

#View(DFc)

DFc2 <- DFc %>%
  mutate(estCfactorS = as.numeric(estCfactorS)) %>%
  mutate(conversion = ifelse(is.na(conversion), "Missing", conversion)) %>%
  mutate(meanCF = ifelse(is.na(meanCF), "25", meanCF)) %>%
  mutate(CF_final = as.numeric(ifelse(conversion == "Missing", meanCF, conversion))) %>%
  mutate(CF_final = ifelse(meanCF == "25", estCfactorS, CF_final)) 

## spot check sites
DF_SOO <- DFc2 %>%
  filter(SiteCode == "SOO") %>%
  select(CoreName_2, estCfactorS, meanCF, CF_final)
View(DF_SOO)

plot(DFc2$meanCF, DFc2$CF_final)
hist(DFc2$CF_final)

# COMPRESSION CORRECTION FACTORS FIXED. Now we have a datafile with no compaction corrections > 1. there are a few discrepancies between Matt's data and Melisa's original sheets, and we estimated values for Port Joli, Port l'Hebert and Taylor's head. 

## for later model prediction, i think we want to get a carbon stock / cm of core depth. the segments are not even values. This is from the methods: Each core was thawed and then extruded in segments (0-2 cm, 2-5 cm, 5-10 cm, 10-20 cm, 20-30 cm, 30-40 cm, 40-50 cm, 50-60 cm) up to 60 cm in depth. So i will divide the c_stock values by corrected core thickness to get a c_stock / cm of core depth.

## Create new corrected core segment values
DF3 <- DFc2 %>%
  mutate(corr_SD = Extracted_IntervalStart_Depth_cm * CF_final) %>%
  mutate(corr_ED = Extracted_IntervalEnd_Depth_cm * CF_final) %>%
  mutate(corr_thickness = corr_ED - corr_SD) %>% 
  mutate(corr_segment_midpoint = (corr_ED + corr_SD)/2) %>%
 # mutate(c_stock = c_dens * corr_thickness) %>% #C (g / cm2)
  select(!c(Elevation, X.y, Watercourse_NEAR_DIST.y, Corrected_Midpoint_cm_rounded, Corrected_IntervalStart_Depth_cm, Corrected_IntervalEnd_Depth_cm)) %>% #REI_Scaled,
  filter(X.1 != "142") %>%
  filter(X.1 != "429") #duplicated rows

View(DF3)
plot(DF3$corr_segment_midpoint, DF3$Corrected_Midpoint_cm)
hist(DF3$corr_thickness)
#DF3$c_stock <- DF3$c_dens * DF3$corr_thickness
#hist(log(DF3$c_stock))
plot(DF3$corr_segment_midpoint, DF3$corr_thickness)

## DF3 is the dataset to use. 
write.csv(DF3, "./data/BCdataforanalysis.csv")

# old troubleshooting code ------------------------------------------------


## start here to see if the corrections make sense (they look pretty good), and re-analyze data. need core lengths from Melisa to add in the atlantic

#between here and models is me figuring out data problems that were resolved in the code above

test <- DF %>%
  filter(grepl("PEI", CoreName_2, ignore.case = TRUE)) %>%
  select(CoreName_2, estCfactorS) %>%
  group_by(CoreName_2) %>% 
  filter(estCfactorS != "0") %>%
  summarise(mean = mean(as.numeric(as.character(estCfactorS))))

# fixed cardigan Bare 1 corrected values, found an error in the original google sheet (the wrong cell was identified for the division to estimate the correction. So this proved to me that unusual numbers (negative, > 1) are likely wrong)

# NEXT: need to recalculate midpoint and thickness once we're happy with the corrected depths; i think that will expose an issue with the BC conversions to end 

# trouble cases
problems <- DF %>%
  filter(estCfactorS > 1) %>%
  select(CoreName_2, estCfactorS, estCfactorE, Extracted_IntervalStart_Depth_cm, Corrected_IntervalStart_Depth_cm)

DF$thickness_new <- DF$Corrected_IntervalEnd_Depth_cm - DF$Corrected_IntervalStart_Depth_cm

problems <- DF %>%
  filter(thickness_new < 0) %>%
  select(CoreName_2, thickness_new, estCfactorS, estCfactorE, Extracted_IntervalStart_Depth_cm, Corrected_IntervalStart_Depth_cm)

problems <- DF %>%
  filter(estCfactorE != estCfactorS) %>%
  select(CoreName_2, thickness_new, estCfactorS, estCfactorE, Extracted_IntervalStart_Depth_cm, Corrected_IntervalStart_Depth_cm)

DFtest <- DF %>%
  filter(estCfactorE > 0.4)
  
testplot <- ggplot(DF3, aes(x = estCfactorE, y = estCfactorS)) + 
  geom_point() +
  geom_text(
    label=rownames(DF3), 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T
  )
testplot

DF_CRB <- DF %>%
  filter(CoreName_2 == "BC-CRB-2-BA-2")





# Modeling Carbon Density -------------------------------------------------

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

mod23 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod24 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + log(Watercourse_NEAR_DIST.x) + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod25 <- lme(log(c_dens) ~1 + sqrt(REI_Scaled)*Type + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")


model.sel(mod0, mod1, mod2,mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod18a, mod19, mod20,  mod23, mod24, mod25)

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

# Example of predicting ---------------------------------------------------

## predicting: example: https://rdrr.io/cran/nlme/man/predict.lme.html
# and https://stats.stackexchange.com/questions/234618/meaning-of-predict-fixed-and-predict-subject-values-in-predict-function-of-nlm

# I need the 'new dataset' with the fixed effects levels for each site, core, and add in depths of 25, 60 and 100, and have the model give me those c_dens values for those cores.

View(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
newOrth <- data.frame(Sex = c("Male","Male","Female","Female","Male","Male"),
                      age = c(15, 20, 10, 12, 2, 4),
                      Subject = c("M01","M01","F30","F30","M04","M04"))
## The 'Orthodont' data has *no* 'F30', so predict  NA  at level 1 :
predict(fm1, newOrth, level = 0:1)



# Predicting C_stocks from C_dens models -------------------------------------------------------

#C_density

## create the new data frame with the depths we want to predict C_dens and eventually C_stock for: 1 - 100 cm
 
new_data4 <- data.frame(CoreName_2 = rep(unique(DF3$CoreName_2), each = 100), corr_segment_midpoint = rep(c(1:100), times = length(unique(DF3$CoreName_2)))) 
                      
#merge new_data with DF3 to get site level variables: Site, REI_Raw, Type, Watercourse_NEAR_DIST.x
#how to get Percent.Silt.Fraction? this is a segment level variable. use mean silt val for core. 

Site_info <- DF3 %>% 
  group_by(Site, CoreName_2, REI_Scaled, Type, Watercourse_NEAR_DIST.x) %>%
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
  mutate(c_stock_cm = c_dens * V_cseg) # c_stock in g per cm core length

new_data6 <- new_data5 %>%
  left_join(pred_vals) %>% 
  group_by(Site, CoreName_2, REI_Scaled, Type, Watercourse_NEAR_DIST.x) %>%
  mutate(c_60 = ifelse(corr_segment_midpoint < 61, c_stock_cm, 0)) %>%
  mutate(c_25 = ifelse(corr_segment_midpoint < 26, c_stock_cm, 0)) %>%
  summarise_at(., c("c_stock_cm", "c_60", "c_25"), sum) %>% #C g / inch core
  rename(c_100 = c_stock_cm) %>%
  mutate(c_60_cm2 = c_60/A_cseg) %>% # c_stock / cm2
  mutate(c_100_cm2 = c_100/A_cseg) %>%
  mutate(c_25_cm2 = c_25/A_cseg)

write.csv(new_data6, file = "predicted.csv")


## checking
data7 <- DF3 %>%
  group_by(CoreName_2) %>%
  mutate(c_stock_est = (c_dens*V_cseg)/A_cseg) %>%
  summarise(Mean_c = mean(c_stock_est)) %>%
  left_join(new_data6)

hist(data7$c_100) # g carbon / 100 cm core of diameter 17 cm2, or
hist(data7$c_100_cm2) 
hist(data7$c_60_cm2) 
hist(data7$c_25_cm2)


c_stock <- ggplot(data7, aes(x = c_dens, y = c_stock_cm)) +
  geom_point() +
  xlab("C_dens (g / cm3))") +
  ylab("C Stock (g)") +
  ggtitle("Carbon Stock per cm of core") +
  geom_abline(intercept = 0, slope = 1)
c_stock
  
checking_100 <- ggplot(data7, aes(x = Mean_c, y = c_100_cm2)) +
  geom_point() +
  xlab("C Stock est by Mean C (gC / cm^2)") +
  ylab("C Stock est by model (gC / cm^2)") +
  ggtitle("Carbon Stock to 100 cm (line is 1:1)") +
  geom_abline(intercept = 0, slope = 1)
checking_100
ggsave("Pred vs mean C_Stock 100cm.pdf", path = "./figures/", width = 4, height = 4)

checking_60 <- ggplot(data7, aes(x = Mean_c*.60, y = c_60_cm2)) +
  geom_point() +
  xlab("C Stock est by Mean C (gC / cm^2)") +
  ylab("C Stock est by model (gC / cm^2)") +
  ggtitle("Carbon Stock to 60 cm (line is 1:1)") +
  geom_abline(intercept = 0, slope = 1)
checking_60
ggsave("Pred vs mean C stock 60cm.pdf", path = "./figures/", width = 4, height = 4)

checking_25 <- ggplot(data7, aes(x = Mean_c*.25, y = c_25_cm2)) +
  geom_point() +
  xlab("C Stock est by Mean C (gC / cm^2)") +
  ylab("C Stock est by model (gC / cm^2)") +
  ggtitle("Carbon Stock to 25 cm (line is 1:1)") +
  geom_abline(intercept = 0, slope = 1)
checking_25
ggsave("Pred vs mean C stock 25cm.pdf", path = "./figures/", width = 4, height = 4)



# Modeling carbon stock ---------------------------------------------------

## i'm not sure if it make sense to analyze carbon stock here, b/c i estimated it above from the c_dens model already. should we do both?

# Modeling Carbon Stock with ranef for depth effect
mod0.1 <- lme(log(c_stock + 0.01) ~ 1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ 1 | CoreName_2, method = "ML")

mod0.3 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ 1 | Site/CoreName_2, method = "ML")

mod0.0 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ corr_segment_midpoint | CoreName_2, method = "ML")

mod0.2 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~ corr_segment_midpoint | Site/CoreName_2, method = "ML")


anova(mod0.2, mod0.3, mod0.0, mod0.1)

##model comparison
mod0 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod1 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

## all site level predictors
mod2 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + log(Watercourse_NEAR_DIST.x) + Coast, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod3 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod4 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + log(Watercourse_NEAR_DIST.x)*Type + Coast, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod5 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + log(Watercourse_NEAR_DIST.x), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod6 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod7 <- lme(log(c_stock + 0.01) ~1, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

## all site + core level predictors

mod8 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + (corr_segment_midpoint), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod9 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod11 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod12 <- lme(log(c_stock + 0.01) ~1 + Type + log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod13 <- lme(log(c_stock + 0.01) ~1 + Type + log(Watercourse_NEAR_DIST.x) + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod14 <- lme(log(c_stock + 0.01) ~1 + Type + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod15 <- lme(log(c_stock + 0.01) ~1 + Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod16 <- lme(log(c_stock + 0.01) ~1 + Type*log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod17 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod18 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

#mod18a <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction) + corr_segment_midpoint, data = DF, random = ~1| Site/CoreName_2, method = "REML")

mod19 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x)  + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod20 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x)  + corr_segment_midpoint, data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod23 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + Coast + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod24 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")

mod25 <- lme(log(c_stock + 0.01) ~1 + log(REI_Raw)*Type + corr_segment_midpoint + sqrt(Percent.Silt.Fraction), data = DF3, random = ~corr_segment_midpoint | Site/CoreName_2, method = "REML")


model.sel(mod0, mod1, mod2,mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod19, mod20,  mod23, mod24, mod25)

## comparing predicted vs observed values
plot(predict(mod11), log(DF3$c_stock))


# Predicting C_stocks from C_stock models -------------------------------------------------------

#this is not quite working; i haven't dealt with the +0.01 yet from the model. i'm too tired to figure this out now, so don't use this yet!

## create the new data frame with the depths we want to predict C_stock for: 1 - 100 cm

new_data4s <- data.frame(CoreName_2 = rep(unique(DF3$CoreName_2), each = 100), corr_segment_midpoint = rep(c(1:100), times = 84)) 

#merge new_data with DF3 to get site level variables: Site, REI_Raw, Type, Watercourse_NEAR_DIST.x
#how to get Percent.Silt.Fraction? this is a segment level variable. use mean silt val for core. 

Site_info_s <- DF3 %>% 
  group_by(Site, CoreName_2, REI_Raw, Type, Coast, Watercourse_NEAR_DIST.x) %>%
  summarise(mud = mean(Percent.Silt.Fraction)) 

#use site level variables from above
new_data5s <- new_data4s %>%
  left_join(Site_info_s) %>%
  rename(Percent.Silt.Fraction = mud) %>%
  mutate(corr_segment_midpoint = as.numeric(corr_segment_midpoint))

predicted_vals <- predict(mod11, new_data5s, level = 0:2)

pred_vals_s <- predicted_vals %>%
  separate(CoreName_2, into = c("Site", "CoreName_2"),
           sep ="/") %>%
  mutate(corr_segment_midpoint = rep(c(1:100), times = 84)) %>%
  mutate(c_stock = exp(predict.CoreName_2))

new_data6s <- new_data5s %>%
  left_join(pred_vals_s) %>% 
  group_by(Site, CoreName_2, REI_Raw, Type, Watercourse_NEAR_DIST.x) %>%
  mutate(c_60s = ifelse(corr_segment_midpoint < 61, c_stock, 0)) %>%
  mutate(c_25s = ifelse(corr_segment_midpoint < 26, c_stock, 0)) %>%
  summarise_at(., c("c_stock", "c_60s", "c_25s"), sum)

write.csv(new_data6, file = "predicted.csv")


## checking
data7s <- DF3 %>%
  group_by(CoreName_2) %>%
  summarise(Mean_s = mean(c_stock)) %>%
  left_join(new_data6s)

checking_100_s <- ggplot(data7s, aes(x = Mean_s*100, y = c_stock)) +
  geom_point() +
  xlab("Mean C gC / cm^2") +
  ylab("C stock 100 cm") +
  geom_abline(intercept = 0, slope = 1)
checking_100_s
ggsave("Pred vs mean C_stock 100cm.pdf", path = "./figures/", width = 4, height = 4)

checking_60_s <- ggplot(data7s, aes(x = Mean_s*60, y = c_60s)) +
  geom_point() +
  xlab("Mean C gC / cm^2") +
  ylab("C stock 60 cm") +
  geom_abline(intercept = 0, slope = 1)
checking_60_s
ggsave("Pred vs mean C_stock 60cm.pdf", path = "./figures/", width = 4, height = 4)

checking_25_s <- ggplot(data7s, aes(x = Mean_s*25, y = c_25s)) +
  geom_point() +
  xlab("Mean C gC / cm^3") +
  ylab("C stock 25 cm") +
  geom_abline(intercept = 0, slope = 1)
checking_25_s
ggsave("Pred vs mean C_stock 25cm.pdf", path = "./figures/", width = 4, height = 4)




# Modeling OC_Per ---------------------------------------------------------
## this is out of date now

#Testing all parameters
mod0 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod1 <- lme(log(OC_Per) ~1 + log(REI_Raw) + Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

## all site level predictors
mod2 <- lme(log(OC_Per) ~1 + log(REI_Raw) + intertidal + log(Watercourse_NEAR_DIST.x) + Coast, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod3 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod4 <- lme(log(OC_Per) ~1 + log(REI_Raw) + intertidal + log(Watercourse_NEAR_DIST.x)*Type + Coast, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod5 <- lme(log(OC_Per) ~1 + log(REI_Raw) + intertidal + log(Watercourse_NEAR_DIST.x), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod6 <- lme(log(OC_Per) ~1 + log(REI_Raw) + intertidal, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod7 <- lme(log(OC_Per) ~1 + log(REI_Raw), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod8 <- lme(log(OC_Per) ~1 + intertidal, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

## all site + core level predictors

mod9 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod10 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3) + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod11 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod12 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod13 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod28 <- lme(log(OC_Per) ~1 + Type + log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod29 <- lme(log(OC_Per) ~1 + Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod30 <- lme(log(OC_Per) ~1 + Type + Coast + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod32 <- lme(log(OC_Per) ~1 + Type + log(Corrected_DBD_g_cm3) + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod31 <- lme(log(OC_Per) ~1 + Type*log(Watercourse_NEAR_DIST.x) + Coast + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod27 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod14 <- lme(log(OC_Per) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x) + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod15 <- lme(log(OC_Per) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x)  + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod16 <- lme(log(OC_Per) ~1 + log(REI_Raw) + Type + log(Watercourse_NEAR_DIST.x)  + Corrected_Midpoint_cm_rounded, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod17 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + Corrected_Midpoint_cm_rounded, data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod18 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod19 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + Coast + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod20 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Watercourse_NEAR_DIST.x) + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod21 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + intertidal + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod22 <- lme(log(OC_Per) ~1 + log(REI_Raw)*Type + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod23 <- lme(log(OC_Per) ~1 + log(REI_Raw) + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod24 <- lme(log(OC_Per) ~1 + log(REI_Raw) + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod25 <- lme(log(OC_Per) ~1 + log(REI_Raw) + log(Watercourse_NEAR_DIST.x) + Coast + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

mod26 <- lme(log(OC_Per) ~1 + log(REI_Raw) + log(Watercourse_NEAR_DIST.x) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF, random = ~Corrected_Midpoint_cm | Site/CoreName_2, method = "REML")

model.sel(mod0,mod1, mod2,mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod19, mod20, mod21, mod22, mod23, mod24, mod25, mod26, mod27, mod28, mod29, mod30, mod31, mod32)


mod0 <- lm(log(OC_Per) ~1 + log(REI_Raw) + Type + intertidal + log(Watercourse_NEAR_DIST.x) + Coast + log(Corrected_DBD_g_cm3) + Corrected_Midpoint_cm_rounded + sqrt(Percent.Silt.Fraction), data = DF)

vif_values <- vif(mod0)

cor_data <- DF[, c("REI_Raw", "Watercourse_NEAR_DIST.x", "Corrected_DBD_g_cm3", "Corrected_Midpoint_cm_rounded", "Percent.Silt.Fraction")]

cor(cor_data)


# Plotting ----------------------------------------------------------------

depth_fun <- function(z) exp(0.35/z) -1

Cdens_coredepth_Ann <- ggplot(subset(DF3, Site == "Annandale"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
  scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2, ncol = 6)
Cdens_coredepth_Ann

Cdens_coredepth_Air <- ggplot(subset(DF3, Site == "Airpark"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
  scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2, ncol = 6)
Cdens_coredepth_Air

Cdens_coredepth_Car <- ggplot(subset(DF3, Site == "Cardigan"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
  scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2, ncol = 6)
Cdens_coredepth_Car

Cdens_coredepth_Cat <- ggplot(subset(DF3, Site == "Cates Park"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
  scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2, ncol = 6)
Cdens_coredepth_Cat

Cdens_coredepth_Cen <- ggplot(subset(DF3, Site == "Centennial Beach"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
  scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2, ncol = 6)
Cdens_coredepth_Cen

Cdens_coredepth_Cre <- ggplot(subset(DF3, Site == "Crescent Beach"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
 scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2)
Cdens_coredepth_Cre

Cdens_coredepth_Ful <- ggplot(subset(DF3, Site == "Fulmore"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
  scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2)
Cdens_coredepth_Ful

Cdens_coredepth_MIS <- ggplot(subset(DF3, Site == "Mason's Island Shallow"), aes(x = corr_segment_midpoint, y = c_dens)) +
  geom_point() +
  scale_x_reverse() +
  geom_function(fun = depth_fun) +
  facet_wrap( ~ CoreName_2)
Cdens_coredepth_MIS



# LOI ---------------------------------------------------------------------

## merge files to see if the data match

# remove duplicated row
DF_LOI <- DF_LOI %>%
  filter(X != "72")

DFF <- DF3 %>%
  left_join(DF_LOI, by = c("CoreName_2", "Extracted_IntervalStart_Depth_cm"))
  
plot(log(DFF$LOI_Percent) ~ log(DFF$c_dens))
plot(log(DFF$LOI_Percent) ~ log(DFF$OC_Per.x))
plot(log(DFF$LOI_Percent) ~ log(DFF$c_stock))
plot(log(DFF$OC_Per.y) ~ log(DFF$OC_Per.x))   #whew!
   
lm(log(DFF$LOI_Percent) ~ log(DFF$c_dens))
lm(log(DFF$LOI_Percent) ~ log(DFF$OC_Per.x))


## facet plot of LOI as possible predictor of c_dens, c_stock and OC_Per

# c_dens logged
mod1 <- lm(log(DFF$c_dens) ~ log(DFF$LOI_Percent))
summary(mod1)

LOI_fun <- function(x) coef(mod1)[2] * (x) + coef(mod1)[1]

LOI_c_dens <- ggplot(DFF, aes(x = log(LOI_Percent), y = log(c_dens))) +
  geom_point() + 
  theme_bw() +
  geom_function(fun = LOI_fun)
  
LOI_c_dens


#c_dens untransformed
mod1 <- lm((DFF$c_dens) ~ (DFF$LOI_Percent))
summary(mod1)

LOI_fun <- function(x) coef(mod1)[2] * (x) + coef(mod1)[1]

LOI_c_dens <- ggplot(DFF, aes(x = (LOI_Percent), y = (c_dens))) +
  geom_point() + 
  theme_bw() +
  geom_function(fun = LOI_fun)

LOI_c_dens

# c_stock logged
mod2 <- lm(log(DFF$c_stock + 0.01) ~ log(DFF$LOI_Percent))
summary(mod2)

LOI_fun <- function(x) coef(mod2)[2] * (x) + coef(mod2)[1]

LOI_c_stock <- ggplot(DFF, aes(x = log(LOI_Percent), y = log(c_stock + 0.01))) +
  geom_point() + 
  theme_bw() +
  geom_function(fun = LOI_fun)

LOI_c_stock


# OC_Per logged
mod3 <- lm(log(DFF$OC_Per.x) ~ log(DFF$LOI_Percent))
summary(mod3)

LOI_fun <- function(x) coef(mod3)[2] * (x) + coef(mod3)[1]

LOI_c_OC <- ggplot(DFF, aes(x = log(LOI_Percent), y = log(OC_Per.x))) +
  geom_point() + 
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) +
  geom_smooth(method=lm , color="black", se=TRUE) +
  #geom_function(fun = LOI_fun) +
  ylab("Percent Organic Carbon (ln(%))") +
  xlab("LOI estimate (ln(%))") +
  ylim(-4, 2) +
  annotate("text", x = 1, y=2, label = paste("R^2 == 0.88"), parse = TRUE) +
  annotate("text", x = 0, y=2, label = ("y = 1.34x - 1.59"))

LOI_c_OC

ggsave("LOI_OC.pdf", path = "./figures", width = 4, height = 4)

## ok next steps: 
# make a clean plot, with appropriate scale - done
# this is an exact repeat of Matt's; is that right?
# methods text
# predicted values (do this first)
     