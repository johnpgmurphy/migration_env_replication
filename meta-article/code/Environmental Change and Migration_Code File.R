###  
### A Meta-Analysis of Country Level Studies on Environmental Change and Migration
###
### AUTHORS: Hoffmann, R., Dimitrova, A., Muttarak, R., Crespo-Cuaresma, J., Peisker, J. 
###

# >>> Set your working directory here
setwd("~/Desktop/Gov52/rep_project/meta-article")
# >>> Save the main RData files in the working directory to access them with this code

### **************************************************************
### OUTLINE ----------------------------------------------------
### **************************************************************
 
# 1. PREPARATION
# 2. OUTPUT TABLES AND FIGURES FOR MAIN PAPER
# 3. OUTPUT TABLES AND FIGURES FOR EXTENDED DATA
# 4. OUTPUT TABLES AND FIGURES FOR SUPPLEMENTARY MATERIAL


#
##
### **************************************************************
### 1. PREPARATIONS ---------------------------------------------- 
### **************************************************************
##
#


##
## DATA AND VARIABLE DESCRIPTIONS  ******************************** 
## 

## SL = study line / individual model

# DATA SETS       | DESCRIPTION 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# meta            | main meta data file, contains information on  1803 study line cases
# meta.plus       | main meta file plus coefficients estimated for environmental changes in *destination* countries
# country         | country data set with information on country characteristics and environmental change observed in countries for predictions
# country_rep     | NUmber of times countries are included in the country samples

# VARIABLE        | DESCRIPTION                       | SCALE   | CATEGORIES 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - # "paper_short"   | paper identifier                  | factor
# "paper"         | paper identifier                  | factor
# "author_id      | author pair identifier            | factor
# "Table"         | table identifier                  | factor
# "year"          | year of publication               | continuous

# stancoef        | standardized coefficient          | continuous
# stancoef.abs    | absolute standardized coefficient | continuous
# stancoef.w      | weighted standardized coefficient | continuous
# stanse          | standardized standard errors      | continuous
# coef            | orig. coefficient from study line | continuous
# serror          | orig. stan. error from study line | continuous
# mig_sd          | stan. deviation of migration var. | continuous
# env_sd          | stan. deviation of env. var.      | continuous

# "fe_time"       | study uses time fixed effects     | dummy   | 0=no fixed effects, 1=fixed effects
# "fe_spatial"    | study uses spatial fixed effects  | dummy   | 0=no fixed effects, 1=fixed effects
# interaction"    | effect size based on interation   | factor  | 0=no interaction, 1=interaction with dummy, effect is main coefficient, 2=interaction with dummy, combined effect, 3=interaction with continuous variable, combined effect
# "tab"           | table count within paper          | continuous
# "countrysample" | number of countries in sample     | continuous 
# "yearscovered"  | number of years in panel          | continuous
# "lin"           | all variables in linear form      | dummy
# "robustse"      | estimates (cluster) robust SE     | dummy  

# "rapidonset"    | SL considers rapid-onset event... | dummy   | 0=other env. factor, 1=rapid-onset event
# "pre_lev"       | ... precipitation level change    | dummy   | 0=other env. factor, 1=prec. level change
# "pre_var"       | ... precipitation variability     | dummy   | 0=other env. factor, 1=prec. variability change/anomaly
# "tem_lev"       | ... temperature level change      | dummy   | 0=other env. factor, 1=temp level change
# "tem_var",      | ... temperature variability       | dummy   | 0=other env. factor, 1=temp. variability change/anomaly
# "precipitation" | SL considers prec. level/variab   | dummy   | 0=other env. factor, 1=precipitation 
# "temperature"   |  SL considers temp. level/variab  | dummy   | 0=other env. factor, 1=temperature 
# "env_lag"       | lag between ind. and dep. var.    | continuous
# "env_lag_dum    | lag between ind. and dep. var.    | dummy   | 0=no lag, 1 = lag > 0
# "env_other_sum" | other env. factor controlled for  | count
# "env_other_dum" | other env. factor controlled for  | dummy
# "env_timespan"  | time interval to measure env. var | factor  | 1=1 year, 2=5 years, 3=10years
# "env_timespan_dum" | interval to measure env. var   | dummy   | 0=1 year, 1=>1 year 
# "env_others"    | other env. var. controlled for    | string
# "o_tem"         | other env. factor: temperature    | dummy   | 0=no other temp. measure, 1=other temp. measure controlled for
# "o_pre"         | other env. factor: precipitation  | dummy   | 0=no other prec. measure, 1=other prec. measure controlled for
# "o_sho"         | other env. factor: rapidonset     | dummy   | 0=no other rapid onset measure, 1=other rapid onset measure controlled for

# "internal",     |  SL considers internal migr. only | dummy   | 0=other migration, 1=internal migration 
# "dest_high"     | SL considers migr. to high-inc. c.| dummy   | 0=other migration, 1=international migration to high-income countries
# "dest_low",     | SL considers migr. to low-inc. c. | dummy   | 0=other migration, 1=international migration to low-income countries
# "dest_world"    | worldwide migr.                   | dummy   | 0=other migration, 1=international migration, worldwide
# "dest_ambi"     | destination cannot be categorized  | dummy   | 0=other migration, 1=international migration, destination ambiguous

# "period start"  | Start year of panel               | continuous
# "period end"    | End year of panel                 | continuous
# "income_channel"    | controls for income channel   | dummy   | 0=no control, 1= control
# "conflict_channel"  | controls for conflict channel | dummy   | 0=no control, 1= control 
# "weights"           | SL uses weights               | dummy   | 0=no weights, 1= weights 
# "robustness_check"  | SL is a robustness check      | dummy   | 0=no robustness check, 1= robustness check 
# "lin"               | SL uses a lin-lin estimation  | dummy   | 0=no lin-lin, 1= lin-lin 
# "gravity"           | SL uses a gravity model       | dummy   | 0=no gravity, 1= gravity 

# "published1"        | study was published           | dummy   | 0=unpublished (working paper, mimeo), 1= published 
# "published2"        | study was published by IF     | factor   | 0=unpublished, 1= published with IF <=2,2 = published with IF >2
# "pubimp"            | impact factor of publication  | continuous with 0 if not published by journal

# * Compositional shares *
# "nonoecd", "H","L" ,"LM" ,"UM", "M" ,"agr" "conflict_mepv_5" , "europena" ,"asia","lac" , "ssa" ,"mena"   
# > non oecd , high-income, low-income, lower-middle-income, upper-middle-income, agriculturally dependent, conflict, Europe/North America, Asia, Latin America/Caribbean, Sub Saharan Africa, Middle East and North Africa

# * Alternative conflict measures * 
# "conflict_mepv_5civ", "conflict_mepv_10", "conflict_mepv_10civ","conflict_ucdp_25_8" , "conflict_ucdp_0_8","conflict_ucdp_100_8", "conflict_ucdp_50_8"
 
# * Further controls included in specification *
# "controls_political" "controls_population" "controls_pastmigr"   "controls_econlevel"  "controls_culture" "controls_geo"   
#  "other_controls"  -> list of all control variables in specification 

# * Further study line model related information *
# "R2", "t.stats" ,  "se_detail" , "sample.size"  , "m_type", "m_detail"

# * Further information on migration outcome *
# "mig_destination", "mig_measure" , "mig_detail"  , "region"  

# * Further information on environmental factor/key independent variable *
# "env_type",  "env_location"  



##
## LOADING DATA AND PACKAGES *************************************
## 

rm(list=ls())

# >>> The following packages have to be installed for the code to run

require(meta)
require(metafor)
require(lfe)
require(lme4)
require(tidyverse)
require(stargazer)
require(ggpubr)
require(jtools)
require(RColorBrewer) 
require(rgeos)
require(rnaturalearth)

# get world map
world <- ne_countries(returnclass = "sf")
# set ggplot theme
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

citation("meta")
citation("metafor")
citation("lfe")
citation("lme4")
citation("tidyverse")
citation("stargazer")
citation("ggpubr")
citation("jtools")
citation("RColorBrewer") 
citation("rgeos")
citation("rnaturalearth")

# load data
load("Environmental Change and Migration_Meta and Country Data.rdata")

##  
## DESCRIPTIVE STATISTICS  ***************************************
##  

prop.table(table(meta$temperature));  table(meta$temperature)
prop.table(table(meta$precipitation));table(meta$precipitation)
prop.table(table(meta$rapidonset)); table(meta$rapidonset)
prop.table(table(meta$internal));table(meta$internal)


## ESTIMATING AVERAGE EFFECTS ACROSS ALL COEFFICIENTS Using Meta and Metafor packages
## see also: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/

m.pm <- metagen(stancoeff,
                stanse,
                data = meta,
                byvar=paper,
                studlab = paste(paper),
                comb.fixed = FALSE,
                comb.random = TRUE,
                method.tau = "PM",
                hakn = T,
                prediction = TRUE,
                sm = "COR")
m.pm


##
##  DEFINING VARIABLES FOR MODELING ******************************
##

## MODELS IN MAIN TEXT

# Baseline controls
con_base <- c("interaction", "fe_time", "fe_spatial", "control_sum") 

# Environmental variables 1
con_clim1 <- c("pre_var", "rapidonset", "tem_lev", "tem_var")

# Environmental variables 2
con_clim2 <- c("env_lag_dum", "env_timespan_dum", "env_other_dum")

# Migration variables
con_mig <- c("internal", "dest_low", "dest_high", "dest_ambi") # , 

# Compositional measures 1: Regional composition
con_context1 <- "nonoecd"

# Compositional measures 2: Composition by country characteristics
con_context2 <-c("L", "LM", "UM", "agr", "conflict_mepv_5")


## MODELS IN SUPPLEMENTARY MATERIALS

# control variables used in models
con_control <- c("income_channel", "conflict_channel", "controls_political", "controls_population", "controls_pastmigr", "controls_econlevel", "controls_culture", "controls_geo")

# Sample size of models
con_samplesize <- c("countrysample", "yearscovered")

# Time and spatial fixed effects
con_fixedeffects <- c("fe_time", "fe_spatial")

# Weighting and further estimation features
con_estimation <- c("weights", "lin", "robustness_check")

# Publication characteristics of studies
con_pub1 <- c( "published1")

# Additional controls for focus on differential time periods
con_time2 <- c("period_start" , "period_end")

# Defining data and weights
dat <- meta
wt <- 1/meta$stanse^2


## 
## FUNCTIONS TO IMPLEMENT MODELS *****************************
##  

func.felm <- function(yname,xnames, fixedeffect.cluster) {
  controls <- paste(xnames, collapse = "+")
  fe.c <- paste(fixedeffect.cluster,  collapse = "|" )
  formula <- as.formula(paste(yname,"~",controls, "|", fe.c, sep = ""))
  return(formula)
}

func.lmer <- function(yname,xnames) {
  controls <- paste(xnames, collapse = "+")
  formula <- as.formula(paste(yname,"~",controls,"+(1|paper)", sep = ""))
  return(formula)
}



#
##
### **************************************************************
### 2. OUTPUT TABLES AND FIGURES FOR MAIN PAPER ------------------
### **************************************************************
##
#


##
## FIGURE 1 - CODE AVAILABLE UPON REQUEST
##


## 
## TABLE 1 WEIGHTED FELM MODEL: STANDARDIZED COEFF AS OUTCOME -----------------
## 

## COLUMN 1
m1a <- func.felm("stancoeff", c(con_base, con_clim1), c("paper", "0", "paper"))
m1b <- felm(m1a , data = dat, weights=wt)
summary(m1b)

## COLUMN 2
m2a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2), c("paper", "0", "paper"))
m2b <- felm(m2a , data = dat, weights=wt)
summary(m2b)

## COLUMN 3
m3a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig), c("paper", "0", "paper"))
m3b <- felm(m3a , data = dat, weights=wt)
summary(m3b)

## COLUMN 4
m4a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig, con_context1), c("paper", "0" , "paper"))
m4b <- felm(m4a , data = dat, weights=wt) 
summary(m4b) 

## COLUMN 5
m5a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2), c("paper", "0" , "paper"))
m5b <- felm(m5a ,  data = dat, weights=wt)
summary(m5b)


## TABLE 1
stargazer(m1b, m2b, m3b, m4b, m5b,
          type="html",
          out="table 1_baseline_m1-m5.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = c("interaction", "fe_time", "fe_spatial", "control_sum", "yearscovered", "countrysample"))

rm(m1a, m2a, m3a, m4a, m5a, m1b, m2b, m3b, m4b, m5b)



##
## FIGURE 2 - VIOLIN PLOTS: DISTRIBUTION OF ENVIRONMENTAL EFFECTS BY ENVIRONMENTAL HAZARD -----
## 


meta %>% 
  mutate(
    climate_rec = as.factor(
      ifelse(rapidonset == 1,"Rapid-onset disaster",
             ifelse(pre_lev ==1, "Precipitation level change",
                    ifelse(pre_var == 1,"Precipitation variability / anomaly",
                           ifelse(tem_lev == 1, "Temperature level change",
                                  ifelse(tem_var == 1, "Temperature variability / anomaly", "NA")))))),
    climate_rec = fct_reorder(.f = climate_rec, .x = stancoeff.w, .fun = median, na.rm = TRUE, .desc = TRUE)
  ) %>% 
  ggplot(aes(y = climate_rec, x = stancoeff.w)) + 
  geom_vline(xintercept = 0) +
  geom_jitter(width = 0.3, size = 0.7, alpha = 0.3) +
  geom_violin(scale="width", bw = 0.5, fill = NA) + 
  geom_boxplot(outlier.shape=NA, width=0.32, fill="grey", alpha=0.5, position = position_dodge2(preserve = "single"))+
  labs(
    y = "",
    x = "Precision-weighted standardized environmental effect"
  ) +
  coord_cartesian(xlim = c(-5,7)) +
  scale_x_continuous(breaks = seq(-5,7, 1))+
  theme(axis.text.y=element_text(size=10))

ggsave("figure 2_violin_differences in environmental effects by hazards.jpg", dpi=600, width=8, height=6)



##
## FIGURE 3 - LINE PLOTS: PREDICTED ENVIRONMENTAL EFFECTS BY SAMPLE COMPOSITION -----
## 

## Estimating predicted values using a mixed random effects model

pred1a <- func.lmer("stancoeff", c(con_base, con_context2)) # including only compositional shares
pred1b <- lmer(pred1a, data = dat , weights=wt)
summary(pred1b)


pred1b_coef <- coefficients(pred1b)

dat.predict <- dat %>% 
  mutate(predict = mean(pred1b_coef$paper$`(Intercept)`)+
           mean(pred1b_coef$paper$fe_time)+ # predictions based on coefficients for models controling for time ...
           mean(pred1b_coef$paper$fe_spatial)+ # ... and spatial fixed effects
           L*mean(pred1b_coef$paper$L)+
           LM*mean(pred1b_coef$paper$LM)+
           UM*mean(pred1b_coef$paper$UM)+
           agr*mean(pred1b_coef$paper$agr)+
           conflict_mepv_5*mean(pred1b_coef$paper$conflict_mepv_5))

summary(dat.predict$predict)


## Plotting predicted effect sizes
plot1 <- dat.predict %>% 
  filter(agr<0.5 & conflict_mepv_5<0.5) %>% 
  ggplot(aes(y=predict, x=M))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.05, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color = "black") +
  labs(
    x = "% middle-income countries in sample", 
    y = "Predicted environmental effect"
    ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
    )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))
plot1

plot2 <- dat.predict %>%
  filter(L>0.8) %>% 
  ggplot(aes(y=predict, x=agr))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm",se=T, color="Black")+
  labs(
    x = "% agriculturally dependent countries in sample",
    y = ""
    )+
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
    )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.3),
        axis.text=element_text(size=11.3))
plot2

plot3 <- dat.predict %>%
  filter(agr<0.5) %>% 
  ggplot(aes(y=predict, x=conflict_mepv_5))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color="black")+
  labs(
    x = "% conflict countries in sample",
    y = ""
    ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
    ) +
  scale_y_continuous(breaks=seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))

  plot3

plot123 <- ggarrange(plot1, plot2, plot3, align="h", nrow = 1, labels = "auto")

plot123    
ggsave("figure 3_line plots_sample composition effects.jpg", width=15, height=5, dpi=600)

rm(pred1a,pred1b, plot1,plot2,plot3,plot123)


##
## FIGURE 4 - MAP: PREDICTED ENVIRONMENTAL EFFECT ON MIGRATION WORLDWIDE ------
## 

countrydata <- 
  countrydata %>%   
  mutate(
    predictedresponse = 
      mean(pred1b_coef$paper$`(Intercept)`)+
      mean(pred1b_coef$paper$fe_time)+ # predictions based on coefficients for models controling for time ...
      mean(pred1b_coef$paper$fe_spatial)+ #... and spatial fixed effects
      L*mean(pred1b_coef$paper$L)+
      LM*mean(pred1b_coef$paper$LM)+
      UM*mean(pred1b_coef$paper$UM)+
      agr*mean(pred1b_coef$paper$agr),
    conflict_mepv_5*mean(pred1b_coef$paper$conflict_mepv_5),
    predictedmig = env_change*predictedresponse,
    predictedmig_cat = cut(predictedmig, breaks = c(-Inf,-0.025,0.025,0.05,0.1,0.15,0.2,Inf))
  )

summary(countrydata$predictedmig)

map_world <- 
  world %>% 
  mutate(adm0_a3 = recode(adm0_a3, "SDS"="SSD")) %>% 
  right_join(countrydata, by = c("adm0_a3" = "iso3c"))

map_world %>% 
  ggplot() +
    geom_sf(aes(fill = predictedmig_cat)) + 
    coord_sf(crs = "+proj=eqearth") + 
    theme_minimal() +
    theme(
      axis.text = element_blank(),
    ) +
    scale_fill_manual(
      name = "Predicted migration", 
      values = c(
        "(-Inf,-0.025]" = "#3f7ad9",
        "(-0.025,0.025]" = "#dee4fa",
        "(0.025,0.05]"="#fae0de",
        "(0.05,0.1]"="#de9999",
        "(0.1,0.15]"="#db6969",
        "(0.15,0.2]"="#de4949",
        "(0.2, Inf]"="#db1414"
      ),
      labels = c(
        "negative [-0.55,-0.025]",
        "none (-0.025,0.025]",
        "very low (0.025,0.05]",
        "low (0.05,0.1]",
        "moderate (0.1,0.15]",
        "high (0.15, 0.20]",
        "very high (0.2, 0.50]"),
      na.translate = FALSE
    )
ggsave("figure 4_map_predicted environmental migration.jpg", width = 10, height=4.5, dpi=600)




#
##
### **************************************************************
### 3. OUTPUT TABLES AND FIGURES FOR EXTENDED DATA ---------------
### **************************************************************
### 
##
#


###
### EXTENDED DATA FIGURE 1 - CODE AVAILABLE UPON REQUEST
### 

###
### EXTENDED DATA FIGURE 2 - Density Plots: Distribution of Coefficients and Between and Within Study Variation ----
### 


dens <- with(
  density(dat$stancoeff, bw = 0.07, n = 1e4), 
  data.frame(x, y))
density <- 
  dens %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_line() +
  geom_ribbon(
    data = dens %>% filter(x <= 0),
    aes(x = x, ymin = 0, ymax = y),
    fill="grey"
  ) +
  geom_ribbon(
    data = dens %>% filter(x >= 0),
    aes(x = x, ymin = 0, ymax = y),
    fill="grey20"
  ) +
  geom_vline(xintercept = 0) +
  coord_cartesian(
    xlim = c(-1, 1),
    ylim = c(0, 4.1)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-1, 1, 0.25)) +
  labs(
    x = "Standardized environmental effects of all model estimates",
    y = "Density"
  )
density


meta %>% 
  group_by(paper) %>% 
  summarize(stancoeff.mean = mean(stancoeff)) %>% 
  summarize(var=var(stancoeff.mean), sd=sqrt(var))

between <- meta %>% 
  group_by(paper) %>% 
  summarize(stancoeff.mean = mean(stancoeff)) %>% 
  ggplot(aes(x=stancoeff.mean))+
  geom_density(fill="grey", color="grey", adjust = 2)+
  geom_vline(xintercept = 0) +
  xlab("Mean study effect")+
  ylab("Density")+
  coord_cartesian(xlim = c(-0.3,0.5))
between

within <- meta %>% 
  group_by(paper) %>% 
  mutate(
    stancoeff.mean = mean(stancoeff), 
    stancoeff.dev = stancoeff - stancoeff.mean
  ) %>% 
  ggplot(aes(x=stancoeff.dev))+
  geom_density(fill="grey", color="grey")+
  geom_vline(xintercept = 0) +
  xlab("Deviations of estimates from mean study effect")+
  ylab("Density")+
  coord_cartesian(xlim = c(-1,1))
within

plot1 <- ggarrange(density, align = "h", labels = c("a"))
plot1

plot2 <-ggarrange(between, within, align = "h", labels = c("b", "c"))
plot2

ggarrange(plot1, plot2, nrow=2, align = "h")

ggsave("figure ED2_density_distribution of coefficients, between and within study variation.jpg", height=4, width = 7.5, dpi=600)



###
### EXTENDED DATA FIGURE 3 - BOXPLOT AND MAP: DIFFERENCES IN SAMPLE COMPOSITION BY STUDY LINES ----
###

plot_comp_box <- 
  meta %>% 
  gather(
    "conflict_mepv_5", "agr", "UM", "LM", "L", "nonoecd", "lac","mena", "ssa", "asia", "europena", 
    key = "group", value="fraction", convert=T
  ) %>% 
  ggplot(aes(y = fct_inorder(group, ordered=TRUE), x = fraction))+
  geom_boxplot()+
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_discrete(
    labels = c(
      'Conflict >= 5 years',
      'Agriculturally dependent', 
      'Upper-Middle-Income',
      'Lower-Middle-Income',
      'Low-Income', 
      'Non-OECD',
      'LAC', 
      'MENA',
      'SSA',
      'Asia', 
      'Europe/North America'
    )
  ) +
  labs(
    x = "% of countries in study-line samples",
    y = ""
  )
plot_comp_box
ggsave("figure ED3a_boxplot_sample composition measures.jpg", plot = plot_comp_box, height=5, width=6.5)

plot_sample_map <- 
  world %>% 
  left_join(country_rep, by = c("adm0_a3" = "iso3c")) %>% 
  filter(sovereignt != "Antarctica") %>% 
  ggplot() + 
  geom_sf(aes(fill = n)) +
  coord_sf(crs = "+proj=eqearth") +
  theme_minimal() +
  theme(
    axis.text = element_blank()
  ) +
  scale_fill_gradient(
    name = paste0("Number of\nSamples"), #
    low = "white",
    high = "#3f7ad9", 
    na.value = "grey50" 
  )
plot_sample_map
ggsave("figure ED3b_map_sample composition.jpg", plot = plot_sample_map, height=5, width=10)



 
#
##
### **************************************************************
### 4. OUTPUT TABLES AND FIGURES FOR SUPPLEMENTARY MATERIAL ------
### **************************************************************
##
#

###
### FIGURE S1 PRISMA CHART
### 

###
### FIGURE S2- BOXPLOTS: DISTRIBUTION OF STANDARDIZED COEFFICIENTS ----
### 


meta %>% 
  ggplot() +
  geom_vline(aes(xintercept = 0), color="black") +
  geom_vline(aes(xintercept = median(stancoeff)), color="grey", linetype="dashed", size =0.8) +
  geom_boxplot(aes(y = fct_reorder(.f = paper, .x = stancoeff, .fun = median, .desc = TRUE), x = stancoeff),
               outlier.shape = NA) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_x_continuous(breaks = seq(-1, 1, 0.25)) +
  labs(
    x = "Standardized environmental effect",
    y = ""
  )+
  theme(axis.text = element_text(size=10))
ggsave("figure s2_boxplot_distribution of standardized effects.svg", width = 7, height=7, dpi=600)


###
### FIGURE S3 - CONCEPTUAL GRAPH
### 

###
### TABLE S1 and S2 - MANUALLY CREATED
###  

###
### TABLE S3 - DESCRIPTIVES ---------------------------------------------------
### 


meta.sum <- meta %>% 
  mutate(robustness_check = as.numeric(robustness_check)) %>% 
  select(
    stancoeff,  stanse, stancoeff.w,
    pre_lev, pre_var, rapidonset, tem_lev, tem_var,
    env_lag_dum, env_other_dum,env_timespan_dum,
    o_tem, o_pre, o_sho,
    internal, dest_world, dest_low, dest_high,dest_ambi,
    nonoecd, L, LM,UM, agr, conflict_mepv_5,
    europena, ssa, mena, lac,asia,
    income_channel, conflict_channel, 
    controls_political, controls_population, controls_pastmigr, controls_econlevel, controls_culture, controls_geo, control_sum,
    countrysample, yearscovered,period_start, period_end,
    fe_spatial, fe_time,
    weights, lin,robustness_check,
    published1
  ) 

  stargazer(
    as.data.frame(meta.sum),
    summary = TRUE, 
    summary.stat = c("mean", "sd", "min", "median", "max"),
    type = "html", 
    out = "table s3_summary statistics.doc"
  )

rm(meta.sum)


###  
### TABLE S4 - ROBUSTNESS: KEEPING ONLY MODELS WITH SPATIAL FE  ---------------
###

aux2 <- meta %>% filter(fe_spatial==1)

## COLUMN 1
r1a <- func.felm("stancoeff", c(con_base, con_clim1 ), c("paper", "0", "paper"))
r1b <- felm(r1a , data = aux2, weights=1/aux2$stanse^2)
summary(r1b)

## COLUMN 2
r2a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2), c("paper", "0", "paper"))
r2b <- felm(r2a , data = aux2, weights=1/aux2$stanse^2)
summary(r2b)

## COLUMN 3
r3a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig), c("paper", "0", "paper"))
r3b <- felm(r3a , data = aux2, weights=1/aux2$stanse^2)
summary(r3b)

## COLUMN 4
r4a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context1), c("paper", "0" , "paper"))
r4b <- felm(r4a , data = aux2, weights=1/aux2$stanse^2) 
summary(r4b) 

## COLUMN 5
r5a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2), c("paper", "0" , "paper"))
r5b <- felm(r5a ,  data = aux2, weights=1/aux2$stanse^2)
summary(r5b)


stargazer(r1b, r2b, r3b, r4b, r5b,
          type="html",
          out="table s4_robustness_only models with spatial fe.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)

###
### TABLE S5 - ROBUSTNESS: KEEPING ONLY MODELS WITH TIME FE  ------------------
###

aux3 <- meta %>% filter(fe_time==1)


## COLUMN 1
r1a <- func.felm("stancoeff", c(con_base, con_clim1 ), c("paper", "0", "paper"))
r1b <- felm(r1a , data = aux3, weights=1/aux3$stanse^2)
summary(r1b)

## COLUMN 2
r2a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2), c("paper", "0", "paper"))
r2b <- felm(r2a , data = aux3, weights=1/aux3$stanse^2)
summary(r2b)

## COLUMN 3
r3a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig), c("paper", "0", "paper"))
r3b <- felm(r3a , data = aux3, weights=1/aux3$stanse^2)
summary(r3b)

## COLUMN 4
r4a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context1), c("paper", "0" , "paper"))
r4b <- felm(r4a , data = aux3, weights=1/aux3$stanse^2) 
summary(r4b) 

## COLUMN 5
r5a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2), c("paper", "0" , "paper"))
r5b <- felm(r5a ,  data = aux3, weights=1/aux3$stanse^2)
summary(r5b)


stargazer(r1b, r2b, r3b, r4b, r5b,
          type="html",
          out="table s5_robustness_only models with time fe.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)


### 
### TABLE S6 - ROBUSTNESS: LMER MODELS ----------------------------------------
###


## COLUMN 1
r1a <- func.lmer("stancoeff", c(con_base, con_clim1 ))
r1b <- lmer(r1a , data = dat, weights=wt)
summary(r1b)

## COLUMN 2
r2a <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2))
r2b <- lmer(r2a , data = dat, weights=wt)
summary(r2b)

## COLUMN 3
r3a <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig))
r3b <- lmer(r3a , data = dat, weights=wt)
summary(r3b)

## COLUMN 4
r4a <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context1))
r4b <- lmer(r4a , data = dat, weights=wt) 
summary(r4b) 

## COLUMN 5
r5a <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2))
r5b <- lmer(r5a ,  data = dat, weights=wt)
summary(r5b)


stargazer(r1b, r2b, r3b, r4b, r5b,
          type="html",
          out="table s6_robustness_random effects-mixed effects models.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)


###  
### TABLE S7 - ROBUSTNESS: ADDING ALL CONTROLS IN ONE MODEL -------------------
### 

## COLUMN 1
r1a <- func.felm("stancoeff", c(con_base,con_control, con_samplesize, con_estimation, con_time2, con_clim1 ), c("paper", "0", "paper"))
r1b <- felm(r1a , data = dat, weights=wt)
summary(r1b)

## COLUMN 2
r2a <- func.felm("stancoeff", c(con_base,con_control, con_samplesize, con_estimation, con_time2,  con_clim1, con_clim2), c("paper", "0", "paper"))
r2b <- felm(r2a , data = dat, weights=wt)
summary(r2b)

## COLUMN 3
r3a <- func.felm("stancoeff", c(con_base,con_control, con_samplesize, con_estimation, con_time2,  con_clim1,  con_clim1, con_clim2, con_mig), c("paper", "0", "paper"))
r3b <- felm(r3a , data = dat, weights=wt)
summary(r3b)

## COLUMN 4
r4a <- func.felm("stancoeff", c(con_base,con_control, con_samplesize, con_estimation, con_time2,  con_clim1,  con_clim1, con_clim2, con_mig, con_context1), c("paper", "0" , "paper"))
r4b <- felm(r4a , data = dat, weights=wt) 
summary(r4b) 

## COLUMN 5
r5a <- func.felm("stancoeff", c(con_base,con_control, con_samplesize, con_estimation, con_time2,  con_clim1,  con_clim1, con_clim2, con_mig, con_context2), c("paper", "0" , "paper"))
r5b <- felm(r5a ,  data = dat, weights=wt)
summary(r5b)


stargazer(r1b, r2b, r3b, r4b, r5b,
          type="html",
          out="table s7_robustness_including all control variables.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = c(con_base,con_control, con_samplesize, con_estimation, con_time2))



###  
### TABLE S8 - ROBUSTNESS: NO WEIGHTING & REMOVE CASES STANSE > 2  ------------
### 

aux4 <- meta %>% filter(abs(stancoeff) < 2)

## COLUMN 1
r1a <- func.felm("stancoeff", c(con_base, con_clim1 ), c("paper", "0", "paper"))
r1b <- felm(r1a , data = aux4, weights=1/aux4$stanse^2)
summary(r1b)

## COLUMN 2
r2a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2), c("paper", "0", "paper"))
r2b <- felm(r2a , data = aux4, weights=1/aux4$stanse^2)
summary(r2b)

## COLUMN 3
r3a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig), c("paper", "0", "paper"))
r3b <- felm(r3a , data = aux4, weights=1/aux4$stanse^2)
summary(r3b)

## COLUMN 4
r4a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context1), c("paper", "0" , "paper"))
r4b <- felm(r4a , data = aux4, weights=1/aux4$stanse^2) 
summary(r4b) 

## COLUMN 5
r5a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2), c("paper", "0" , "paper"))
r5b <- felm(r5a ,  data = aux4, weights=1/aux4$stanse^2)
summary(r5b)


stargazer(r1b, r2b, r3b, r4b, r5b,
          type="html",
          out="table s8_robustness_removing study lines with large effects 2sd.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)


###  
### TABLE S9 - ROBUSTNESS: EXCLUDING INTERNAL MIGRATION -----------------------
###

aux5 <- meta %>% filter(internal==0)

## COLUMN 1
r1a <- func.felm("stancoeff", c(con_base, con_clim1 ), c("paper", "0", "paper"))
r1b <- felm(r1a , data = aux5, weights=1/aux5$stanse^2)
summary(r1b)

## COLUMN 2
r2a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2), c("paper", "0", "paper"))
r2b <- felm(r2a , data = aux5, weights=1/aux5$stanse^2)
summary(r2b)

## COLUMN 3
r3a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig), c("paper", "0", "paper"))
r3b <- felm(r3a , data = aux5, weights=1/aux5$stanse^2)
summary(r3b)

## COLUMN 4
r4a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context1), c("paper", "0" , "paper"))
r4b <- felm(r4a , data = aux5, weights=1/aux5$stanse^2) 
summary(r4b) 

## COLUMN 5
r5a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2, "conflict_mepv_5"), c("paper", "0" , "paper"))
r5b <- felm(r5a ,  data = aux5, weights=1/aux5$stanse^2)
summary(r5b)


stargazer(r1b, r2b, r3b, r4b, r5b,
          type="html",
          out="table s9_robustness_only international migration.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)

###  
### TABLE S10 - ROBUSTNESS: ALTERNATIVE CONFLICT MEASURES  --------------------
###


con_context2 <-c("L", "LM", "UM", "agr")

r1a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig, con_context2, "conflict_mepv_5"), c("paper", "0" , "paper"))
r1b <- felm(r1a , data = dat, weights=wt)
summary(r1b)

r2a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig, con_context2,  "conflict_mepv_5civ"), c("paper", "0" , "paper"))
r2b <- felm(r2a , data = dat, weights=wt)
summary(r2b)

r3a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig, con_context2, "conflict_mepv_10"), c("paper", "0" , "paper"))
r3b <- felm(r3a , data = dat, weights=wt)
summary(r3b)

r4a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig, con_context2, "conflict_ucdp_0_8"), c("paper", "0" , "paper"))
r4b <- felm(r4a , data = dat, weights=wt)
summary(r4b)

r5a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig, con_context2, "conflict_ucdp_25_8"), c("paper", "0" , "paper"))
r5b <- felm(r5a , data = dat, weights=wt)
summary(r5b)

stargazer(r1b, r2b, r3b, r4b, r5b,
          type="html",
          out="table s10_robustness_alternative conflict measures.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)



### 
### TABLE S11 - EXTENDED FELM MODEL: DIFFERENCES BY CONTROL FOR OTHER ENV. VARIABLES -----
###


aux1 <- meta %>% filter(rapidonset == 0 & temperature==0)
clim1a <- func.felm("stancoeff", c(con_base, "o_pre","o_tem", "o_sho"), c("paper", "0", "paper"))
clim1b <- felm(clim1a , data = aux1, weights=1/aux1$stanse^2)
summary(clim1b)

clim2a <- func.felm("stancoeff.abs", c(con_base, "o_pre","o_tem", "o_sho"), c("paper", "0", "paper"))
clim2b <- felm(clim2a , data = aux1, weights=1/aux1$stanse^2)
summary(clim2b)

aux2 <- meta %>% filter(rapidonset == 0 & precipitation==0)
clim3a <- func.felm("stancoeff", c(con_base, "o_pre","o_tem", "o_sho"), c("paper", "0", "paper"))
clim3b <- felm(clim3a , data = aux2, weights=1/aux2$stanse^2)
summary(clim3b)

clim4a <- func.felm("stancoeff.abs", c(con_base, "o_pre","o_tem", "o_sho"), c("paper", "0", "paper"))
clim4b <- felm(clim4a , data = aux2, weights=1/aux2$stanse^2)
summary(clim4b)

aux3 <- meta %>% filter(temperature == 0 & precipitation==0)
clim5a <- func.felm("stancoeff", c(con_base, "o_pre", "o_tem",  "o_sho"), c("paper", "0", "paper"))
clim5b <- felm(clim5a , data = aux3, weights=1/aux3$stanse^2)
summary(clim5b)

clim6a <- func.felm("stancoeff.abs", c(con_base, "o_pre","o_tem", "o_sho"), c("paper", "0", "paper"))
clim6b <- felm(clim6a , data = aux3, weights=1/aux3$stanse^2)
summary(clim6b)

stargazer(clim1b, clim3b, clim5b, clim2b, clim4b,clim6b,
          type="html",
          out="table s11_extended_testing for interactions between env. variables.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)

### 
### TABLE S12 - EXTENDED FELM MODEL: TESTING FOR PUSH AND PULL FACTORS --------
###

m1a <- func.felm("stancoeff", c(con_base, "env_location"), c("paper", "0" , "paper"))
m1b <- felm(m1a , data = meta.plus, weights=1/meta.plus$stanse^2)
summary(m1b)

m2a <- func.felm("stancoeff", c(con_base, "env_location", con_clim1, con_clim2, con_mig), c("paper", "0" , "paper"))
m2b <- felm(m2a , data = meta.plus, weights=1/meta.plus$stanse^2)
summary(m2b)


stargazer(m1b, m2b, 
          type="html",
          out="table s12_extended_environment as push and pull.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)

# On average, environmental changes in the destination region reduce migraion
aux4 <- meta.plus %>% filter(env_location=="destination")
summary(lm(stancoeff~1, data=aux4, weights = 1/aux4$stanse^2))


## 
## TABLE S13 - EXTENDED FELM MODELS: TESTING FOR INFLUENCE OF DIFFERENT REGIONAL SAMPLE COMPOSITIONS ----
## 

## COLUMN 1
m1a <- func.felm("stancoeff", c(con_base,con_clim1, con_clim2, con_mig, "europena"), c("paper", "0", "paper"))
m1b <- felm(m1a , data = dat, weights=wt)
summary(m1b)

## COLUMN 2
m2a <- func.felm("stancoeff", c(con_base,con_clim1, con_clim2, con_mig, "ssa"), c("paper", "0", "paper"))
m2b <- felm(m2a , data = dat, weights=wt)
summary(m2b)

## COLUMN 3
m3a <- func.felm("stancoeff", c(con_base,con_clim1, con_clim2, con_mig, "mena"), c("paper", "0", "paper"))
m3b <- felm(m3a , data = dat, weights=wt)
summary(m3b)

## COLUMN 4
m4a <- func.felm("stancoeff", c(con_base,con_clim1, con_clim2, con_mig,  "lac"), c("paper", "0" , "paper"))
m4b <- felm(m4a , data = dat, weights=wt) 
summary(m4b) 

## COLUMN 5
m5a <- func.felm("stancoeff", c(con_base,con_clim1, con_clim2, con_mig, "asia"), c("paper", "0" , "paper"))
m5b <- felm(m5a ,  data = dat, weights=wt)
summary(m5b)


stargazer(m1b, m2b, m3b, m4b, m5b,
          type="html",
          out="table s13_extended_testing for influence of different regional sample compositions.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base )
          


###  
### TABLE S14 - EXTENDED FELM MODEL: ABSOLUTE EFFECT SIZES AS OUTCOME ---------
###

abs1a <- func.felm("stancoeff.abs", c("interaction", con_control, "control_sum"), c("paper", "0", "paper"))
abs1b <- felm(abs1a , data = dat, weights=wt)
summary(abs1b)

abs2a <- func.felm("stancoeff.abs", c("interaction", con_fixedeffects ), c("paper", "0", "paper"))
abs2b <- felm(abs2a , data = dat, weights=wt)
summary(abs2b)

abs3a <- func.felm("stancoeff.abs", c("interaction", con_samplesize ), c("paper", "0", "paper"))
abs3b <- felm(abs3a , data = dat, weights=wt)
summary(abs3b)

abs4a <- func.felm("stancoeff.abs", c("interaction", con_estimation), c("paper", "0", "paper"))
abs4b <- felm(abs4a , data = dat, weights=wt)
summary(abs4b)

abs5a <- func.felm("stancoeff", c("interaction", con_time2), c("paper", "0", "paper"))
abs5b <- felm(abs5a , data = dat, weights=wt)
summary(abs5b)

stargazer(abs1b, abs2b, abs3b, abs4b,abs5b,
          type="html",
          out="table s14_extended_absolute coefficients as outcome.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = c("interaction"))



###
### TABLE S15 - EXTENDED LMER MODELS: SIMPLIFIED MAIN MODEL FOR PREDICTIONS ----
###

con_context2 <-c("L", "LM", "UM", "agr", "conflict_mepv_5")

pred1a <- func.lmer("stancoeff", c(con_base, con_context2)) # including only compositional shares
pred1b <- lmer(pred1a, data = dat , weights=wt)
summary(pred1b)


stargazer(pred1b,
          type="html",
          out="table s15_extended_simplified main model for predictions.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)


###
### FIGURE S4 - MAP: ENVIRONMENTAL CHANGE WORLDWIDE ---------------------------
###

map_world %>% 
  ggplot() +
  geom_sf(aes(fill=env_change))+
  coord_sf(crs = "+proj=eqearth") + 
  theme_minimal() +
  theme(
    axis.text = element_blank(),
  ) +
  scale_fill_gradient(
    name = "Environmental change",
    low="white", 
    high="red", 
    breaks = 1:6, 
    na.value = "grey50" 
    )

ggsave("figure s4_map_observed environmental change.jpg", dpi=600)


###
### FIGURE S5 - FUNNEL PLOT: TESTING FOR PUBLICATION EFFECTS ------------------
### 

### THE CODE FOR THE FUNNEL PLOT WAS CREATED BY John K. Sakaluk (2016) : https://sakaluk.wordpress.com/2016/02/16/7-make-it-pretty-plots-for-meta-analysis/

summary(meta$stancoeff)
summary(meta$stanse)

#Store the meta-analytic estimate and its standard error from whatever model you run (substitute your own values)
estimate = median(meta$stancoeff) # median values
se = median(meta$stanse)

#Store a vector of values that spans the range from 0
#to the max value of impression (standard error) in your dataset.
#Make the increment (the final value) small enough (I choose 0.001)
#to ensure your whole range of data is captured
se.seq=seq(0, max(meta$stanse), 0.001)

#Now, compute vectors of the lower-limit and upper limit values for
#the 95% CI region, using the range of SE that you generated in the previous step, and the stored value of your meta-analytic estimate.
ll95 = estimate-(1.96*se.seq)
ul95 = estimate+(1.96*se.seq)

#You can do this for a 99% CI region too
ll99 = estimate-(3.29*se.seq)
ul99 = estimate+(3.29*se.seq)

#And finally, do the same thing except now calculating the confidence interval
#for your meta-analytic estimate based on the stored value of its standard error
meanll95 = estimate-(1.96*se)
meanul95 = estimate+(1.96*se)

#Now, smash all of those calculated values into one data frame (called 'dfCI').
#You might get a warning about '...row names were found from a short variable...'
#You can ignore it.
dfCI = data.frame(ll95, ul95, ll99, ul99, se.seq, estimate, meanll95, meanul95)


 
fp <- meta %>% 
  filter(abs(stancoeff)<=10) %>% 
  ggplot(aes(y = stanse, x = stancoeff)) +
  geom_point(shape = 1, size=1.5) +
  ylab('Standard Error') + xlab('Standardized Coefficient')+
  geom_line(aes(y = se.seq, x = ll95), linetype = 'dotted', data = dfCI) +
  geom_line(aes(y = se.seq, x = ul95), linetype = 'dotted', data = dfCI) +
   geom_segment(aes(y = min(se.seq), x = meanll95, yend = max(se.seq), xend = meanll95), linetype='dotted', data=dfCI) +
  geom_segment(aes(y = min(se.seq), x = meanul95, yend = max(se.seq), xend = meanul95), linetype='dotted', data=dfCI) +
  coord_cartesian(xlim = c(-5,5), ylim = c(2.1,0))+
  scale_x_continuous()+
  theme_bw()+
  scale_y_reverse(limits=c(790,0), expand=c(0,0))+
  theme(axis.text = element_text(size=11))
  

fp

ggsave('figure s5_funnel plot_Standard errors against coefficient size.jpg', width=8, height=6,
       dpi=600)


### 
### Table S16 - EXTENDED MODEL: TESTING FOR PUBLICATION EFFECTS ---------------
###

pu1a <- func.lmer("stancoeff", c("published1"))
pu1b <- lmer(pu1a , data = meta, weights=wt)

pu2a <- func.lmer("stancoeff.abs", c("published1"))
pu2b <- lmer(pu2a , data = meta, weights=wt)

pu3a <- func.lmer("stancoeff", c("published1", con_base))
pu3b <- lmer(pu3a , data = dat, weights=wt)

pu4a <- func.lmer("stancoeff.abs", c("published1", con_base))
pu4b <- lmer(pu4a ,  data = dat, weights=wt)

 

stargazer(pu1b, pu2b, pu3b, pu4b, 
          type="html",
          out="table s16_extended_testing for publication effects.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = con_base)







 