library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)
library(drjacoby)
library(matrixStats)
library(purrr)
library(zoo)

setwd(dirname(getActiveDocumentContext()$path)) 


################################################################################
### TABLE S4 = Table 2 for the alternative severity model ####################################################################
################################################################################

source("R/vx_profile_1.R")
load("../Model Fitting/Sensitivity - Alternative Severity/UKHSA_v6pn_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=TRUE_mcmc_chain.Rdata")  


chain <- mcmc$output %>%
  filter(phase == "sampling") %>%
  subset(select=-c(chain,phase,iteration,logprior,loglikelihood))

posterior_median <- chain %>%
  summarise( 
    across(where(is.numeric), median)
  ) %>%
  mutate(sample=0)
draws <- sample_chains(mcmc,1000)


df_raw <- bind_rows(posterior_median, draws) %>%
          mutate(AZ_ns_off = 0)     # turn off the AZ_ns_off parameter 

## calculate combinations but since the dose 3 is independent of dose 2 in the main model, can just use three for the table 

df_AZ_AZ <- df_raw %>%
  mutate(vaccine = "AZ_AZ",
         d1 = d1_AZ,
         d2 = d2_AZ,
         d3 = bst_AZ,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df_PF_PF <- df_raw %>%
  mutate(vaccine = "PF_PF",
         d1 = d1_PF,
         d2 = d2_PF,
         d3 = bst_PF,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df_MD_MD <- df_raw %>%
  mutate(vaccine = "MD_MD",
         d1 = d1_MD,
         d2 = d2_MD,
         d3 = bst_MD,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df <- bind_rows(df_AZ_AZ, df_PF_PF, df_MD_MD) 

r1 <- 
  # Create input options
  expand_grid(
  vfr = c(0,1),
  vaccine = c("AZ_AZ", "PF_PF", "MD_MD"),
  sample = unique(df$sample)) %>%
  # Join with MCMC samples
  left_join(df, by = c("vaccine", "sample")) %>%
  # Apply vaccine profile function to each row
  mutate(profile = pmap(., vx_profile_altsev)) %>%
  # Format
  unnest(cols = c(profile)) %>%
  pivot_longer(cols = c(Titre_d1, Titre_d2, Titre_d3, Efficacy_d1, Efficacy_d2, Efficacy_d3, Severe_Efficacy_d1,Severe_Efficacy_d2, Severe_Efficacy_d3, Death_Efficacy_d1, Death_Efficacy_d2,Death_Efficacy_d3), names_to = "group", values_to = "Value") %>% 
  mutate(dose = case_when(group == "Titre_d1" | group == "Efficacy_d1" | group == "Severe_Efficacy_d1" | group == "Death_Efficacy_d1"~ 1,
                          group == "Titre_d2" | group == "Efficacy_d2" | group == "Severe_Efficacy_d2" | group == "Death_Efficacy_d2" ~ 2,
                          group == "Titre_d3" | group == "Efficacy_d3" | group == "Severe_Efficacy_d3" | group == "Death_Efficacy_d3"~ 3 )) %>%
  mutate(group = substr(group, 1, nchar(group) - 3))



# posterior median
r1a <- r1 %>%
  filter(sample==0) %>%
  mutate(median = Value)

# mcmc samples
r1b <- r1 %>%
  filter(sample>0) %>%
  filter(dose>1) %>%
  group_by(vfr,vaccine,group,dose,t) %>%
  summarise(upper = quantile(Value, 0.975),
            lower = quantile(Value, 0.025),
  ) 

summary_stats <- left_join(r1a,r1b)


table1 <- summary_stats %>%
  filter(group=="Efficacy"| group=="Severe_Efficacy" | group=="Death_Efficacy") %>%
  mutate(median = round(median,digits=1), lower = round(lower,digits=1), upper = round(upper,digits=1)) %>%
  mutate(vaccine_efficacy = paste0(median, " (",lower,"-",upper,")")) 

tb_dose_2 <- table1 %>%
  filter(dose==2) %>%
  filter(t==1+90 | t==1+180 ) %>%
  mutate(period = case_when(t==91 ~ "90d pd2",
                            t==181 ~ "180d pd2") )

tb_dose_3 <- table1 %>%
  filter(dose==3) %>% 
  filter(t==1+30 | t==1+60  | t==1+90 | t==1+120 | t==1+150 | t==1+180 | t==1+365) %>%
  mutate(period = case_when(t==31 ~ "30d pb",
                            t==61 ~ "60d pb",
                            t==91 ~ "90d pb",
                            t==121 ~ "120d pb",
                            t==151 ~ "150d pb",
                            t==181 ~ "180d pb",
                            t==366 ~ "365d pb") )

tb <- rbind(tb_dose_2,tb_dose_3)
tb <- subset(tb,select=c(vfr,vaccine,group,vaccine_efficacy,period))


tb.inf <- tb %>%
  filter(group == "Efficacy") %>%
  pivot_wider(id_cols = c(vaccine, vfr), names_from = period, values_from = vaccine_efficacy) %>%
  arrange(vaccine, vfr) %>%
  rename(Vaccine = vaccine, "0 = delta, 1 = omicron" = vfr) %>% 
  flextable()

tb.sev <- tb %>%
  filter(group == "Severe_Efficacy") %>%
  pivot_wider(id_cols = c(vaccine, vfr), names_from = period, values_from = vaccine_efficacy) %>%
  arrange(vaccine, vfr) %>%
  rename(Vaccine = vaccine, "0 = delta, 1 = omicron" = vfr) %>% 
  flextable()

tb.death <- tb %>%
  filter(group == "Death_Efficacy") %>%
  pivot_wider(id_cols = c(vaccine, vfr), names_from = period, values_from = vaccine_efficacy) %>%
  arrange(vaccine, vfr) %>%
  rename(Vaccine = vaccine, "0 = delta, 1 = omicron" = vfr) %>% 
  flextable()

save_as_docx("Efficacy against infection" = tb.inf,
             "Efficacy against severe disease" = tb.sev,
             "Efficacy against death" = tb.death,
             path = "../Tables/TableS4.docx")

################################################################################
### TABLE S5 = Table 2 for the additive model ####################################################################
################################################################################

source("R/vx_profile.R")
load("../Model Fitting/Sensitivity - Additive Boost/UKHSA_v6pn_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=TRUE_AltSev=FALSE_mcmc_chain.Rdata")  


chain <- mcmc$output %>%
  filter(phase == "sampling") %>%
  subset(select=-c(chain,phase,iteration,logprior,loglikelihood))

posterior_median <- chain %>%
  summarise( 
    across(where(is.numeric), median)
  ) %>%
  mutate(sample=0)
draws <- sample_chains(mcmc,1000)


df_raw <- bind_rows(posterior_median, draws) %>%
  mutate(AZ_ns_off = 0)     # turn off the AZ_ns_off parameter 

## calculate combinations but since the dose 3 is independent of dose 2 in the main model, can just use three for the table 

df_AZ_AZ <- df_raw %>%
  mutate(vaccine = "AZ_AZ",
         d1 = d2_AZ + d1_AZ,
         d2 = d2_AZ,
         d3 = d2_AZ + bst_AZ,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df_PF_PF <- df_raw %>%
  mutate(vaccine = "PF_PF",
         d1 = d2_PF + d1_PF,
         d2 = d2_PF,
         d3 = d2_PF + bst_PF,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df_MD_MD <- df_raw %>%
  mutate(vaccine = "MD_MD",
         d1 = d2_MD + d1_MD,
         d2 = d2_MD,
         d3 = d2_MD + bst_MD,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df <- bind_rows(df_AZ_AZ, df_PF_PF, df_MD_MD) 

r1 <- 
  # Create input options
  expand_grid(
    vfr = c(0,1),
    vaccine = c("AZ_AZ", "PF_PF", "MD_MD"),
    sample = unique(df$sample)) %>%
  # Join with MCMC samples
  left_join(df, by = c("vaccine", "sample")) %>%
  # Apply vaccine profile function to each row
  mutate(profile = pmap(., vx_profile)) %>%
  # Format
  unnest(cols = c(profile)) %>%
  pivot_longer(cols = c(Titre_d1, Titre_d2, Titre_d3, Efficacy_d1, Efficacy_d2, Efficacy_d3, Severe_Efficacy_d1,Severe_Efficacy_d2, Severe_Efficacy_d3, Death_Efficacy_d1, Death_Efficacy_d2,Death_Efficacy_d3), names_to = "group", values_to = "Value") %>% 
  mutate(dose = case_when(group == "Titre_d1" | group == "Efficacy_d1" | group == "Severe_Efficacy_d1" | group == "Death_Efficacy_d1"~ 1,
                          group == "Titre_d2" | group == "Efficacy_d2" | group == "Severe_Efficacy_d2" | group == "Death_Efficacy_d2" ~ 2,
                          group == "Titre_d3" | group == "Efficacy_d3" | group == "Severe_Efficacy_d3" | group == "Death_Efficacy_d3"~ 3 )) %>%
  mutate(group = substr(group, 1, nchar(group) - 3))



# posterior median
r1a <- r1 %>%
  filter(sample==0) %>%
  mutate(median = Value)

# mcmc samples
r1b <- r1 %>%
  filter(sample>0) %>%
  filter(dose>1) %>%
  group_by(vfr,vaccine,group,dose,t) %>%
  summarise(upper = quantile(Value, 0.975),
            lower = quantile(Value, 0.025),
  ) 

summary_stats <- left_join(r1a,r1b)


table1 <- summary_stats %>%
  filter(group=="Efficacy"| group=="Severe_Efficacy" | group=="Death_Efficacy") %>%
  mutate(median = round(median,digits=1), lower = round(lower,digits=1), upper = round(upper,digits=1)) %>%
  mutate(vaccine_efficacy = paste0(median, " (",lower,"-",upper,")")) 

tb_dose_2 <- table1 %>%
  filter(dose==2) %>%
  filter(t==1+90 | t==1+180 ) %>%
  mutate(period = case_when(t==91 ~ "90d pd2",
                            t==181 ~ "180d pd2") )

tb_dose_3 <- table1 %>%
  filter(dose==3) %>% 
  filter(t==1+30 | t==1+60  | t==1+90 | t==1+120 | t==1+150 | t==1+180 | t==1+365) %>%
  mutate(period = case_when(t==31 ~ "30d pb",
                            t==61 ~ "60d pb",
                            t==91 ~ "90d pb",
                            t==121 ~ "120d pb",
                            t==151 ~ "150d pb",
                            t==181 ~ "180d pb",
                            t==366 ~ "365d pb") )

tb <- rbind(tb_dose_2,tb_dose_3)
tb <- subset(tb,select=c(vfr,vaccine,group,vaccine_efficacy,period))


tb.inf <- tb %>%
  filter(group == "Efficacy") %>%
  pivot_wider(id_cols = c(vaccine, vfr), names_from = period, values_from = vaccine_efficacy) %>%
  arrange(vaccine, vfr) %>%
  rename(Vaccine = vaccine, "0 = delta, 1 = omicron" = vfr) %>% 
  flextable()

tb.sev <- tb %>%
  filter(group == "Severe_Efficacy") %>%
  pivot_wider(id_cols = c(vaccine, vfr), names_from = period, values_from = vaccine_efficacy) %>%
  arrange(vaccine, vfr) %>%
  rename(Vaccine = vaccine, "0 = delta, 1 = omicron" = vfr) %>% 
  flextable()

tb.death <- tb %>%
  filter(group == "Death_Efficacy") %>%
  pivot_wider(id_cols = c(vaccine, vfr), names_from = period, values_from = vaccine_efficacy) %>%
  arrange(vaccine, vfr) %>%
  rename(Vaccine = vaccine, "0 = delta, 1 = omicron" = vfr) %>% 
  flextable()


save_as_docx("Efficacy against infection" = tb.inf,
             "Efficacy against severe disease" = tb.sev,
             "Efficacy against death" = tb.death,
            path = "../Tables/TableS5.docx")

