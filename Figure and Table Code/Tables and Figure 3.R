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

source("R/vx_profile.R")


##################################
##### LOAD THE PARAMETERS 
##################################

load("../Model Fitting/Main/UKHSA_v6pn_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain.Rdata")  
#load("../Model Fitting/Main/UKHSA_v6pn_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain.Rdata") 


chain <- mcmc$output %>%
  filter(phase == "sampling") %>%
  subset(select=-c(chain,phase,iteration,logprior,loglikelihood))

## Calculate parameter estimates and bounds of transformed parameters from 10,000 MCMC samples

draws <- sample_chains(mcmc, 10000)

draws_transform <- draws %>%
  select(-sample, -AZ_ns_off ) %>%
  mutate( ab50 = 10^ni50, 
          ab50_s = 10^ns50, 
          ab50_d = 10^nd50, 
          d1_AZ = 10^(d2_AZ + d1_AZ),
          d1_PF = 10^(d2_PF + d1_PF),
          d1_MD = 10^(d2_MD + d1_MD),
          d2_PF = 10^(d2_PF),
          d2_AZ = 10^(d2_AZ),
          d2_MD = 10^(d2_MD),
          d3_PF = 10^(bst_PF),
          d3_AZ = 10^(bst_AZ),
          d3_MD = 10^(bst_MD),
          om_red = 10^(om_red)) %>%
  select(-ni50, -ns50, -nd50, -bst_PF, -bst_AZ, -bst_MD) 
  
posterior_median_transform <-draws_transform %>%
  summarise( 
    across(where(is.numeric), median)
  )%>%
  mutate(measure = "median")

posterior_upper <- draws_transform %>%
  summarise( 
    across(where(is.numeric), quantile, 0.975)
  )%>%
  mutate(measure = "upper")

posterior_lower <- draws_transform %>%
  summarise( 
    across(where(is.numeric), quantile, 0.025)
  ) %>%
  mutate(measure = "lower")

params_est <- posterior_median_transform %>%
  rbind(posterior_upper) %>%
  rbind(posterior_lower) %>%
  pivot_longer(cols = c(d1_AZ, d1_PF, d1_MD,d2_PF, d2_AZ,d2_MD, d3_AZ, d3_PF, d3_MD, 
                        om_red, ab50, ab50_s, ab50_d, k, hl_s, hl_l, period_s, period_l)) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(median = round(median,digits=3), lower = round(lower,digits=3), upper = round(upper,digits=3)) %>%
  flextable()


save_as_docx("Posterior Parameter Estimates" = params_est,
               path = "../Tables/Table1.docx")

## Calculate Table 2 and 3 estimates from full chain and bounds from 1,000 MCMC samples

posterior_median <- chain %>%
  summarise( 
    across(where(is.numeric), median)
  ) %>%
  mutate(sample=0)
draws <- sample_chains(mcmc,1000)

################################################################################
### TABLE 2 ####################################################################
################################################################################

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

df_AZ_PF <- df_raw %>%
  mutate(vaccine = "AZ_PF",
         d1 = d1_AZ,
         d2 = d2_AZ,
         d3 = bst_PF,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df_AZ_MD <- df_raw %>%
  mutate(vaccine = "AZ_MD",
         d1 = d1_AZ,
         d2 = d2_AZ,
         d3 = bst_MD,
         ab50 = ni50,
         ab50_s =  ns50, 
         ab50_d = nd50) %>% 
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

df_PF_MD <- df_raw %>%
  mutate(vaccine = "PF_MD",
         d1 = d1_PF,
         d2 = d2_PF,
         d3 = bst_MD,
         ab50 =  ni50,
         ab50_s =  ns50, 
         ab50_d =  nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df_MD_PF <- df_raw %>%
  mutate(vaccine = "MD_PF",
         d1 = d1_MD,
         d2 = d2_MD,
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

#df <- bind_rows(df_AZ_AZ, df_AZ_PF, df_AZ_MD, df_PF_PF, df_PF_MD, df_MD_PF, df_MD_MD) 
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
               path = "../Tables/Table2.docx")

save_as_docx("Efficacy against infection" = tb.inf,
             "Efficacy against severe disease" = tb.sev,
             "Efficacy against death" = tb.death,
             path = "../Tables/TableS2.docx")
  
################################################################################
### TABLE 3 ####################################################################
################################################################################


# Moderna bivalent vaccine against BA1/2

ratio <- 2372/1473  # ratio of day 29 NAT against BA1/2 for the bivalent vaccine compared to monovalent vaccine = 1.61
bivalent_scaling <- 1/ratio

df_raw_bivalent <- df_raw %>%
  mutate(om_red = om_red*bivalent_scaling)


df_MD_MD <- df_raw_bivalent %>%
  mutate(vaccine = "MD_MD",
         d1 = d1_MD,
         d2 = d2_MD,
         d3 = bst_MD,
         ab50 = d2_PF + ni50,
         ab50_s = d2_PF + ns50, 
         ab50_d = d2_PF + nd50) %>% 
  select(-c(d1_PF, d1_AZ, d1_MD, d2_PF, d2_AZ, d2_MD, bst_PF,bst_AZ, bst_MD,
            AZ_ns_off,ni50,ns50,nd50, period_l)) 

df <- df_MD_MD

r2 <- 
  # Create input options
  expand_grid(
    vfr = c(1),
    vaccine = c("MD_MD"),
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
r2a <- r2 %>%
  filter(sample==0) %>%
  mutate(median = Value)

# mcmc samples
r2b <- r2 %>%
  filter(sample>0) %>%
  filter(dose>1) %>%
  group_by(vfr,vaccine,group,dose,t) %>%
  summarise(upper = quantile(Value, 0.975),
            lower = quantile(Value, 0.025),
  ) 

summary_stats <- left_join(r2a,r2b)


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
             path = "../Tables/Table3.docx")

################################################################################
### DATA FOR FIGURE 3 ####################################################################
################################################################################

df_ob <- r1a %>%
  filter((vaccine=="MD_MD") & (vfr==1) & (dose==3)) %>%
  filter(t<=365) %>%
  mutate(VE_original_boost = median) %>%
  select(c(t,group,VE_original_boost))

df_nob <- r1a %>%
  filter((vaccine=="MD_MD") & (vfr==1) & (dose==3)) %>%
  mutate(VE_no_boost = median) %>%
  mutate(t = t - 365) %>%
  filter(t>=0) %>%
  select(c(t,group,VE_no_boost))

df_bb <- r2a %>%
  filter((vaccine=="MD_MD") & (vfr==1) & (dose==3)) %>%
  filter(t<=365) %>%
  mutate(VE_bivalent_boost = median) %>%
  select(c(t,group,VE_bivalent_boost))

combined <- left_join(df_nob, df_ob)
combined <- left_join(combined, df_bb)

combined <- combined %>%
  filter(group != "Titre") %>%
  mutate(outcome = if_else(group == "Efficacy", "mild disease", group),
         outcome = if_else(group == "Severe_Efficacy", "hospitalisation", outcome),
         outcome = if_else(group == "Death_Efficacy", "death", outcome)) %>%
  arrange(outcome,t) %>%
  select (-c(group)) %>%
  mutate(no_boost_to_original = VE_original_boost - VE_no_boost,
         original_to_bivalent = VE_bivalent_boost - VE_original_boost,
         no_boost_to_bivalent = VE_bivalent_boost - VE_no_boost) %>%
  mutate(outcome = factor(outcome, levels = c("mild disease", "hospitalisation", "death"))) 

#################################################################################################

# trajectories

df_trajectories <- combined %>%
  select(t, outcome, VE_original_boost, VE_bivalent_boost, VE_no_boost) %>%
  pivot_longer(cols = c(VE_original_boost, VE_bivalent_boost, VE_no_boost)) %>%
  mutate(name = factor(name, levels = c("VE_no_boost", "VE_original_boost", "VE_bivalent_boost"), labels = c("3 doses only", "4th dose monovalent", "4th dose bivalent"))) 

p2 <- ggplot(df_trajectories, aes(x = t, y = value, col = name)) +
  geom_line() +
  facet_wrap(~outcome) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  labs(x = "time following dose 4 (days)", y = "effectiveness (%)", col = "")

p2

#################################################################################################

# barplots


day_outcomes <- combined %>%
  select(t, outcome, no_boost_to_original, original_to_bivalent) %>%
  pivot_longer(col = c(no_boost_to_original, original_to_bivalent)) %>%
  group_by(t, outcome) %>%
  mutate(effect_proportion = value / sum(value)) %>%
  ungroup() %>%
  mutate(name = factor(name, levels = c("original_to_bivalent", "no_boost_to_original"), labels = c("bivalent vs original","any 4th dose vs no 4th dose")))

overall <- combined %>%
  group_by(outcome) %>%
  summarise(no_boost_to_original = sum(no_boost_to_original),
            original_to_bivalent = sum(original_to_bivalent)) %>%
  pivot_longer(col = c(no_boost_to_original, original_to_bivalent)) %>%
  group_by(outcome) %>%
  mutate(effect_proportion = value / sum(value)) %>%
  ungroup() %>%
  mutate(t = "overall")

df_out <- rbind(day_outcomes, overall) %>%
  mutate(name = factor(name, levels = c("original_to_bivalent", "no_boost_to_original"), labels = c("bivalent vs original","any 4th dose vs no 4th dose")))

p1 <- ggplot(data = day_outcomes, aes(x = t, y = effect_proportion * 100, fill = name)) +
  geom_area(alpha = 0.8) +
  facet_wrap(~ outcome, nrow = 1) +
  labs(x = "time after boost (days)", y = "effect proportion (%)", fill = "") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_manual(values = c("#7570b3","#d95f02" ))
p1

################################################

library(patchwork)

layout <- "
A
B
"
combined <- p2 + p1 +
  plot_annotation(tag_levels = "A") + 
  plot_layout(design = layout)
combined

ggsave("../Figures/Figure 3.png", combined, height = 6, width = 8)



