library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(emdbook)
library(drjacoby)

setwd(dirname(getActiveDocumentContext()$path)) 


source("R/vx_profile.R")
source("R/vx_profile_1.R")


##################################
##### SUMMARY FITS
##################################

# main model

load("../Model Fitting/Main/UKHSA_v6_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain.Rdata")  
draws <- sample_chains(mcmc, 5000)
draws_transform <- draws %>%
  select(-sample, -AZ_ns_off ) %>%
  mutate( ab50 = 10^(d2_PF + ni50),
          ab50_s = 10^(d2_PF + ns50),
          ab50_d = 10^(d2_PF + nd50), 
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

main_posterior_median_transform <-draws_transform %>%
  summarise( 
    across(where(is.numeric), median)
  )%>%
  mutate(measure = "median")

chain_main <- mcmc$output %>%
  filter(phase=="sampling") %>%
  mutate(fit = logprior + loglikelihood )

chain_main %>%
  summarise(fit = mean(fit), LnL = mean(loglikelihood), LnL_SE = sd(loglikelihood))

load("../Model Fitting/Sensitivity - Additive Boost/UKHSA_v6_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=TRUE_AltSev=FALSE_mcmc_chain.Rdata")  

chain_addboost <- mcmc$output %>%
  filter(phase=="sampling") %>%
  mutate(fit = logprior + loglikelihood )

chain_addboost %>%
  summarise(fit = mean(fit), LnL = mean(loglikelihood), LnL_SE = sd(loglikelihood))

load("../Model Fitting/Sensitivity - Alternative Severity/UKHSA_v6_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=TRUE_mcmc_chain.Rdata")  

chain_altsev <- mcmc$output %>%
  filter(phase=="sampling") %>%
  mutate(fit = logprior + loglikelihood )

chain_altsev %>%
  summarise(fit = mean(fit), LnL = mean(loglikelihood), LnL_SE = sd(loglikelihood))

##################################
##### ALTERNATIVE SEVERITY
##################################

load("../Model Fitting/Sensitivity - Alternative Severity/UKHSA_v6_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=TRUE_mcmc_chain.Rdata")  
draws <- sample_chains(mcmc, 5000)

draws_transform <- draws %>%
  select(-sample, -AZ_ns_off ) %>%
  mutate( ab50 = 10^(d2_PF + ni50),
          ab50_s = ns50,
          ab50_d = nd50,
      #    ab50_s = 10^(d2_PF + ns50), # same function as infection - modify relative risk later
      #    ab50_d = 10^(d2_PF + nd50), # same function as infection - modify relative risk later
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

altsev_posterior_median_transform <- posterior_median_transform

n<-lseq(0.01,100,length.out=100)

ni50v        <- draws_transform$ab50
ns50v        <- draws_transform$ab50_s
nd50v        <- draws_transform$ab50_d
kv           <- draws_transform$k
om_redv       <- draws_transform$om_red

ru <-NULL

for (i in (1:length(kv))){
  
  ni50         <- ni50v[i] # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ns50         <- ns50v[i]
  nd50         <- nd50v[i]
  k            <- kv[i] # shape parameter of efficacy curve
  
  ef_infection <- 1 / (1 + exp(-k * (log10(n) - log10(ni50))))
  ef_severe <- 1 - (1-ef_infection)*exp(ns50)
  ef_death <- 1 - (1-ef_infection)*exp(nd50)
  
#  ef_severe    <- 1 / (1 + exp(-k * (log10(n) - log10(ns50))))
#  ef_death    <- 1 / (1 + exp(-k * (log10(n) - log10(nd50))))
  
  r <- data.frame(n=n, run=rep(i,length(n)), ef_infection=ef_infection, ef_severe=ef_severe, ef_death=ef_death)  
  ru<-rbind(ru,r)
  
}

ru <- ru %>%
  pivot_longer(cols = c("ef_infection", "ef_severe", "ef_death"), names_to = "type") %>%
  mutate(type = factor(type, levels = c("ef_infection", "ef_severe", "ef_death")))

ru_summary <- ru %>%
  group_by(type, n) %>%
  summarise(median = median(value),
            upper = quantile(value, 0.975),
            lower = quantile(value, 0.025),)

ni50      <- posterior_median_transform$ab50 
ns50      <- posterior_median_transform$ab50_s 
nd50      <- posterior_median_transform$ab50_d 
k         <- posterior_median_transform$k


ef_infection <- 1 / (1 + exp(-k * (log10(n) - log10(ni50))))
ef_severe <- 1 - (1-ef_infection)*exp(ns50)
ef_death <- 1 - (1-ef_infection)*exp(nd50)

#ef_severe <- 1 / (1 + exp(-k * (log10(n) - log10(ns50))))
#ef_death <- 1 / (1 + exp(-k * (log10(n) - log10(nd50))))

r1 <- data.frame(n,ef_infection,ef_severe, ef_death, model = "Alternative Severity")


# include main model curves

ni50_2      <- main_posterior_median_transform$ab50 
ns50_2      <- main_posterior_median_transform$ab50_s 
nd50_2      <- main_posterior_median_transform$ab50_d 
k_2      <- main_posterior_median_transform$k


ef_infection <- 1 / (1 + exp(-k_2 * (log10(n) - log10(ni50_2))))
ef_severe <- 1 / (1 + exp(-k_2 * (log10(n) - log10(ns50_2))))
ef_death <- 1 / (1 + exp(-k_2 * (log10(n) - log10(nd50_2))))

r2 <- data.frame(n,ef_infection,ef_severe, ef_death , model = "Main")

df <- rbind(r1, r2)

####

g1 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_infection"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = df, aes(x = n, y = 100*ef_infection, linetype = model), size = 1,col = "darkblue") +
  labs(x = "immunity level", y = "vaccine effectiveness mild (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1)) +
  scale_linetype_manual(values = c("dotted", "solid"))

g1

g2 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_severe"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = df, aes(x = n, y = 100*ef_severe, linetype = model), size = 1,col = "darkblue") +
  labs(x = "immunity level", y = "vaccine effectiveness hospitalisation (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))+
  scale_linetype_manual(values = c("dotted", "solid"))


g2

g3 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_death"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = df, aes(x = n, y = 100*ef_death, linetype = model), size = 1,col = "darkblue") +
  labs(x = "immunity level", y = "vaccine effectiveness death (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))+
  scale_linetype_manual(values = c("dotted", "solid"))

g3




### generate VE estimates for Alternative Severity and Main model over 1 year following boosting 

df_raw <- main_posterior_median_transform

df_AZ_AZ <- df_raw %>%
  mutate(vaccine = "AZ_AZ",
         d2 = log10(d2_AZ),
         d1 = log10(d1_AZ) - d2,
         d3 = log10(d3_AZ),
         ab50 = log10(ab50),
         ab50_s = log10(ab50_s),
         ab50_d = log10(ab50_d),
         om_red = log10(om_red),
         period_s_2 = period_l,
         sample = 1) %>% 
  select(-c(measure, d1_AZ, d1_PF, d1_MD, d2_AZ, d2_PF, d2_MD, d3_AZ, d3_PF, d3_MD, period_l)) 

df_PF_PF <- df_raw %>%
  mutate(vaccine = "PF_PF",
         d2 = log10(d2_PF),
         d1 = log10(d1_PF) - d2,
         d3 = log10(d3_PF),
         ab50 = log10(ab50),
         ab50_s = log10(ab50_s),
         ab50_d = log10(ab50_d),
         om_red = log10(om_red),
         period_s_2 = period_l,
         sample = 1) %>% 
  select(-c(measure, d1_AZ, d1_PF, d1_MD, d2_AZ, d2_PF, d2_MD, d3_AZ, d3_PF, d3_MD, period_l)) 

df_MD_MD <- df_raw %>%
  mutate(vaccine = "MD_MD",
         d2 = log10(d2_MD),
         d1 = log10(d1_MD) - d2,
         d3 = log10(d3_MD),
         ab50 = log10(ab50),
         ab50_s = log10(ab50_s),
         ab50_d = log10(ab50_d),
         om_red = log10(om_red),
         period_s_2 = period_l,
         sample = 1) %>% 
  select(-c(measure, d1_AZ, d1_PF, d1_MD, d2_AZ, d2_PF, d2_MD, d3_AZ, d3_PF, d3_MD, period_l)) 

df_main <- bind_rows(df_AZ_AZ, df_PF_PF, df_MD_MD) 

df_raw <- altsev_posterior_median_transform

df_AZ_AZ <- df_raw %>%
  mutate(vaccine = "AZ_AZ",
         d2 = log10(d2_AZ),
         d1 = log10(d1_AZ) - d2,
         d3 = log10(d3_AZ),
         ab50 = log10(ab50),
         om_red = log10(om_red),
         period_s_2 = period_l,
         sample = 1) %>% 
  select(-c(measure, d1_AZ, d1_PF, d1_MD, d2_AZ, d2_PF, d2_MD, d3_AZ, d3_PF, d3_MD, period_l)) 

df_PF_PF <- df_raw %>%
  mutate(vaccine = "PF_PF",
         d2 = log10(d2_PF),
         d1 = log10(d1_PF) - d2,
         d3 = log10(d3_PF),
         ab50 = log10(ab50),
         om_red = log10(om_red),
         period_s_2 = period_l,
         sample = 1) %>% 
  select(-c(measure, d1_AZ, d1_PF, d1_MD, d2_AZ, d2_PF, d2_MD, d3_AZ, d3_PF, d3_MD, period_l)) 

df_MD_MD <- df_raw %>%
  mutate(vaccine = "MD_MD",
         d2 = log10(d2_MD),
         d1 = log10(d1_MD) - d2,
         d3 = log10(d3_MD),
         ab50 = log10(ab50),
          om_red = log10(om_red),
         period_s_2 = period_l,
         sample = 1) %>% 
  select(-c(measure, d1_AZ, d1_PF, d1_MD, d2_AZ, d2_PF, d2_MD, d3_AZ, d3_PF, d3_MD, period_l)) 

df_altsev <- bind_rows(df_AZ_AZ, df_PF_PF, df_MD_MD) 

r1_main <- 
  # Create input options
  expand_grid(
    vfr = 1,
    vaccine = c("AZ_AZ", "PF_PF", "MD_MD"),
    sample = 1) %>%
  # Join with vaccine parameters
  left_join(df_main, by = c("vaccine", "sample")) %>%
  # Apply vaccine profile function to each row
  mutate(profile = pmap(., vx_profile)) %>%
  # Format
  unnest(cols = c(profile)) %>%
  pivot_longer(cols = c(Titre_d1, Titre_d2, Titre_d3, Efficacy_d1, Efficacy_d2, Efficacy_d3, Severe_Efficacy_d1,Severe_Efficacy_d2, Severe_Efficacy_d3, Death_Efficacy_d1, Death_Efficacy_d2,Death_Efficacy_d3), names_to = "group", values_to = "Value") %>% 
  mutate(dose = case_when(group == "Titre_d1" | group == "Efficacy_d1" | group == "Severe_Efficacy_d1" | group == "Death_Efficacy_d1"~ 1,
                          group == "Titre_d2" | group == "Efficacy_d2" | group == "Severe_Efficacy_d2" | group == "Death_Efficacy_d2" ~ 2,
                          group == "Titre_d3" | group == "Efficacy_d3" | group == "Severe_Efficacy_d3" | group == "Death_Efficacy_d3"~ 3 )) %>%
  mutate(group = substr(group, 1, nchar(group) - 3))

r2_main <- r1_main %>%
  filter(dose==3 & vaccine=="MD_MD" & (t<=365) & group != "Titre") %>%
  mutate(model="Main")

r1_altsev <- 
  # Create input options
  expand_grid(
    vfr = 1,
    vaccine = c("AZ_AZ", "PF_PF", "MD_MD"),
    sample = 1) %>%
  # Join with vaccine parameters
  left_join(df_altsev, by = c("vaccine", "sample")) %>%
  # Apply vaccine profile function to each row
  mutate(profile = pmap(., vx_profile_altsev)) %>%
  # Format
  unnest(cols = c(profile)) %>%
  pivot_longer(cols = c(Titre_d1, Titre_d2, Titre_d3, Efficacy_d1, Efficacy_d2, Efficacy_d3, Severe_Efficacy_d1,Severe_Efficacy_d2, Severe_Efficacy_d3, Death_Efficacy_d1, Death_Efficacy_d2,Death_Efficacy_d3), names_to = "group", values_to = "Value") %>% 
  mutate(dose = case_when(group == "Titre_d1" | group == "Efficacy_d1" | group == "Severe_Efficacy_d1" | group == "Death_Efficacy_d1"~ 1,
                          group == "Titre_d2" | group == "Efficacy_d2" | group == "Severe_Efficacy_d2" | group == "Death_Efficacy_d2" ~ 2,
                          group == "Titre_d3" | group == "Efficacy_d3" | group == "Severe_Efficacy_d3" | group == "Death_Efficacy_d3"~ 3 )) %>%
  mutate(group = substr(group, 1, nchar(group) - 3))


r2_altsev <- r1_altsev %>%
  filter(dose==3 & vaccine=="MD_MD" & (t<=365) & group != "Titre") %>%
  mutate(model="Alternative Severity")

r2 <-bind_rows(r2_main, r2_altsev) %>% 
  mutate(outcome = if_else(group == "Efficacy", "mild disease", group),
         outcome = if_else(group == "Severe_Efficacy", "hospitalisation", outcome),
         outcome = if_else(group == "Death_Efficacy", "death", outcome)) %>%
  mutate(outcome = factor(outcome, levels = c("mild disease", "hospitalisation", "death")))

p1 <- ggplot(r2,aes(x = t, y = Value,  linetype = model) ) +
  geom_line(col = "darkblue", size=1) +
  facet_wrap(~outcome) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
 #       panel.grid = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none") +
 labs(x = "time after boost (days)", y = "vaccine effectiveness (%)", col = "") +
  scale_linetype_manual(values = c("dotted", "solid"))
p1

layout <- "
ABC
DDD "
  
combined <- g1 + g2  + g3 + p1 + plot_annotation(tag_levels = "A") +  plot_layout(design = layout, guides = "collect") & theme(legend.position = 'bottom')

combined

ggsave("../Figures/Figure S2 Alternative Severity.png",combined, height = 10, width = 10)  

##################################
##### ADDITIVE BOOST
##################################

load("../Model Fitting/Sensitivity - Additive Boost/UKHSA_v6_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=TRUE_AltSev=FALSE_mcmc_chain.Rdata")  
draws <- sample_chains(mcmc, 5000)

draws_transform <- draws %>%
  select(-sample, -AZ_ns_off ) %>%
  mutate( ab50 = 10^(d2_PF + ni50),
          ab50_s = 10^(d2_PF + ns50), 
          ab50_d = 10^(d2_PF + nd50),
          d1_AZ = 10^(d2_AZ + d1_AZ),
          d1_PF = 10^(d2_PF + d1_PF),
          d1_MD = 10^(d2_MD + d1_MD),
          d3_PF = 10^(d2_PF + bst_PF),
          d3_AZ = 10^(d2_AZ + bst_AZ),
          d3_MD = 10^(d2_MD + bst_MD),
          d2_PF = 10^(d2_PF),
          d2_AZ = 10^(d2_AZ),
          d2_MD = 10^(d2_MD),
          om_red = 10^(om_red)) %>%
  select(-ni50, -ns50, -nd50, -bst_PF, -bst_AZ, -bst_MD) 

posterior_median_transform <-draws_transform %>%
  summarise( 
    across(where(is.numeric), median)
  )%>%
  mutate(measure = "median")


n<-lseq(0.01,100,length.out=100)

ni50v        <- draws_transform$ab50
ns50v        <- draws_transform$ab50_s
nd50v        <- draws_transform$ab50_d
kv           <- draws_transform$k
om_redv       <- draws_transform$om_red

ru <-NULL

for (i in (1:length(kv))){
  
  ni50         <- ni50v[i] # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ns50         <- ns50v[i]
  nd50         <- nd50v[i]
  k            <- kv[i] # shape parameter of efficacy curve
  
  ef_infection <- 1 / (1 + exp(-k * (log10(n) - log10(ni50))))
  ef_severe    <- 1 / (1 + exp(-k * (log10(n) - log10(ns50))))
  ef_death    <- 1 / (1 + exp(-k * (log10(n) - log10(nd50))))
  
  r <- data.frame(n=n, run=rep(i,length(n)), ef_infection=ef_infection, ef_severe=ef_severe, ef_death=ef_death)  
  ru<-rbind(ru,r)
  
}

ru <- ru %>%
  pivot_longer(cols = c("ef_infection", "ef_severe", "ef_death"), names_to = "type") %>%
  mutate(type = factor(type, levels = c("ef_infection", "ef_severe", "ef_death")))

ru_summary <- ru %>%
  group_by(type, n) %>%
  summarise(median = median(value),
            upper = quantile(value, 0.975),
            lower = quantile(value, 0.025),)

ni50      <- posterior_median_transform$ab50 
ns50      <- posterior_median_transform$ab50_s 
nd50      <- posterior_median_transform$ab50_d 
k         <- posterior_median_transform$k


ef_infection <- 1 / (1 + exp(-k * (log10(n) - log10(ni50))))
ef_severe <- 1 / (1 + exp(-k * (log10(n) - log10(ns50))))
ef_death <- 1 / (1 + exp(-k * (log10(n) - log10(nd50))))

r1 <- data.frame(n,ef_infection,ef_severe, ef_death, model = "Alternative Boost")


# include main model curves

ni50_2      <- main_posterior_median_transform$ab50 
ns50_2      <- main_posterior_median_transform$ab50_s 
nd50_2      <- main_posterior_median_transform$ab50_d 
k_2      <- main_posterior_median_transform$k


ef_infection <- 1 / (1 + exp(-k_2 * (log10(n) - log10(ni50_2))))
ef_severe <- 1 / (1 + exp(-k_2 * (log10(n) - log10(ns50_2))))
ef_death <- 1 / (1 + exp(-k_2 * (log10(n) - log10(nd50_2))))

r2 <- data.frame(n,ef_infection,ef_severe, ef_death , model = "Main")
df <- rbind(r1, r2)

####

g4 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_infection"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = df, aes(x = n, y = 100*ef_infection, linetype = model), size = 1,col = "darkblue") +
  labs(x = "immunity level", y = "vaccine effectiveness mild (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1)) +
  scale_linetype_manual(values = c("dotted", "solid"))

g4

g5 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_severe"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = df, aes(x = n, y = 100*ef_severe, linetype = model), size = 1,col = "darkblue") +
  labs(x = "immunity level", y = "vaccine effectiveness hospitalisation (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))+
  scale_linetype_manual(values = c("dotted", "solid"))


g5

g6 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_death"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = df, aes(x = n, y = 100*ef_death, linetype = model), size = 1,col = "darkblue") +
  labs(x = "immunity level", y = "vaccine effectiveness death (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))+
  scale_linetype_manual(values = c("dotted", "solid"))

g6


combined <-  g4 + g5 + g6 + plot_annotation(tag_levels = "A") +  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

combined

ggsave("../Figures/Figure S3 Additive Boost.png",combined, height = 6, width = 12)  








