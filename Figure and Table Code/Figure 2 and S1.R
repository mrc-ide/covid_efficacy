library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(drjacoby)

setwd(dirname(getActiveDocumentContext()$path)) 

source("R/create_profile.R")

##################################
##### LOAD THE PARAMETERS AND DATA
##################################

## main fits
load("../Model Fitting/Main/UKHSA_v6pn_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain.Rdata")  

all_data <- read.csv("../Model Fitting/Data/UKHSA_VE_Jun22_65+.csv")

## fits using 18-64 age group

#load("../Model Fitting/Main/UKHSA_v6pn_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain.Rdata") 
#all_data <- read.csv("../Model Fitting/Data/UKHSA_VE_Jun22_18_64.csv")



chain <- sample_chains(mcmc, 10000)

posterior_median <- chain %>%
  summarise( 
  across(where(is.numeric), median)
  )

name <- c("AZ-AZ","AZ-PF","AZ-MD",
          "PF-PF","PF-MD",
          "MD-PF","MD-MD")

d1_name <- c("AZ", "AZ", "AZ",
             "PF", "PF", 
             "MD", "MD")


log10_d2_PF <- posterior_median$d2_PF
log10_d2_AZ <- posterior_median$d2_AZ
log10_d2_MD <- posterior_median$d2_MD


d1_AZ       <- log10_d2_AZ + posterior_median$d1_AZ
d1_PF       <- log10_d2_PF + posterior_median$d1_PF
d1_MD       <- log10_d2_MD + posterior_median$d1_MD

d2_AZ    <- posterior_median$d2_AZ
d2_PF    <- posterior_median$d2_PF
d2_MD    <- posterior_median$d2_MD

d3_AZ    <- posterior_median$bst_AZ
d3_PF    <- posterior_median$bst_PF
d3_MD    <- posterior_median$bst_MD

ab_50       <- posterior_median$ni50 
ab_50_severe <- posterior_median$ns50
ab_50_death  <- posterior_median$nd50

ab_50 <-rep(ab_50,7) 
ab_50_severe <-rep(ab_50_severe,7)
ab_50_death <-rep(ab_50_death,7)

k           <- posterior_median$k
hl_s        <- posterior_median$hl_s
hl_l        <- posterior_median$hl_l
period_s    <- posterior_median$period_s

# fixed parameter
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 scale

om_red <- posterior_median$om_red
vfr <- c(0,om_red)

mu_ab_d1 <- c(d1_AZ, d1_AZ, d1_AZ, d1_PF, d1_PF, d1_MD, d1_MD)
mu_ab_d2 <- c(d2_AZ, d2_AZ, d2_AZ, d2_PF, d2_PF, d2_MD, d2_MD)
mu_ab_d3 <- c(d3_AZ, d3_PF, d3_MD, d3_PF, d3_MD, d3_PF, d3_MD)   


# transforms
dr_s <- -log(2)/hl_s  # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

# Timing of doses
max_t <- 365*2 # number of days to model
t <- 0:max_t
t_d2 <- 84 # timing of second dose relative to first
t_d3 <- 180 # timing of third dose relative to second dose

param_list <-
  data.frame(
    name,
    mu_ab_d1,
    mu_ab_d2,
    mu_ab_d3,
    t_d2,
    t_d3,
    dr_s,
    dr_l,
    period_s,
    ab_50,
    ab_50_severe,
    ab_50_death
  ) 
  
## manipulate data for plotting
## remove data points not used in fitting

all_data <- all_data %>%
            mutate(t_orig=floor((t_max+t_min)/2)) %>%
            filter(t_orig>=14) %>%
            mutate(t_fit= case_when( dose==1 ~ t_orig-21,
                                     dose>1 ~ t_orig-14)) %>%
            filter((dose==1 & t_fit<t_d2)|(dose==2 & t_fit<t_d3)|dose==3) %>%
            mutate(t= case_when( dose==1 ~ t_fit,
                                 dose==2 ~ t_d2+t_fit,
                                 dose==3 ~ t_d2+t_d3+t_fit )) %>%
            mutate(type=case_when( endpoint==1 ~ "Efficacy",
                                   endpoint==2 ~ "Efficacy_Severe",
                                   endpoint==3 ~ "Death"))

# initialise other parameters
r1_summary <- NULL
summary_stats <- NULL

##########################
###### OUTPUT OPTIONS
##########################

#####################################################################################
### loop through the vaccines to calculate the profiles of NAT and efficacy over time
#####################################################################################


for (m in 1:2){ #1 for delta, 2 for omicron
for (j in 1:7){
  
          mu_ab_d1 <- param_list$mu_ab_d1[j] - vfr[m]
          mu_ab_d2 <- param_list$mu_ab_d2[j] - vfr[m]
          mu_ab_d3 <- param_list$mu_ab_d3[j] - vfr[m]
          dr_s <- param_list$dr_s[j]
          dr_l <- param_list$dr_l[j]
          period_s <- param_list$period_s[j]
          t_d2 <- param_list$t_d2[j]
          t_d3 <- param_list$t_d3[j]
          t_d4 <- param_list$t_d4[j]
          ab_50 <- param_list$ab_50[j]
          ab_50_severe <- param_list$ab_50_severe[j]
          ab_50_death <- param_list$ab_50_death[j]

          out0 <-
            draw(
              mu_ab_d1 = mu_ab_d1,
              mu_ab_d2 = mu_ab_d2,
              mu_ab_d3 = mu_ab_d3,
              std10 = 0,
              ab_50 = ab_50,
              ab_50_severe = ab_50_severe,
              ab_50_death = ab_50_death,
              dr_s = dr_s,
              dr_l = dr_l,
              period_s = period_s,
              t = t,
              t_d2 = t_d2,
              t_d3 = t_d3,
              k = k
            )
          
          sub0 <-
            data.frame(
              t = t,
              run = 0,
              variant = m,
              Titre = out0$titre,
              Efficacy = out0$VE,
              Efficacy_Severe = out0$VE_severe,
              Death = out0$VE_death
            )
        
          sub0 <- sub0 %>%
              mutate(Efficacy = Efficacy * 100,
              Efficacy_Severe = Efficacy_Severe * 100,
              Death = Death * 100) %>%
              pivot_longer(cols = c("Titre", "Efficacy", "Efficacy_Severe", "Death"), names_to = "type") %>%
              mutate(type = factor(type, levels = c("Titre", "Efficacy", "Efficacy_Severe", "Death"))) 
          

  
  if(j==1){
    uk_data <- all_data %>%
    #filter(vaccine_num==1 & (booster_num==0 | booster_num==1) ) 
    filter((vaccine == d1_name[j] & dose %in% c(1,2) ) | (dose == 3 & booster_num == 1))
    
  }
  if(j==2){
    uk_data <- all_data %>%
      filter((vaccine == d1_name[j] & dose %in% c(1,2) ) | (dose == 3 & booster_num == 2))
    #filter(vaccine_num==1 & (booster_num==0 | booster_num ==2) ) 
  }
  if(j==3){
    uk_data <- all_data %>%
      filter((vaccine == d1_name[j] & dose %in% c(1,2) ) | (dose == 3 & booster_num == 3))
    
    #filter(vaccine_num==1 & (booster_num==0 | booster_num ==3) ) 
          }
  if(j==4){
    uk_data <- all_data %>%
      filter((vaccine == d1_name[j] & dose %in% c(1,2) ) | (dose == 3 & booster_num == 2))

  }
  if(j==5){
    uk_data <- all_data %>%
    #filter(vaccine_num==2 & (booster_num==0 | booster_num ==3) ) 
      filter((vaccine == d1_name[j] & dose %in% c(1,2) ) | (dose == 3 & booster_num == 3))
    
  }
  if(j==6){
    uk_data <- all_data %>%
    #filter(vaccine_num==2 & (booster_num==0 | booster_num ==2) ) 
    filter((vaccine == d1_name[j] & dose %in% c(1,2) ) | (dose == 3 & booster_num == 2))
    
  }
  if(j==7){
    uk_data <- all_data %>%
    #filter(vaccine_num==3 & (booster_num==0 | booster_num ==3) ) 
      filter((vaccine == d1_name[j] & dose %in% c(1,2) ) | (dose == 3 & booster_num == 3))
    
  }
          
  sub0 <- sub0 %>%
    left_join(uk_data,by=c("t","type","variant"))%>%
    # dont allow data CIs to be <0
    mutate(L95 = if_else(L95 < 0, 0, L95)) %>%
    mutate(j = j)

  summary_stats <-rbind(summary_stats,sub0)
}}

df <- summary_stats %>%
  mutate(type = factor(type, levels = c("Titre", "Efficacy", "Efficacy_Severe", "Death"), labels = c("IL", "mild disease", "hospitalisation", "death"))) %>%
  mutate(variant = factor(variant, levels = c(1,2), labels = c("delta", "omicron")))

data_to_output <- df %>%
  mutate(vaccines = case_when(j==1 ~ "AZ-AZ",
                       j==2 ~ "AZ-PF",
                       j==3 ~ "AZ-MD",
                       j==4 ~ "PF-PF",
                       j==5 ~ "PF-MD",
                       j==6 ~ "MD-PF",
                       j==7 ~ "MD-MD")
  ) %>%
  select(-Label, -vaccine, -vaccine_num, -booster_num, -dose, -endpoint, -t_min, -t_max, -control_I, -case_I, -VE, -N,
         - L95, -U95, -Age, -t_orig, -t_fit, -j)
saveRDS(data_to_output, "../Figures/vaccine_profiles.rds")

plots <- NULL

plot_efficacy_function <- function(d, i){
  plot_out <- ggplot(data = filter(d, type != "IL"), aes(x = t, y = value, col = factor(variant))) +
    geom_line() +
    geom_point(aes(x = t, y=VE), size = 1.5) +
    geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    lims(y = c(0,100)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    labs(x = "time (days)", y = "efficacy (%)", col = "variant")
  return(plot_out)
}

plot_nat_function <- function(d, i){
  plot_out <- ggplot(data = filter(d, type == "IL"), aes(x = t, y = value, col = factor(variant))) +
    geom_line() +
    geom_point(aes(x = t, y=VE), size = 1.5) +
    geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    scale_y_log10(limits = c(1e-3,1e1)) +
    labs(x = "time (days)", y = "IL", col = "variant", title = name[i])
  return(plot_out)
}

plot_efficacy_function2 <- function(d, i){
  plot_out <- ggplot(data = filter(d, type != "IL"), aes(x = t, y = value, col = factor(variant))) +
    geom_line() +
    geom_point(aes(x = t, y=VE), size = 1.5) +
    geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    lims(y = c(0,100)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.title.x=element_blank(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    labs(x = "time (days)", y = "effectiveness (%)", col = "variant")
  return(plot_out)
}

plot_nat_function2 <- function(d, i){
  plot_out <- ggplot(data = filter(d, type == "IL"), aes(x = t, y = value, col = factor(variant))) +
    geom_line() +
    geom_point(aes(x = t, y=VE), size = 1.5) +
    geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.title.x=element_blank(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    scale_y_log10(limits = c(1e-3,1e1)) +
    labs(x = "time (days)", y = "IL", col = "variant")
  return(plot_out)
}

for (i in 1:7){
  d1 <- filter(df, j==i )
  plots[[i*2]] <- plot_efficacy_function2(d1, i)
  plots[[i*2-1]] <- plot_nat_function2(d1, i)
}

layout <- "
ABCDE
FGHHH
IJKKK
LMNNN
OPQQQ
RSTTT
UVWWW
XYZZZ
XYZZZ
"
plot_spacer() + 
  wrap_elements(grid::textGrob('immunity level')) +
  wrap_elements(grid::textGrob('mild disease')) +
  wrap_elements(grid::textGrob('hospitalisation')) +
  wrap_elements(grid::textGrob('death')) +
  wrap_elements(grid::textGrob('AZ-AZ')) + plots[[1]] + plots[[2]] + 
  wrap_elements(grid::textGrob('AZ-PF'))  + plots[[3]] + plots[[4]] +
  wrap_elements(grid::textGrob('AZ-MD'))  + plots[[5]] + plots[[6]] +
  wrap_elements(grid::textGrob('PF-PF')) + plots[[7]] + plots[[8]] + 
  wrap_elements(grid::textGrob('PF-MD')) + plots[[9]] + plots[[10]] + 
  wrap_elements(grid::textGrob('MD-PF')) + plots[[11]] + plots[[12]] + 
  wrap_elements(grid::textGrob('MD-MD')) + plots[[13]] + plots[[14]]+
  plot_spacer() + grid::textGrob('time (days)') + grid::textGrob('time (days)') +
  plot_layout(guides = "collect", design = layout, widths = c(0.3,1,1,1,1), heights = c(0.2, 1, 1, 1, 1, 1,1,1,0.2)) & theme(legend.position = 'bottom')

ggsave("../Figures/Figure 2.png", height = 16, width = 10)
#ggsave("../Figures/Figure S1.png", height = 16, width = 10)







