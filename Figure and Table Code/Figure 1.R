library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(emdbook)
library(drjacoby)

setwd(dirname(getActiveDocumentContext()$path)) 

##################################
##### LOAD THE PARAMETERS 
##################################

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

r1 <- data.frame(n,ef_infection,ef_severe, ef_death, model = "refitted")


# include original Khoury et al Nat Med curves
ni50_2 <- 0.2
ns50_2 <- 0.03
k_2 <- 2.94
vfr <- 3.9


ef_infection <- 1 / (1 + exp(-k_2 * (log10(n/vfr) - log10(ni50_2))))
ef_severe <- 1 / (1 + exp(-k_2 * (log10(n/vfr) - log10(ns50_2))))

r2 <- data.frame(n,ef_infection,ef_severe, ef_death = NA, model = "Khoury et al")

df <- rbind(r1, r2)

####

g1 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_infection"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = df, aes(x = n, y = 100*ef_infection, linetype = model), size = 1,col = "darkblue") +
  labs(x = "immunity level", y = "vaccine effectiveness mild (%)", col = "model (Delta VFR)", linetype = "model (Delta VFR)") +
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
  labs(x = "immunity level", y = "vaccine effectiveness hospitalisation (%)", col = "model (Delta VFR)", linetype = "model (Delta VFR)") +
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
  labs(x = "immunity level", y = "vaccine effectiveness death (%)", col = "model (Delta VFR)", linetype = "model (Delta VFR)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))+
  scale_linetype_manual(values = c("dotted", "solid"))

g3

om_redv       <- draws_transform$om_red

ru <-NULL
om_scalingv <- om_redv



for (i in (1:length(kv))){
  
  ni50         <- ni50v[i] # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ns50         <- ns50v[i]
  nd50         <- nd50v[i]
  k            <- kv[i] # shape parameter of efficacy curve
  
  ef_infection <- 1 / (1 + exp(-k * (log10(n/om_scalingv[i]) - log10(ni50))))
  ef_severe    <- 1 / (1 + exp(-k * (log10(n/om_scalingv[i]) - log10(ns50))))
  ef_death    <- 1 / (1 + exp(-k * (log10(n/om_scalingv[i]) - log10(nd50))))
  
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

om_scaling <- posterior_median_transform$om_red
ni50      <- posterior_median_transform$ab50 
ns50      <- posterior_median_transform$ab50_s 
nd50      <- posterior_median_transform$ab50_d 
k         <- posterior_median_transform$k



ef_infection <- 1 / (1 + exp(-k * (log10(n/om_scaling) - log10(ni50))))
ef_severe <- 1 / (1 + exp(-k * (log10(n/om_scaling) - log10(ns50))))
ef_death <- 1 / (1 + exp(-k * (log10(n/om_scaling) - log10(nd50))))

r1 <- data.frame(n,ef_infection,ef_severe, ef_death, model = "refitted")


# include original Khoury et al Nat Med curves
ni50_2 <- 0.2
ns50_2 <- 0.03
k_2 <- 2.94
vfr <- 9.7

ef_infection <- 1 / (1 + exp(-k_2 * (log10(n/vfr) - log10(ni50_2))))
ef_severe <- 1 / (1 + exp(-k_2 * (log10(n/vfr) - log10(ns50_2))))

r2 <- data.frame(n,ef_infection,ef_severe, ef_death = NA, model = "Khoury et al")

df <- rbind(r1, r2)

####

g4 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_infection"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkred") +
  geom_line(data = df, aes(x = n, y = 100*ef_infection, linetype = model), size = 1,col = "darkred") +
  labs(x = "immunity level", y = "vaccine effectiveness mild (%)", col = "model (Delta VFR)", linetype = "model (Omicron VFR)") +
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
  geom_ribbon(data = filter(ru_summary, type == "ef_severe"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkred") +
  geom_line(data = df, aes(x = n, y = 100*ef_severe, linetype = model), size = 1,col = "darkred") +
  labs(x = "immunity level", y = "vaccine effectiveness hospitalisation (%)", col = "model (Delta VFR)", linetype = "model (Omicron VFR)") +
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
  geom_ribbon(data = filter(ru_summary, type == "ef_death"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkred") +
  geom_line(data = df, aes(x = n, y = 100*ef_death, linetype = model), size = 1,col = "darkred") +
  labs(x = "immunity level", y = "vaccine effectiveness death (%)", col = "model (Delta VFR)", linetype = "model (Omicron VFR)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))+
  scale_linetype_manual(values = c("dotted", "solid"))

g6




combined <- g1 + g2  + g3 + g4 + g5 + g6 + plot_annotation(tag_levels = "A") +  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

combined

ggsave("../Figures/Figure1.png",combined, height = 7, width = 10)  
