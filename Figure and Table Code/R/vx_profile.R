vx_profile <- function(vfr,
                       vaccine,
                       sample,
                       d1,
                       om_red,
                       ab50,
                       ab50_s,
                       ab50_d,
                       k,
                       hl_s,
                       hl_l,
                       period_s,
                       period_s_2,
                       d2,
                       d3
                       ) {
  
  mu_ab_d1 <- d2 + d1 - vfr*om_red
  mu_ab_d2 <- d2 - vfr*om_red
  mu_ab_d3 <- d3 - vfr*om_red
  
  lg10 <- log(10)
  
  t <- 0:731  #vaccinated on day 0
  dr_s <- -log(2)/hl_s  # Corresponding decay rate in days for half life above
  dr_l <- -log(2)/hl_l 
  period_s <- period_s
  denom=log10(exp(dr_l*period_s)+exp(dr_s*period_s))
  cum_dr_vec=log10(exp(dr_s*t+dr_l*period_s)+exp(dr_l*t+dr_s*period_s))-denom
  
  n1 = mu_ab_d1+cum_dr_vec
  
  efficacy_dose1_infection <- 1/(1+exp(-k*(n1 -ab50)))
  efficacy_dose1_severe <- 1/(1+exp(-k*(n1 -ab50_s)))
  efficacy_dose1_death <- 1/(1+exp(-k*(n1 -ab50_d)))
  
  n2=mu_ab_d2 + cum_dr_vec
 
  efficacy_dose2_infection <- 1/(1+exp(-k*(n2 -ab50)))
  efficacy_dose2_severe <- 1/(1+exp(-k*(n2 -ab50_s)))
  efficacy_dose2_death <- 1/(1+exp(-k*(n2 -ab50_d)))
  
  denom=log10(exp(dr_l*period_s)+exp(dr_s*period_s))
  cum_dr_vec=log10(exp(dr_s*t+dr_l*period_s)+exp(dr_l*t+dr_s*period_s))-denom
  
  n3=mu_ab_d3 + cum_dr_vec
  
  efficacy_dose3_infection <- 1/(1+exp(-k*(n3 -ab50)))
  efficacy_dose3_severe <- 1/(1+exp(-k*(n3 -ab50_s)))
  efficacy_dose3_death <- 1/(1+exp(-k*(n3 -ab50_d)))
  
  sub <- data.frame(t=t, mu_ab_d1=10^(mu_ab_d1),
                    Titre_d1=10^n1, Efficacy_d1=efficacy_dose1_infection*100, Severe_Efficacy_d1=efficacy_dose1_severe*100,Death_Efficacy_d1=efficacy_dose1_death*100,
                    Titre_d2=10^n2, Efficacy_d2=efficacy_dose2_infection*100, Severe_Efficacy_d2=efficacy_dose2_severe*100,Death_Efficacy_d2=efficacy_dose2_death*100,
                    Titre_d3=10^n3, Efficacy_d3=efficacy_dose3_infection*100, Severe_Efficacy_d3=efficacy_dose3_severe*100,Death_Efficacy_d3=efficacy_dose3_death*100)
  return(sub)
}