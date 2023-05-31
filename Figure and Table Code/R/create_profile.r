## function to create titre and efficacy profiles

draw <-
  function(mu_ab_d1 = mu_ab_d1,
           mu_ab_d2 = mu_ab_d2,
           mu_ab_d3,
           std10 = std10,
           ab_50 = ab_50,
           ab_50_severe = ab_50_severe,
           ab_50_death = ab_50_death,
           dr_s = dr_s,
           dr_l = dr_l,
           period_s = period_s,
           t = t,
           t_d2 = t_d2,
           t_d3,
           k = k
  ){
  
  # decay rates over time
  ## simple biphasic decay on log10 scale
  denom=log10(exp(dr_l*period_s)+exp(dr_s*period_s))
  cum_dr_vec=log10(exp(dr_s*t+dr_l*period_s)+exp(dr_l*t+dr_s*period_s))-denom
 
  # initiate titre vector
  nt <- rep(0, length(t))
  nt[1] <- mu_ab_d1
  nt[2:t_d2] <- nt[1] + cum_dr_vec[2:t_d2]
  nt[t_d2+1] <- mu_ab_d2
  nt[(t_d2+2):(t_d2+t_d3)] <- nt[t_d2+1] + cum_dr_vec[2:t_d3]
  nt[t_d2+t_d3+1] <- mu_ab_d3
  nt[(t_d2+t_d3+2):length(t)] <-  nt[t_d2+t_d3+1] + cum_dr_vec[2:(length(t)-(t_d2+t_d3))]
 
  
  # relate titre to efficacy over time - using log-10 parameters
  ef_infection <- 1 / (1 + exp(-k * (nt - ab_50)))
  ef_severe <- 1 / (1 + exp(-k * (nt - ab_50_severe)))
  ef_death <- 1 / (1 + exp(-k * (nt - ab_50_death)))
  
  # return titres on normal scale
  nt <- 10^(nt)
  
  # output list containing titre and vaccine efficacy
  ret <- list(titre = nt, VE = ef_infection, VE_severe = ef_severe, VE_death = ef_death)
  return(ret)
  }
